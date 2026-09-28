#!/usr/bin/env python3
"""LLM design classification: project packets (JSONL) -> structured verdicts (JSONL). On re-run, studies already classified are skipped.

Usage:
  python3 03_classify_design_llm.py --packets data/design_packets/labeled_v1.jsonl --out results/design_llm/labeled_v1.pred.jsonl [--model claude-opus-5] [--effort medium] [--workers 4] [--dry-run] [--limit N]
Auth: ANTHROPIC_API_KEY environment variable or an `ant auth login` profile.
"""
import argparse, json, os, sys, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Literal
from pydantic import BaseModel, Field

SYSTEM = """You are auditing public RNA-seq datasets (GEO/SRA) for a study on covariate (sex/age) confounding in differential expression analysis.
For each dataset you receive structured metadata: title, abstract, tissue, sample-attribute keys with their value levels and counts, and heuristic flags.
Decide whether the dataset supports a clean BETWEEN-SUBJECT two-group differential expression comparison on bulk tissue or primary cells from separate individuals/animals, and if so, define the comparison from the sample attributes.

Rules:
- usable_design = true only if: bulk RNA-seq (not single-cell), tissue or freshly isolated/sorted primary cells (not cultured, not stimulated in vitro, not cell lines, not organoids/iPSC), and the two groups are different subjects/animals (not the same subject before/after, not paired sites of the same subject, not time series of the same individuals, not technical replicates).
- The comparison must be a biological contrast between subjects: disease vs control, tumor vs normal tissue from different subjects, disease subtype A vs B, genotype KO vs WT, in-vivo treatment/diet/infection vs control, severity or outcome strata. Batch, site, date, cohort name, sequencing run, RNA quality, ancestry alone, sex alone, or age alone are NOT acceptable comparison variables.
- Single-sex studies (all male or all female) ARE usable designs; do not exclude for that reason. Report sex composition in sex_composition instead.
- Studies that mix usable and unusable samples (in-vitro treated subsets, extra time points, repeated visits, other tissues) are usable if a subset filter yields a clean between-subject comparison with at least 3 samples per group; describe the filter in subset_filters. If a subject identifier exists, name it in repeated_measures_subject_key so one sample per subject can be kept.
- When several biological contrasts are possible, list up to 3 in contrasts ordered by priority: prefer the study's primary disease/condition contrast (as stated in title/abstract) over secondary strata, treatment arms, or genotype modifiers; prefer healthy/wild-type/untreated as control.
- Copy comparison keys and level strings EXACTLY as they appear in sample_attributes (lowercase).
- Set scope_flags for: tumor tissue comparisons, cell-free/extracellular RNA, FFPE or degraded material, targeted (non-whole-transcriptome) assays suspected from the title, pooled samples, mixed organisms in one project.
- Be conservative: when the design cannot be determined from the provided metadata, set usable_design = false with reason "insufficient_metadata".
"""

class SubsetFilter(BaseModel):
    action: Literal["keep", "exclude"]
    key: str
    value: str

class Contrast(BaseModel):
    comparison_key: str = Field(..., description="sample attribute key, exactly as given")
    case_level: str
    control_level: str
    note: Optional[str] = None

class DesignVerdict(BaseModel):
    usable_design: bool
    experimental_setting: Literal["in_vivo_tissue", "ex_vivo_primary_cells", "in_vitro_culture_or_stimulation", "cell_line", "organoid_or_stem_cell", "unclear"]
    design_type: Literal["between_subject_case_control", "in_vivo_treatment_arms", "genotype_comparison", "repeated_measures_same_subjects", "time_course", "dose_response", "cell_type_or_region_comparison", "technical_or_batch_comparison", "single_group_no_contrast", "other"]
    exclusion_reason: Optional[Literal["single_cell", "in_vitro_or_cultured", "cell_line", "repeated_measures", "no_biological_contrast", "comparison_is_batch_or_site", "groups_too_small", "insufficient_metadata", "other"]] = None
    contrasts: List[Contrast] = Field(default_factory=list, description="up to 3 candidate contrasts, best first; empty if unusable")
    subset_filters: List[SubsetFilter] = Field(default_factory=list)
    repeated_measures_subject_key: Optional[str] = None
    sex_composition: Literal["both_sexes", "male_only", "female_only", "unknown"]
    scope_flags: List[Literal["tumor_tissue", "cell_free_rna", "ffpe_or_degraded", "targeted_assay_suspected", "pooled_samples", "mixed_organisms"]] = Field(default_factory=list)
    confidence: Literal["high", "medium", "low"]
    rationale: str = Field(..., description="one or two sentences")

def strict_schema(model) -> dict:
    """JSON schema for structured output: additionalProperties=false set on every object (API requirement)."""
    sc = model.model_json_schema()
    def walk(n):
        if isinstance(n, dict):
            if n.get("type") == "object": n.setdefault("additionalProperties", False)
            for v in n.values(): walk(v)
        elif isinstance(n, list):
            for v in n: walk(v)
    walk(sc); return sc

def build_user(p: dict) -> str:
    return "Dataset metadata (JSON):\n" + json.dumps(p, ensure_ascii=False, indent=1) + "\n\nReturn the verdict."

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packets", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="claude-opus-5"); ap.add_argument("--effort", default="medium", choices=["low","medium","high","xhigh","max"])
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--limit", type=int, default=0); ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    packets = [json.loads(l) for l in open(a.packets)]
    if a.limit: packets = packets[:a.limit]
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            try: done.add(json.loads(l)["study"])
            except Exception: pass
    todo = [p for p in packets if p["study"] not in done]
    est_in = sum(len(json.dumps(p)) for p in todo) / 3.5 + len(SYSTEM) / 3.5
    print(f"packets={len(packets)} done={len(done)} todo={len(todo)}  est input tokens≈{est_in:,.0f} (+ ~300 output/each)", file=sys.stderr)
    if a.dry_run:
        price = {"claude-opus-5": (5, 25), "claude-sonnet-5": (2, 10), "claude-haiku-4-5": (1, 5)}.get(a.model, (5, 25))
        print(f"est cost @{a.model}: ${est_in/1e6*price[0] + len(todo)*300/1e6*price[1]:.2f} (system prompt cached after first call → lower)", file=sys.stderr); return
    import anthropic
    client = anthropic.Anthropic(max_retries=4)
    lock = threading.Lock(); os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True); fh = open(a.out, "a")
    def one(p):
        t0 = time.time()
        r = client.messages.parse(model=a.model, max_tokens=4000,
            system=[{"type": "text", "text": SYSTEM, "cache_control": {"type": "ephemeral"}}],
            thinking={"type": "adaptive"}, output_config={"effort": a.effort},
            messages=[{"role": "user", "content": build_user(p)}], output_format=DesignVerdict)
        served_by = a.model
        if r.stop_reason == "refusal":
            # safety-classifier false positive (e.g. infectious-disease metadata): retry the same request with the server-side fallback model
            cat = getattr(r.stop_details, "category", None)
            rb = client.beta.messages.create(model=a.model, max_tokens=4000, betas=["server-side-fallback-2026-07-01"], fallbacks="default",
                system=[{"type": "text", "text": SYSTEM, "cache_control": {"type": "ephemeral"}}],
                thinking={"type": "adaptive"}, output_config={"effort": a.effort, "format": {"type": "json_schema", "schema": strict_schema(DesignVerdict)}},
                messages=[{"role": "user", "content": build_user(p)}])
            if rb.stop_reason == "refusal": raise RuntimeError(f"refusal (fallback too): {cat}")
            text = next(b.text for b in rb.content if b.type == "text"); v = DesignVerdict.model_validate_json(text); r = rb; served_by = f"{rb.model} (fallback after refusal:{cat})"
        else:
            v = r.parsed_output
        d = v.model_dump(); first = d["contrasts"][0] if d["contrasts"] else {}
        rec = {"study": p["study"], "organism": p.get("organism"), **d, "comparison_key": first.get("comparison_key"), "case_level": first.get("case_level"), "control_level": first.get("control_level"), "model": served_by, "effort": a.effort, "prompt_version": "v2",
               "usage": {"in": r.usage.input_tokens, "out": r.usage.output_tokens, "cache_read": r.usage.cache_read_input_tokens, "cache_write": r.usage.cache_creation_input_tokens}, "sec": round(time.time() - t0, 1)}
        with lock: fh.write(json.dumps(rec, ensure_ascii=False) + "\n"); fh.flush()
        return rec
    n_ok = n_err = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one, p): p for p in todo}
        for f in as_completed(futs):
            p = futs[f]
            try: rec = f.result(); n_ok += 1; print(f"[{n_ok+n_err}/{len(todo)}] {rec['study']} usable={rec['usable_design']} {rec['design_type']} key={rec.get('comparison_key')} ({rec['sec']}s)", file=sys.stderr)
            except anthropic.RateLimitError as e: n_err += 1; print(f"rate limit on {p['study']}: {e}", file=sys.stderr); time.sleep(20)
            except Exception as e: n_err += 1; print(f"ERROR {p['study']}: {e}", file=sys.stderr)
    print(f"done ok={n_ok} err={n_err}", file=sys.stderr)

if __name__ == "__main__": main()
