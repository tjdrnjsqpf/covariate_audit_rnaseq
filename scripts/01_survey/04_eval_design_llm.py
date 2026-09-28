#!/usr/bin/env python3
"""Evaluate LLM verdicts against manual labels. Usage: python3 04_eval_design_llm.py --pred results/design_llm/labeled_v1.pred.jsonl --labels config/design_labels_v1.tsv"""
import argparse, csv, json, re
from collections import Counter, defaultdict
def norm(s): return re.sub(r"\s+", " ", (s or "").strip().lower())
def level_match(pred, truth):
    p, t = norm(pred), norm(truth)
    if not p or not t: return False
    if t.startswith("^"):   # ground truth written as a regex
        try: return re.search(t, p) is not None
        except re.error: return False
    return p == t or p.startswith(t) or t.startswith(p)
ap = argparse.ArgumentParser(); ap.add_argument("--pred", required=True); ap.add_argument("--labels", required=True); a = ap.parse_args()
lab = {r["study"]: r for r in csv.DictReader(open(a.labels), delimiter="\t")}
# scope exclusions: design is valid but outside the pilot scope -> counted as include in the design-only evaluation
# Match keys tested against the free-text 'reason' column of the original (Korean) hand-label sheet; the Korean literals are kept so that the validation reproduces exactly.
SCOPE = ("all male", "all female", "전원 남성", "전원 여성", "여성 2%", "tumor", "종양", "exrna", "cfrna", "serum", "ffpe", "rna quality", "extracellular", "targeted assay", "sarcoma", "glioma", "t-all", "leukemia", "neoplasm", "melanoma", "grade", "subtype")
for r in lab.values():
    r["label_design"] = "include" if (r["label"] == "include" or any(k in r["reason"].lower() for k in SCOPE)) else "exclude"
pred = {}
for l in open(a.pred):
    r = json.loads(l); pred[r["study"]] = r
rows = [(lab[s], pred[s]) for s in lab if s in pred]
print(f"labeled={len(lab)} predicted={len(pred)} evaluated={len(rows)}")
def contrasts_of(p):
    cs = p.get("contrasts") or []
    if not cs and p.get("comparison_key"): cs = [{"comparison_key": p.get("comparison_key"), "case_level": p.get("case_level"), "control_level": p.get("control_level")}]
    return cs
def any_match(p, l):
    return any(norm(c.get("comparison_key")) == norm(l["group_key"]) and level_match(c.get("case_level"), l["case_level"]) and level_match(c.get("control_level"), l["control_level"]) for c in contrasts_of(p))
for LAB in ["label", "label_design"]:
  print(f"\n######## evaluation against: {LAB} ({'original pilot labels' if LAB == 'label' else 'design-only: scope exclusions counted as include'})")
  for org in ["human", "mouse", "all"]:
      sub = [(l, p) for l, p in rows if org == "all" or l["organism"] == org]
      if not sub: continue
      tp = sum(1 for l, p in sub if l[LAB] == "include" and p["usable_design"]); fn = sum(1 for l, p in sub if l[LAB] == "include" and not p["usable_design"])
      fp = sum(1 for l, p in sub if l[LAB] == "exclude" and p["usable_design"]); tn = sum(1 for l, p in sub if l[LAB] == "exclude" and not p["usable_design"])
      acc = (tp + tn) / len(sub); sens = tp / max(1, tp + fn); spec = tn / max(1, tn + fp)
      inc = [(l, p) for l, p in sub if l["label"] == "include" and p["usable_design"]]
      key_ok = sum(1 for l, p in inc if norm(p.get("comparison_key")) == norm(l["group_key"]))
      lev_ok = sum(1 for l, p in inc if norm(p.get("comparison_key")) == norm(l["group_key"]) and level_match(p.get("case_level"), l["case_level"]) and level_match(p.get("control_level"), l["control_level"]))
      lev_swap = sum(1 for l, p in inc if norm(p.get("comparison_key")) == norm(l["group_key"]) and level_match(p.get("case_level"), l["control_level"]) and level_match(p.get("control_level"), l["case_level"]))
      top3 = sum(1 for l, p in inc if any_match(p, l))
      print(f"== {org}: n={len(sub)}  accuracy={acc:.2f}  sensitivity(include)={sens:.2f} ({tp}/{tp+fn})  specificity(exclude)={spec:.2f} ({tn}/{tn+fp})")
      print(f"   among true includes predicted usable ({len(inc)}): key match={key_ok} ({key_ok/max(1,len(inc)):.2f}), key+levels match(top1)={lev_ok} ({lev_ok/max(1,len(inc)):.2f}), match within top-3 contrasts={top3} ({top3/max(1,len(inc)):.2f}), levels swapped={lev_swap}")
print("\n== disagreements ==")
for l, p in rows:
    if (l["label_design"] == "include") != bool(p["usable_design"]):
        print(f"  {l['study']} [{l['organism']}] truth={l['label_design']} ({l['reason'][:50]}) | pred usable={p['usable_design']} {p['design_type']} {p.get('exclusion_reason')} conf={p['confidence']} :: {p['rationale'][:120]}")
print("\n== key/level mismatches among agreed includes ==")
for l, p in rows:
    if l["label"] == "include" and p["usable_design"] and not any_match(p, l):
        print(f"  {l['study']}: truth [{l['group_key']}] {l['case_level']} vs {l['control_level']} | pred " + " || ".join(f"[{c.get('comparison_key')}] {c.get('case_level')} vs {c.get('control_level')}" for c in contrasts_of(p)))
print("\n== predicted exclusion reasons =="); print(Counter(p.get("exclusion_reason") for _, p in rows if not p["usable_design"]))
