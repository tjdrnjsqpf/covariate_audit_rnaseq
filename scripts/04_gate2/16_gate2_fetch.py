#!/usr/bin/env python3
"""Gate 2-a: locate each study's original paper (GEO→PubMed, Europe PMC accession search) and extract the Methods section from the OA full text.
Usage: python 16_gate2_fetch.py --org human --out results/gate2/fulltext_human.jsonl [--limit N]
Studies already processed are skipped."""
import argparse, json, os, re, sys, time, csv
import requests
import xml.etree.ElementTree as ET
PROJ = "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
S = requests.Session(); S.headers["User-Agent"] = "covariate-audit/1.0 (research; contact via GitHub)"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"; EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
def get(url, params=None, tries=4, sleep=0.4):
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=60); time.sleep(sleep)
            if r.status_code == 200: return r
            if r.status_code == 429: time.sleep(5 * (i + 1)); continue
            if r.status_code == 404: return None
        except requests.RequestException: time.sleep(2 * (i + 1))
    return None
def geo_pmids(gsm):
    r = get(f"{EUTILS}/esearch.fcgi", {"db": "gds", "term": f"{gsm}[ACCN]", "retmode": "json"})
    ids = (r.json().get("esearchresult", {}).get("idlist", []) if r else [])
    if not ids: return [], None
    r = get(f"{EUTILS}/esummary.fcgi", {"db": "gds", "id": ",".join(ids[:10]), "retmode": "json"})
    if not r: return [], None
    pm, gse = [], None
    for k, v in r.json().get("result", {}).items():
        if k == "uids": continue
        pm += [str(p) for p in v.get("pubmedids", [])]
        if str(v.get("accession", "")).startswith("GSE"): gse = v["accession"]
    return sorted(set(pm)), gse
def epmc_search(q):
    r = get(f"{EPMC}/search", {"query": q, "format": "json", "pageSize": 5, "resultType": "lite"})
    return (r.json().get("resultList", {}).get("result", []) if r else [])
def pmid_to_pmcid(pmid):
    res = epmc_search(f"EXT_ID:{pmid} AND SRC:MED")
    for x in res:
        if x.get("pmcid"): return x["pmcid"], x.get("title"), x.get("isOpenAccess") == "Y"
    return None, (res[0].get("title") if res else None), False
def fulltext_sections(pmcid):
    r = get(f"{EPMC}/{pmcid}/fullTextXML")
    if not r or not r.text.strip().startswith("<"): return None, None
    try: root = ET.fromstring(r.content)
    except ET.ParseError: return None, None
    def txt(el): return re.sub(r"\s+", " ", " ".join(el.itertext())).strip()
    abstract = " ".join(txt(a) for a in root.iter("abstract"))
    keep, other = [], []
    for sec in root.iter("sec"):
        t = sec.find("title"); title = txt(t) if t is not None else ""
        body = txt(sec)
        if re.search(r"method|material|statistic|analys|differential|expression|bioinformatic|rna.?seq|sequencing|processing", title, re.I): keep.append(f"[{title}] {body}")
        else: other.append(f"[{title}] {body}")
    methods = "\n\n".join(keep) if keep else "\n\n".join(other)
    return abstract[:3000], methods[:30000]
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--org", default="human"); ap.add_argument("--out", required=True); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); os.chdir(PROJ); os.makedirs(os.path.dirname(a.out), exist_ok=True)
    summ = {r["study"]: r for r in csv.DictReader(open(f"results/main_{a.org}_summary_v2.tsv"), delimiter="\t")}
    ps_file = "results/project_summary_all.tsv" if a.org == "human" else "results/mouse_project_summary_sexinferred.tsv"
    ps = {r["study"]: r for r in csv.DictReader(open(ps_file), delimiter="\t")}
    studies = [s for s, r in summ.items() if r.get("valid") == "TRUE"]
    if a.limit: studies = studies[:a.limit]
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            try: done.add(json.loads(l)["study"])
            except Exception: pass
    fh = open(a.out, "a"); n = 0
    for s in studies:
        if s in done: continue
        n += 1; p = ps.get(s, {}); gsm = p.get("gsm_first", ""); title = p.get("study_title", "")
        rec = {"study": s, "gsm": gsm, "study_title": title, "route": None, "pmids": [], "pmcid": None, "paper_title": None, "open_access": False, "abstract": None, "methods": None}
        try:
            pm, gse = ([], None)
            if gsm and gsm.startswith("GSM"): pm, gse = geo_pmids(gsm); rec["gse"] = gse
            if pm: rec["route"] = "geo_pubmed"
            else:
                res = epmc_search(f'"{s}"' + (f' OR "{gse}"' if gse else ""))
                pm = [x["pmid"] for x in res if x.get("pmid")]
                if pm: rec["route"] = "epmc_accession"
                elif title and len(title) > 25:
                    res = epmc_search(f'TITLE:"{title[:200]}"'); pm = [x["pmid"] for x in res[:1] if x.get("pmid")]
                    if pm: rec["route"] = "epmc_title"
            rec["pmids"] = pm[:5]
            for pmid in rec["pmids"]:
                pmcid, ptitle, oa = pmid_to_pmcid(pmid); rec["paper_title"] = rec["paper_title"] or ptitle
                if pmcid:
                    ab, me = fulltext_sections(pmcid)
                    if me: rec.update(pmcid=pmcid, open_access=True, abstract=ab, methods=me, pmid_used=pmid); break
        except Exception as e: rec["error"] = str(e)[:200]
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n"); fh.flush()
        print(f"[{n}] {s} route={rec['route']} pmids={len(rec['pmids'])} fulltext={'Y' if rec['methods'] else 'N'}", file=sys.stderr)
    print("done", file=sys.stderr)
if __name__ == "__main__": main()
