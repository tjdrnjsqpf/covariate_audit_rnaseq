#!/usr/bin/env python3
"""Fetch SRA runinfo for human and mouse bulk RNA-seq runs released 2020-2026 (NCBI E-utilities; no API key needed, 3 req/s).
Output: results/response/sra_runinfo_{human,mouse}_<suffix>.csv. Usage: python3 55_sra_post2020_runinfo.py <human|mouse> [from_date to_date suffix] (default 2020/01/01 2026/12/31 2020plus)"""
import sys, time, urllib.request, urllib.parse, os, re
ORG = sys.argv[1]; name = {"human": "Homo sapiens", "mouse": "Mus musculus"}[ORG]
FROM, TO, SUF = (sys.argv[2], sys.argv[3], sys.argv[4]) if len(sys.argv) > 4 else ("2020/01/01", "2026/12/31", "2020plus")
OUT = f"results/response/sra_runinfo_{ORG}_{SUF}.csv"; os.makedirs("results/response", exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
term = f'"{name}"[Organism] AND "rna seq"[Strategy] AND "transcriptomic"[Source] AND ("{FROM}"[PDAT] : "{TO}"[PDAT])'
def get(url, tries=10):
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r: return r.read()
        except Exception as e:
            time.sleep(5 * (t + 1)); err = e
    raise err
q = get(E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "sra", "term": term, "usehistory": "y", "retmax": 0})).decode()
count = int(re.search(r"<Count>(\d+)</Count>", q).group(1)); web = re.search(r"<WebEnv>(\S+)</WebEnv>", q).group(1); key = re.search(r"<QueryKey>(\d+)</QueryKey>", q).group(1)
print(ORG, "runs:", count, flush=True)
step = 10000; first = True; start0 = 0; mode = "wb"
if os.path.exists(OUT) and os.path.getsize(OUT) > 0:   # resume: keep whole 10k chunks already written (runs are counted as lines starting with an accession)
    with open(OUT, "rb") as f: done = sum(1 for l in f if l[:3] in (b"SRR", b"ERR", b"DRR"))
    keep = (done // step) * step
    if keep > 0:
        with open(OUT, "rb") as f: lines = f.readlines()
        hdr = [l for l in lines if l.startswith(b"Run,")][:1]; body = [l for l in lines if l[:3] in (b"SRR", b"ERR", b"DRR")][:keep]
        with open(OUT, "wb") as f: f.write(b"".join(hdr + body))
        start0 = keep; first = False; mode = "ab"; print(f"resuming at {start0} (kept {keep} runs)", flush=True)
with open(OUT, mode) as f:
    for start in range(start0, count, step):
        b = get(E + "efetch.fcgi?" + urllib.parse.urlencode({"db": "sra", "query_key": key, "WebEnv": web, "retstart": start, "retmax": step, "rettype": "runinfo", "retmode": "text"}))
        lines = b.split(b"\n")
        if not first: lines = [l for l in lines if l and not l.startswith(b"Run,")]
        f.write(b"\n".join(l for l in lines if l) + b"\n"); first = False
        print(f"  {start + step}/{count}", flush=True); time.sleep(0.4)
print("done", OUT, flush=True)
