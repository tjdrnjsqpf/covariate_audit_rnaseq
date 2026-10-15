#!/usr/bin/env python3
"""Advisor comment Q5, step 2: fetch BioSample attributes (E-utilities) for 2020+ studies with sex recorded (>=90% of samples, 8-200 samples) and write a long TSV.
Input: results/response/M_projects_{sp}_2020plus.tsv (58), slim runinfo (57). Output: results/response/biosample_attrs_{sp}_2020plus.tsv.gz
Usage: python3 59_biosample_attrs.py <human|mouse>  (3 requests/s, batches of 200 BioSamples)"""
import sys, time, urllib.request, urllib.parse, os, gzip, csv, xml.etree.ElementTree as ET
sp = sys.argv[1]; E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
elig = {r["study"] for r in csv.DictReader(open(f"results/response/M_projects_{sp}_2020plus.tsv"), delimiter="\t") if r["elig"] == "TRUE"}
bs = {}
for r in csv.DictReader(open(f"results/response/sra_runinfo_{sp}_2020plus.slim.tsv"), delimiter="\t", quoting=csv.QUOTE_NONE):
    if r["SRAStudy"] in elig and r["BioSample"].startswith("SAM"): bs.setdefault(r["BioSample"], r["SRAStudy"])
ids = sorted(bs, key=lambda a: (not a.startswith("SAMN"), a)); OUT = f"results/response/biosample_attrs_{sp}_2020plus.tsv.gz"; TMP = OUT + ".tmp"
print(sp, "studies", len(elig), "biosamples", len(ids), flush=True)
def get(url, data, tries=6):
    for t in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data.encode()), timeout=300) as r: return r.read()
        except Exception as e: time.sleep(5 * (t + 1)); err = e
    raise err
step = 200; nattr = 0
with gzip.open(TMP, "wt") as o:
    o.write("study\tbiosample\tkey\tvalue\n")
    for i in range(0, len(ids), step):
        chunk = ids[i:i+step]
        if chunk[0].startswith("SAMN"):   # NCBI accessions can be passed as id=; SAMD/SAMEA numerics would be misread as SAMN uids, so resolve them through esearch first
            x = get(E + "efetch.fcgi", urllib.parse.urlencode({"db": "biosample", "id": ",".join(chunk), "retmode": "xml"}))
        else:
            q = get(E + "esearch.fcgi", urllib.parse.urlencode({"db": "biosample", "term": " OR ".join(a + "[accn]" for a in chunk), "usehistory": "y", "retmax": 0})).decode()
            import re; w = re.search(r"<WebEnv>(\S+)</WebEnv>", q); k = re.search(r"<QueryKey>(\d+)</QueryKey>", q)
            if not (w and k): print("  esearch failed at", i, flush=True); time.sleep(2); continue
            time.sleep(0.35); x = get(E + "efetch.fcgi", urllib.parse.urlencode({"db": "biosample", "query_key": k.group(1), "WebEnv": w.group(1), "retmax": 1000, "retmode": "xml"}))
        try: root = ET.fromstring(x)
        except ET.ParseError: print("  parse error at", i, flush=True); time.sleep(2); continue
        for b in root.iter("BioSample"):
            acc = b.get("accession") or ""
            if acc not in bs: continue
            for a in b.iter("Attribute"):
                k = (a.get("harmonized_name") or a.get("attribute_name") or "").strip().lower(); v = (a.text or "").strip()
                if k and v: o.write(f"{bs[acc]}\t{acc}\t{k}\t" + v.replace("\t", " ").replace("\n", " ") + "\n"); nattr += 1
        if (i // step) % 20 == 0: print(f"  {i+len(chunk)}/{len(ids)} attrs {nattr}", flush=True)
        time.sleep(0.35)
os.replace(TMP, OUT); print("done", OUT, nattr, flush=True)
