#!/usr/bin/env python3
"""Reduce the runinfo CSV (which contains malformed lines) to the needed columns as a slim TSV. Usage: python3 57_slim_runinfo.py human|mouse [suffix]"""
import csv, sys
csv.field_size_limit(1<<30)
sp = sys.argv[1]; suf = sys.argv[2] if len(sys.argv) > 2 else "2020plus"
src = f"results/response/sra_runinfo_{sp}_{suf}.csv"; dst = f"results/response/sra_runinfo_{sp}_{suf}.slim.tsv"
keep = ["Run","ReleaseDate","LibraryStrategy","LibrarySource","LibraryLayout","Platform","SRAStudy","BioProject","Sample","BioSample","TaxID","ScientificName","Sex","CenterName"]
seen=set(); n=0; bad=0; hdr=None
with open(src, newline="", errors="replace") as f, open(dst, "w") as o:
    o.write("\t".join(keep)+"\n")
    for row in csv.reader(f):
        if not row: continue
        if row[0]=="Run": hdr={k:i for i,k in enumerate(row)}; continue
        if hdr is None or len(row)!=len(hdr) or not row[0][:3] in ("SRR","ERR","DRR"): bad+=1; continue
        if row[0] in seen: continue
        seen.add(row[0]); n+=1
        o.write("\t".join(row[hdr[k]].replace("\t"," ").replace("\n"," ") for k in keep)+"\n")
print(sp, "kept", n, "bad/dup", bad)
