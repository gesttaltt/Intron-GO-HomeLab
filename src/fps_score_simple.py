#!/usr/bin/env python
import argparse, polars as pl, gzip
from collections import defaultdict
from tqdm import tqdm
import math

def load_exons(path):
    ex=pl.read_parquet(path)
    ex=ex.with_columns((pl.col("end")-pl.col("start")+1).alias("length"))
    return ex

def load_go(path):
    go=pl.read_parquet(path)
    return go

def load_loeuf(path):
    if path is None:
        return {}
    df=None
    if path.endswith(".gz"):
        df=pl.read_csv(path, separator="\t")
    else:
        df=pl.read_csv(path, separator="\t")
    cols=df.columns
    cand=[c for c in cols if "gene" in c.lower() and "symbol" in c.lower()]
    gene_col=cand[0] if cand else "gene"
    loeuf_col=None
    for k in ["loeuf","lof.oe_ci.upper","oe_lof_upper"]:
        if k in cols:
            loeuf_col=k
            break
    if loeuf_col is None:
        raise SystemExit("Could not find LOEUF column in file.")
    m={}
    for row in df.select([gene_col, loeuf_col]).iter_rows():
        g, v = row[0], row[1]
        try:
            m[str(g)] = float(v)
        except:
            continue
    return m

class SimpleVCF:
    """Simple VCF reader using gzip"""
    def __init__(self, path):
        self.path = path

    def __iter__(self):
        with gzip.open(self.path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                yield SimpleVariant(parts)

class SimpleVariant:
    def __init__(self, parts):
        self.CHROM = parts[0]
        self.POS = int(parts[1])
        self.ID = parts[2]
        self.REF = parts[3]
        self.ALT = parts[4]
        self.QUAL = parts[5]
        self.FILTER = parts[6]
        self.INFO = parts[7]

def exonic_density(vcf_path, exons_df):
    # Build in-memory index per chrom for quick scan
    exon_tbl = exons_df.to_dicts()
    total_len_by_gene=defaultdict(int)
    for r in exon_tbl:
        total_len_by_gene[r["gene_name"]] += int(r["length"])
    counts=defaultdict(int)

    # Use simple VCF reader instead of cyvcf2
    vcf=SimpleVCF(vcf_path)

    # Count lines for progress bar
    print("Counting variants...")
    total_vars = 0
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                total_vars += 1

    print(f"Processing {total_vars} variants...")
    for var in tqdm(vcf, total=total_vars, desc="VCF variants"):
        pos=var.POS
        # naive O(N) over exons per line would be slow; prefilter by chrom
        # here we just restrict to chr21 files, so acceptable for sprint
        for r in exon_tbl:
            if pos >= r["start"] and pos <= r["end"]:
                counts[r["gene_name"]] += 1
    density={g: (counts[g]/(L/1000.0)) if L>0 else 0.0 for g,L in total_len_by_gene.items()}
    return density

def main(args):
    ex=load_exons(args.exons)
    go=load_go(args.go)
    lo=load_loeuf(args.loeuf) if not args.no_constraint else {}
    dens=exonic_density(args.vcf, ex)
    # Map gene->GO terms present
    go_map={}
    for row in go.select(["symbol","go_id"]).unique().iter_rows():
        sym, gid = row
        go_map.setdefault(sym, set()).add(gid)
    # Score
    fps=defaultdict(float)
    genes_by_go=defaultdict(set)
    for gene, d in dens.items():
        terms = go_map.get(gene, [])
        w = 1.0
        if gene in lo:
            try:
                w = 1.0/float(lo[gene]) if lo[gene]>0 else 1.0
            except: pass
        for t in terms:
            fps[t] += d * w
            genes_by_go[t].add(gene)
    # Build frame with names
    go_names={}
    if "go_name" in go.columns:
        for r in go.select(["go_id","go_name"]).unique().iter_rows():
            go_names[r[0]] = r[1]
    rows=[]
    for gid,score in fps.items():
        gs = genes_by_go[gid]
        avg = sum(dens.get(g,0.0) for g in gs)/max(1,len(gs))
        med_lo = None
        if lo:
            vals=[lo[g] for g in gs if g in lo]
            med_lo = float(sorted(vals)[len(vals)//2]) if vals else None
        rows.append({"go_id":gid,"go_name":go_names.get(gid,"?"),"n_genes":len(gs),"fps":score,"mean_density":avg,"median_loeuf":med_lo})
    pl.DataFrame(rows).sort("fps", descending=True).write_parquet(args.out)
    print(f"Wrote {args.out}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--go", required=True)
    ap.add_argument("--loeuf", default=None)
    ap.add_argument("--no-constraint", action="store_true")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()
    main(args)
