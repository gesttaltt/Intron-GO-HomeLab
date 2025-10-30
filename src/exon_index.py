#!/usr/bin/env python
import argparse, polars as pl, gzip, io

def parse_gtf(gtf_path, chrom="21"):
    cols = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    # Read as csv with tab separator; skip comments
    df = pl.read_csv(
        gtf_path, 
        separator="\t", 
        has_header=False, 
        comment_prefix="#", 
        new_columns=cols, 
        null_values="."
    )
    df = df.filter((pl.col("feature")=="exon") & (pl.col("seqname").str.replace("^chr","").eq(chrom)))
    # Extract gene_name from attribute field
    gene = (pl.col("attribute")
            .str.extract(r'gene_name "([^"]+)"', 1)
            .alias("gene_name"))
    exon_id = (pl.col("attribute")
            .str.extract(r'exon_id "([^"]+)"', 1)
            .alias("exon_id"))
    return df.select([pl.col("seqname").alias("chrom"),
                      pl.col("start").cast(pl.Int64),
                      pl.col("end").cast(pl.Int64),
                      pl.col("strand"),
                      gene,
                      exon_id])

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chrom", default="21")
    args=ap.parse_args()
    df=parse_gtf(args.gtf, args.chrom)
    df.write_parquet(args.out)
    print(df.shape)
