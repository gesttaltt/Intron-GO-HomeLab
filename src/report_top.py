#!/usr/bin/env python
import argparse, polars as pl
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--fps", required=True)
    ap.add_argument("--top", type=int, default=25)
    args=ap.parse_args()
    df=pl.read_parquet(args.fps).sort("fps", descending=True).head(args.top)
    print(df)
