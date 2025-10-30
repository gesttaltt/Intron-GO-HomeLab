#!/usr/bin/env python
import argparse, sys, requests, polars as pl, time
from tqdm import tqdm

# QuickGO annotations endpoint
API = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"

def fetch_go(aspect: str, batch=1000):
    params = {
        "taxonId": "9606",
        "aspect": aspect,           # BP/MF/CC
        "limit": str(batch),
        "includeFields": "goId,goName,geneProductId,symbol,evidenceCode"
    }
    url = API
    page=1
    rows=[]
    with tqdm(desc=f"QuickGO {aspect}", unit="rows") as pbar:
        while True:
            params["page"] = str(page)
            r = requests.get(url, params=params, headers={"Accept":"application/json"})
            if r.status_code != 200:
                print(f"Error {r.status_code}: {r.text}")
            r.raise_for_status()
            j = r.json()
            results = j.get("results", [])
            if not results:
                break
            for it in results:
                rows.append({
                    "go_id": it.get("goId"),
                    "go_name": it.get("goName"),
                    "symbol": it.get("symbol"),
                    "gene_product": it.get("geneProductId"),
                    "evidence": it.get("evidenceCode"),
                })
            pbar.update(len(results))
            page += 1
            # Quick throttle
            time.sleep(0.05)
    return pl.DataFrame(rows).unique()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--aspect", default="BP", choices=["BP","MF","CC"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    df = fetch_go(args.aspect)
    df.write_parquet(args.out)
    print(df.shape)
