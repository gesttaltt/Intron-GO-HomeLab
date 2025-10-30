#!/usr/bin/env python
"""
Improved GO annotation fetcher with multiple fallback strategies:
1. Try QuickGO REST API
2. Fall back to GOA file download
3. Parse GAF format
"""
import argparse, sys, requests, polars as pl, time, gzip, io
from tqdm import tqdm
from pathlib import Path

# QuickGO annotations endpoint
API = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
GOA_URL = "http://geneontology.org/gene-associations/goa_human.gaf.gz"

def fetch_go_api(aspect: str, batch=100, max_results=50000):
    """Try to fetch from QuickGO API with improved parameters"""
    params = {
        "taxonId": "9606",
        "aspect": aspect,
        "limit": batch,
    }
    url = API
    page = 1
    rows = []

    print(f"Attempting QuickGO API for aspect={aspect}...")

    with tqdm(desc=f"QuickGO {aspect}", unit="rows", total=max_results) as pbar:
        while len(rows) < max_results:
            params["page"] = page
            try:
                r = requests.get(url, params=params, headers={"Accept": "application/json"}, timeout=10)

                if r.status_code != 200:
                    print(f"\nAPI Error {r.status_code}: {r.text[:200]}")
                    return None

                j = r.json()
                results = j.get("results", [])

                if not results:
                    print(f"\nNo more results at page {page}")
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
                time.sleep(0.1)  # Be nice to the API

            except requests.exceptions.RequestException as e:
                print(f"\nNetwork error: {e}")
                return None
            except Exception as e:
                print(f"\nUnexpected error: {e}")
                return None

    if rows:
        print(f"\nSuccessfully fetched {len(rows)} annotations from API")
        return pl.DataFrame(rows).unique()
    return None

def parse_gaf(gaf_path):
    """Parse GAF (Gene Association File) format"""
    rows = []

    print(f"Parsing GAF file: {gaf_path}")

    with gzip.open(gaf_path, 'rt') as f:
        for line in tqdm(f, desc="Parsing GAF"):
            if line.startswith('!'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 15:
                continue

            # GAF 2.1 format:
            # 0: DB, 1: DB_Object_ID, 2: DB_Object_Symbol, 3: Qualifier
            # 4: GO_ID, 5: DB_Reference, 6: Evidence_Code, 7: With/From
            # 8: Aspect, 9: DB_Object_Name, 10: DB_Object_Synonym
            # 11: DB_Object_Type, 12: Taxon, 13: Date, 14: Assigned_By

            taxon = parts[12]
            if 'taxon:9606' not in taxon.lower():
                continue  # Only human

            rows.append({
                "go_id": parts[4],
                "go_name": parts[9] if len(parts[9]) > 0 else "?",
                "symbol": parts[2],
                "gene_product": f"{parts[0]}:{parts[1]}",
                "evidence": parts[6],
            })

    print(f"Parsed {len(rows)} human annotations from GAF")
    return pl.DataFrame(rows).unique()

def download_goa(output_dir="data/go"):
    """Download GOA file"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gaf_path = output_dir / "goa_human.gaf.gz"

    if gaf_path.exists():
        print(f"GOA file already exists: {gaf_path}")
        return gaf_path

    print(f"Downloading GOA file from {GOA_URL}...")

    try:
        r = requests.get(GOA_URL, stream=True, timeout=30)
        r.raise_for_status()

        total_size = int(r.headers.get('content-length', 0))

        with open(gaf_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                    pbar.update(len(chunk))

        print(f"Downloaded to {gaf_path}")
        return gaf_path

    except Exception as e:
        print(f"Download failed: {e}")
        return None

def fetch_go(aspect: str, method="auto"):
    """
    Fetch GO annotations with fallback strategy

    method: 'api', 'goa', or 'auto' (try API then GOA)
    """

    if method in ["auto", "api"]:
        df = fetch_go_api(aspect)
        if df is not None:
            return df

        if method == "api":
            raise SystemExit("API method failed and no fallback allowed")

    # Fallback to GOA file
    print("Falling back to GOA file download...")
    gaf_path = download_goa()

    if gaf_path is None:
        raise SystemExit("Both API and GOA download failed")

    df = parse_gaf(gaf_path)

    # Filter by aspect if specified
    if aspect != "ALL":
        # Note: GAF aspect codes are F (MF), P (BP), C (CC)
        aspect_map = {"BP": "P", "MF": "F", "CC": "C"}
        # This is approximate - GAF column 8 has aspect, but we don't parse it above
        # For now, return all and let user filter
        print(f"Warning: GAF parsing returns all aspects. Filtering by aspect not implemented yet.")

    return df

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--aspect", default="BP", choices=["BP", "MF", "CC", "ALL"])
    ap.add_argument("--method", default="auto", choices=["auto", "api", "goa"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = fetch_go(args.aspect, args.method)
    df.write_parquet(args.out)
    print(f"Final output: {df.shape}")
    print(df.head(10))
