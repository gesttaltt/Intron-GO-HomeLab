#!/usr/bin/env python
"""Create sample GO annotations for testing when QuickGO is unavailable"""
import polars as pl

# Sample GO annotations for some chr21 genes (real annotations would come from QuickGO)
# These are illustrative examples
sample_data = [
    # Gene APP (Amyloid precursor protein) - chr21 gene
    {"symbol": "APP", "go_id": "GO:0006915", "go_name": "apoptotic process",
     "gene_product": "UniProtKB:P05067", "evidence": "IDA"},
    {"symbol": "APP", "go_id": "GO:0007155", "go_name": "cell adhesion",
     "gene_product": "UniProtKB:P05067", "evidence": "IEA"},
    {"symbol": "APP", "go_id": "GO:0016485", "go_name": "protein processing",
     "gene_product": "UniProtKB:P05067", "evidence": "TAS"},

    # Gene SOD1 (Superoxide dismutase) - chr21 gene
    {"symbol": "SOD1", "go_id": "GO:0006979", "go_name": "response to oxidative stress",
     "gene_product": "UniProtKB:P00441", "evidence": "IDA"},
    {"symbol": "SOD1", "go_id": "GO:0055114", "go_name": "oxidation-reduction process",
     "gene_product": "UniProtKB:P00441", "evidence": "IDA"},

    # Gene RUNX1 (Runt-related transcription factor 1) - chr21 gene
    {"symbol": "RUNX1", "go_id": "GO:0006355", "go_name": "regulation of transcription",
     "gene_product": "UniProtKB:Q01196", "evidence": "IDA"},
    {"symbol": "RUNX1", "go_id": "GO:0030099", "go_name": "myeloid cell differentiation",
     "gene_product": "UniProtKB:Q01196", "evidence": "IDA"},

    # Gene DYRK1A - chr21 gene
    {"symbol": "DYRK1A", "go_id": "GO:0006468", "go_name": "protein phosphorylation",
     "gene_product": "UniProtKB:Q13627", "evidence": "IDA"},
    {"symbol": "DYRK1A", "go_id": "GO:0007399", "go_name": "nervous system development",
     "gene_product": "UniProtKB:Q13627", "evidence": "TAS"},

    # Gene HLCS - chr21 gene
    {"symbol": "HLCS", "go_id": "GO:0009102", "go_name": "biotin metabolic process",
     "gene_product": "UniProtKB:P50747", "evidence": "IDA"},
    {"symbol": "HLCS", "go_id": "GO:0018271", "go_name": "protein-biotin ligase activity",
     "gene_product": "UniProtKB:P50747", "evidence": "IDA"},

    # Add more generic terms that might apply to multiple genes
    {"symbol": "APP", "go_id": "GO:0008150", "go_name": "biological process",
     "gene_product": "UniProtKB:P05067", "evidence": "IEA"},
    {"symbol": "SOD1", "go_id": "GO:0008150", "go_name": "biological process",
     "gene_product": "UniProtKB:P00441", "evidence": "IEA"},
    {"symbol": "RUNX1", "go_id": "GO:0008150", "go_name": "biological process",
     "gene_product": "UniProtKB:Q01196", "evidence": "IEA"},
]

df = pl.DataFrame(sample_data)
df.write_parquet("data/go/go_bp_human.parquet")
print(f"Created sample GO data: {df.shape}")
print(df)
