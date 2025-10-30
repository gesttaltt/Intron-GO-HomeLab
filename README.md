# GO‑Domino‑HomeLab (chr21 sprint)

**Goal (today):** compute a **GO‑aware Functional Perturbation Score (FPS)** on 1000G chr21 using GO annotations + exon‑overlap SNP density + LOEUF weights. Runs on 16 GB RAM + Ryzen 5 + RTX 3050 (GPU optional).

## Quickstart
```bash
# 1) Conda env
mamba env create -f environment.yml && conda activate go-domino

# 2) Fetch small refs & chr21 VCF (~170MB)
bash scripts/01_get_chr21.sh

# 3) Build exon map from GENCODE GTF
python src/exon_index.py --gtf data/ref/gencode.v49lift37.annotation.gtf.gz --out data/ref/exons_chr21.parquet

# 4) Pull human GO annotations (QuickGO) & build gene↔GO map
python src/go_fetch.py --aspect BP --out data/go/go_bp_human.parquet

# 5) Compute GO Functional Perturbation Score (FPS)
python src/fps_score.py   --vcf data/vcf/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz   --exons data/ref/exons_chr21.parquet   --go data/go/go_bp_human.parquet   --loeuf data/gnomad/loeuf_v4_1.tsv.gz   --out outputs/fps_chr21.parquet

# 6) View top‑fragile GO terms
python src/report_top.py --fps outputs/fps_chr21.parquet --top 25
```
**Note:** Place a LOEUF/constraint TSV in `data/gnomad/loeuf_v4_1.tsv.gz` (see README Sources). You can run without LOEUF using `--no-constraint`.

## What it computes
- Exonic SNP density per gene (variants per kb), chr21 only (fast).
- Aggregates by GO term using the GO DAG (descendants included, dedup by gene).
- **FPS = Σ_gene ( density * w_gene )**, with `w_gene = 1/LOEUF` (or 1 if not provided).

Outputs: `outputs/fps_chr21.parquet` with columns: `go_id,go_name,n_genes,fps,mean_density,median_loeuf`.

## Sources (fetch via scripts/ docs)
- QuickGO REST API for GO annotations.
- GENCODE v49 **lift37** GTF (mapped to GRCh37/hg19) for exons.
- 1000 Genomes Phase 3 chr21 VCF (+tbi).
- gnomAD constraint (LOEUF) table v4.1 (optional).

## Scale‑up
- Change `01_get_chr21.sh` to other chromosomes or the full set.
- Switch `--aspect` to MF/CC, or filter GO terms by depth/info‑content.
- Later: swap SNP density → VEP/SnpEff consequences; add splice‑site windows; ClinVar overlay.
