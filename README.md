# GO-Domino-HomeLab (chr21 sprint)

**Goal:** compute a **GO-aware Functional Perturbation Score (FPS)** on 1000G chr21 using GO annotations + exon-overlap SNP density + LOEUF weights. Runs on 16 GB RAM + Ryzen 5 + RTX 3050 (GPU optional).

## ✅ Optimization Complete (2025-10-30)

Successfully implemented **three major optimizations** achieving:
- **120× faster VCF reading** (330 → 40,000 variants/sec)
- **13-18× faster total runtime** (55 min → 3-4 min)
- **135× more GO terms scored** (12 → 1,622)
- **Complete LOEUF constraint weighting** (17,941 genes)

**See:** \`local-reports/optimization-complete.md\` for detailed benchmarks and results.

## Performance Summary

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| VCF Reading | 330 var/sec | 40,000 var/sec | **120× faster** |
| Total Runtime | ~55 minutes | ~3-4 minutes | **13-18× faster** |
| GO Annotations | 14 (sample) | 539,784 (full) | **38,556× more** |
| GO Terms Scored | 12 | 1,622 | **135× more** |
| LOEUF Coverage | 0 genes | 17,941 genes | **New capability** |

## Quickstart (Optimized Pipeline)

\`\`\`bash
# 1) Python venv (conda/mamba optional)
python -m venv venv
source venv/Scripts/activate  # Windows Git Bash
# or: source venv/bin/activate  # Linux/Mac
pip install -r requirements.txt

# 2) Fetch small refs & chr21 VCF (~170MB)
bash scripts/01_get_chr21.sh

# 3) Build exon map from GENCODE GTF
python src/exon_index.py \
  --gtf data/ref/gencode.v49lift37.annotation.gtf.gz \
  --out data/ref/exons_chr21.parquet

# 4) Pull human GO annotations (with GOA fallback)
python src/go_fetch_improved.py \
  --aspect BP \
  --method auto \
  --out data/go/go_bp_full.parquet

# 5) Download gnomAD LOEUF constraint data
# Download from: https://gnomad.broadinstitute.org/downloads
# Place at: data/gnomad/gnomad.v4.1.constraint_metrics.tsv

# 6) Compute GO Functional Perturbation Score (OPTIMIZED)
python src/fps_score_final.py \
  --vcf data/vcf/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  --exons data/ref/exons_chr21.parquet \
  --go data/go/go_bp_full.parquet \
  --loeuf data/gnomad/gnomad.v4.1.constraint_metrics.tsv \
  --chrom 21 \
  --out outputs/fps_chr21_final.parquet

# 7) View results
python -c "import polars as pl; print(pl.read_parquet('outputs/fps_chr21_final.parquet').head(10))"
\`\`\`

**Runtime:** ~3-4 minutes total (includes VCF loading, exon scanning, GO aggregation)

## What it computes

- **Exonic SNP density per gene** (variants per kb), chr21 only (fast).
- Aggregates by GO term using full GO annotations (539,784 annotations).
- **FPS = Σ_gene ( density × w_gene )**, with \`w_gene = 1/LOEUF\` (evolutionary constraint weighting).

**Outputs:** \`outputs/fps_chr21_final.parquet\` with columns:
- \`go_id\`: Gene Ontology term ID
- \`go_name\`: GO term name/description
- \`n_genes\`: Number of genes in this GO term
- \`fps\`: Functional Perturbation Score (higher = more fragile)
- \`mean_density\`: Average variant density across genes
- \`median_loeuf\`: Median LOEUF constraint score

## Results (chr21 Final)

**Data Processed:**
- 1,105,538 variants on chr21
- 41,262 exons from 1,061 genes
- 539,784 GO annotations covering 42,597 genes
- 17,941 genes with LOEUF constraint scores

**Top 10 GO Terms by Functional Fragility Score:** (see \`local-reports/optimization-complete.md\` for full table)

## Technical Optimizations

### 1. Complete GO Annotations
**File:** \`src/go_fetch_improved.py\`
- QuickGO REST API with GOA file fallback
- Successfully parsed 539,784 human GO annotations
- Coverage: 42,597 genes (vs 14 sample annotations)

### 2. gnomAD LOEUF Constraint Weighting
**Data:** \`data/gnomad/gnomad.v4.1.constraint_metrics.tsv\`
- Downloaded gnomAD v4.1 constraint metrics (92 MB)
- Filtered to canonical transcripts (36,315 total)
- Loaded LOEUF scores for 17,941 genes
- Weight formula: \`weight = 1 / LOEUF\`

### 3. Fast Spatial Indexing
**File:** \`src/fps_score_final.py\`
- Binary search on sorted numpy arrays using \`np.searchsorted\`
- O(log N) lookups per exon instead of O(N×M) naive scan
- 120× faster VCF reading (40,000 var/sec)
- Memory efficient: single sorted array reused for all exons

## Available Scripts

**Data Acquisition:**
- \`src/exon_index.py\` - Extract exons from GENCODE GTF
- \`src/go_fetch_improved.py\` - Fetch GO annotations with fallback

**FPS Computation (choose one):**
- \`src/fps_score_simple.py\` - Naive implementation (baseline, ~55 min)
- \`src/fps_score_optimized.py\` - PyRanges approach (partial, compatibility issues)
- \`src/fps_score_final.py\` - **RECOMMENDED: Complete optimized solution (~3-4 min)**

## Sources

- **QuickGO REST API / GOA** for GO annotations
- **GENCODE v49 lift37** GTF (mapped to GRCh37/hg19) for exons
- **1000 Genomes Phase 3** chr21 VCF (+tbi)
- **gnomAD v4.1** constraint (LOEUF) table for evolutionary constraint

## Scale-up

- Change \`01_get_chr21.sh\` to other chromosomes or the full genome
- Switch \`--aspect\` to MF/CC for different GO aspects
- Later: swap SNP density → VEP/SnpEff consequences; add splice-site windows; ClinVar overlay

## Requirements

See \`requirements.txt\` for dependencies. Core: Python 3.9+, polars, numpy, tqdm.

## Cross-Platform Compatibility

- **Windows-compatible:** Pure Python VCF parser, no C compilation required
- **No build tools needed:** Bypassed cyvcf2/htslib dependencies
- **Robust data fetching:** Multiple fallback strategies for GO annotations

---

**Last updated:** 2025-10-30
**Status:** Production-ready for chr21, ready for whole-genome extension
