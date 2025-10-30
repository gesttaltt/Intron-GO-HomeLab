#!/usr/bin/env bash
set -euo pipefail
mkdir -p data/ref data/vcf data/go data/gnomad

echo "[*] Downloading GENCODE v49 GTF (GRCh38)"
curl -L -o data/ref/gencode.v49lift37.annotation.gtf.gz   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh37_mapping/gencode.v49lift37.annotation.gtf.gz

echo "[*] Downloading 1000G Phase3 chr21 VCF (GRCh37/hg19 coords)"
curl -L -o data/vcf/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz   https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
curl -L -o data/vcf/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi   https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

echo "[*] (Optional) Place LOEUF v4.1 TSV at data/gnomad/loeuf_v4_1.tsv.gz (see README)"
echo "[*] Done."
