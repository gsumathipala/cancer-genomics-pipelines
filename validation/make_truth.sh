#!/usr/bin/env bash
# Created by Brainstorm, 2026.
#
# Build every known-answer dataset used to verify the pipeline end to end.
#
#   validation/make_truth.sh OUT_DIR
#
# Run inside the cancer_pipeline environment (it needs samtools and wgsim):
#   conda activate cancer_pipeline && validation/make_truth.sh ~/truth
#
# Produces, under OUT_DIR:
#   dna/fastq/        TUMOUR_S1 and NORMAL_S2 -- a matched pair
#   dna/lanesplit/    LANE_S3: the same tumour reads split into four lanes
#   dna/panel.bed     the five target regions
#   dna/truth.tsv     the calls a correct run must make
#   rna/fastq/        RNA_T1_S1 carrying EML4::ALK and BCR::ABL1
#   rna/truth.tsv     the fusions and their exact breakpoints
#
# Every allele fraction is exact BY CONSTRUCTION, from how the haplotypes
# are mixed -- which is what lets check_truth.py grade a run rather than
# describe it.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:?usage: make_truth.sh OUT_DIR}"
mkdir -p "$OUT/dna/fastq" "$OUT/dna/lanesplit" "$OUT/rna/fastq"

# Paired-end 150 bp, ~350 bp fragments, a low error rate, and NO random
# mutations of wgsim's own (-r 0 -R 0): the only variants are the planted
# ones.
sim() {  # fasta pairs seed prefix fragment
  wgsim -N "$2" -1 150 -2 150 -d "${5:-350}" -s 35 -e 0.001 -r 0 -R 0 -X 0 \
        -S "$3" "$1" "$4.1.fq" "$4.2.fq" >/dev/null 2>&1
  # Unique names across haplotypes, or pairs from two files collide.
  sed -i "1~4s/^@/@$4_/" "$4.1.fq" "$4.2.fq"
}

echo "== DNA: haplotypes, panel and truth set"
python3 "$HERE/build_dna_truth.py" "$OUT/dna"
cd "$OUT/dna"
# Diploid mixing:
#   normal = 50% reference + 50% germline         -> germline SNP at 50%
#   tumour = 50% reference + 30% tumour + 20% germ -> somatic 30%, germline 50%
sim ref_hap.fa  10000 11 tR; sim tum_hap.fa 6000 12 tT; sim germ_hap.fa 4000 13 tG
sim ref_hap.fa   7000 21 nR; sim germ_hap.fa 7000 22 nG
cat tR.1.fq tT.1.fq tG.1.fq | gzip > fastq/TUMOUR_S1_R1_001.fastq.gz
cat tR.2.fq tT.2.fq tG.2.fq | gzip > fastq/TUMOUR_S1_R2_001.fastq.gz
cat nR.1.fq nG.1.fq | gzip > fastq/NORMAL_S2_R1_001.fastq.gz
cat nR.2.fq nG.2.fq | gzip > fastq/NORMAL_S2_R2_001.fastq.gz
rm -f ./*.fq

echo "== DNA: the same tumour, split into four lanes"
for R in 1 2; do
  zcat "fastq/TUMOUR_S1_R${R}_001.fastq.gz" | split -l 20000 -d -a 1 - "lanesplit/part_R${R}_"
  for i in 0 1 2 3; do
    gzip -c "lanesplit/part_R${R}_$i" > "lanesplit/LANE_S3_L00$((i + 1))_R${R}_001.fastq.gz"
  done
  rm -f lanesplit/part_R${R}_*
done

echo "== RNA: transcripts, fusions and truth set"
python3 "$HERE/build_rna_truth.py" "$OUT/rna"
cd "$OUT/rna"
sim transcripts_bg.fa     600000 31 bg 300
sim transcripts_fusion.fa   4000 32 fu 300
cat bg.1.fq fu.1.fq | gzip > fastq/RNA_T1_S1_R1_001.fastq.gz
cat bg.2.fq fu.2.fq | gzip > fastq/RNA_T1_S1_R2_001.fastq.gz
rm -f ./*.fq

echo "== done: $OUT"
