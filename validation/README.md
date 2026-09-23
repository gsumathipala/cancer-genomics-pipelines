<!-- Created by Brainstorm, 2026. -->

# Known-answer validation

**Research use only.** Passing this proves the pipeline finds what was
planted in *simulated* reads. It is not a clinical validation: that needs
real specimens with orthogonally confirmed results. What this does and does
not prove is set out in [KNOWN_LIMITATIONS.md](../KNOWN_LIMITATIONS.md),
L-01 to L-05.

The unit tests (`run_tests.py`) check what this bundle decides. This checks
the whole thing — every tool, every handoff, every report — by running it on
reads whose correct answer is known in advance, and grading what comes out.

```bash
conda activate cancer_pipeline
validation/make_truth.sh ~/truth            # ~1 minute; needs samtools, wgsim
# ... run any pathway on ~/truth/dna or ~/truth/rna, then:
python3 validation/check_truth.py dna ~/truth/dna/truth.tsv  <run>/mutect2/<T>.mutect2.filtered.vcf.gz
python3 validation/check_truth.py dna ~/truth/dna/truth.tsv  <run>/mutect2/<T>.mutect2.filtered.vcf.gz --tumour-only
python3 validation/check_truth.py rna ~/truth/rna/truth.tsv  <run>/fusions/<S>.fusions.tsv
```

`check_truth.py` exits 0 only on a clean pass, and reports every miss and
every extra call. It finds the tumour's column from Mutect2's own
`##tumor_sample` header rather than assuming it — a swapped tumour and
normal is one of the failures it exists to catch.

## What is planted

| DNA — built from real hg38 | Constructed VAF |
|---|---|
| BRAF V600E, KRAS G12D, TP53 R273H, EGFR L858R, PIK3CA H1047R | 0.30 |
| a 3 bp deletion in TP53 (intronic: tests indel calling) | 0.30 |
| a heterozygous germline SNP near KRAS, in tumour **and** normal | 0.50 |

Every coordinate is checked against the reference base before anything is
written. A matched-normal run must **not** call the germline SNP; a
tumour-only run will, at ~0.50, because nothing distinguishes a novel
germline variant from a somatic one without a normal.

`dna/lanesplit/` holds the same tumour reads split into four lanes. A run
on it must reach the same depth as the unsplit reads — a run on one lane
reaches a quarter of it.

| RNA — real GENCODE transcripts on real hg38 | Breakpoints |
|---|---|
| EML4::ALK variant 1 (EML4 ex13 → ALK ex20) | chr2:42295516 → chr2:29223528 |
| BCR::ABL1 e14a2 / p210 (BCR ex14 → ABL1 ex2) | chr22:23290413 → chr9:130854064 |

on a background of ~300 canonical transcripts plus the wild-type partners,
because Arriba refuses a library with no normal reads.

## Results, 23 September 2026

Every pathway, run for real on this machine, graded by `check_truth.py`:

| Pathway | Mode | Result |
|---|---|---|
| Web, DNA | matched normal, auto-discover, tumour named | PASS · PCGR: all five drivers tier 1 |
| Web, DNA | custom panel BED, no vendor profile | PASS · coverage made, TMB over 0.0023 Mb |
| Command line, DNA | tumour-only, every resource, PCGR at step 13 | PASS |
| Command line, DNA | matched normal | PASS |
| Orchestrator, DNA | matched normal, stage 2 → stage 3 | PASS |
| Web, RNA | auto-discover | PASS · both fusions, exact breakpoints |
| Command line, RNA | auto-discover | PASS |
| Orchestrator, RNA | from the QC manifest | PASS |
| Web, hybrid | DNA + RNA, one combined PCGR report | PASS · both halves |
| Web, worksheet | DNA row | PASS |
| Web, worksheet | lane-split row (R1 = lane 1) | PASS · full depth, all four lanes |
| Web, worksheet | RNA row | PASS |

Running them found eleven defects, all fixed and each pinned by a named test
in `tests/test_regressions.py` — see [TESTING.md](../TESTING.md). The worst:
with several samples present the tumour had been taken to be whichever
sorted first, which put **the normal in the tumour's place** for ordinary
names like PT01-N / PT01-T.

Two behaviours of third-party tools were also found, and are documented
rather than fixed here: SnpEff's first annotation can name a non-clinical
isoform (BRAF V600E appeared as `p.Val640Glu`; the PCGR report is correct),
and PCGR 2.3.2 tiers fusions without a primary-site match (EML4::ALK in lung
came out tier 2). See `RNA.md` and the man page's CAVEATS.
