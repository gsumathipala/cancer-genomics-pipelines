<!-- Created by Brainstorm, 2026. -->
# RNA — what was built, what was not, and why

The RNA branch detects **gene fusions**. That is the whole of it. This
document says what it deliberately does not do, what each omission would
cost to add, and what has to be true before any of it is used clinically.

It follows [CNV_SCOPE.md](CNV_SCOPE.md), and for the same reason: the way a
pipeline acquires a half-finished capability is by somebody adding the easy
half of it without writing down which half that was.

---

## The shape

```
fastq_qc_clean.py ──┬── align_reads.py ── comprehensive_variant_calling.py   (DNA)
   (shared)         │
                    └── fusion_calling.py                                     (RNA)
                         ├── align_rna.py      STAR, chimeric detection on
                         ├── rna_qc_report.py  is a negative interpretable?
                         └── fusion_report.py  what the calls mean
```

Stage 1 and the shared genome are the join. Everything downstream forks,
and `pipeline_orchestrator.py --assay rna` picks the branch.

**Stage 2 is not in the RNA path, and its absence is not an optimisation.**
bwa-mem2 cannot align a read that crosses an exon–exon junction: there is no
contiguous genomic match, so the read is soft-clipped back to one exon or
discarded. Both outcomes are silent — the BAM opens, the mapping rate looks
only somewhat low — and the reads lost are precisely the ones that could
evidence a fusion. The orchestrator forces `--skip-stage2` on `--assay rna`
rather than trusting anyone to remember.

---

## What is built

| | |
|---|---|
| Splice-aware alignment | STAR, with chimeric detection configured explicitly |
| Fusion calling | Arriba |
| Clinical interpretation | PCGR 2.x, which takes fusions as a first-class input — alone, or combined with the DNA library's VCF |
| Library adequacy | mapping rates, junction counts, mitochondrial and rRNA fractions, and a verdict on whether a negative is interpretable |
| Interpretation | fusion calls ranked by clinical salience, same-gene splice events separated out; PCGR tiers them against the tumour site |
| Reports | two standalone HTML/JSON/TSV reports, plus a patient-attributed PDF |
| Profiles | six RNA panel profiles in the same registry as the DNA ones |

---

## What is not built

### Expression quantification

**Not built.** No salmon, no kallisto, no featureCounts, no counts matrix.

*Why not:* expression is only interpretable against a reference
distribution — a cohort of the same tissue, on the same assay, in the same
laboratory. A TPM on its own says nothing, and a TPM compared against a
public cohort sequenced differently says something misleading. It is the
same constraint that blocks copy number in `CNV_SCOPE.md`, for the same
reason, and it is not solved by installing a tool.

*What it would take:* the tool is one conda package. The **reference cohort**
is the work: 20+ samples of the same tumour type on the same assay, or an
explicit decision to report only within-sample ratios and say so.

### RNA variant calling

**Not built.** No SNVs or indels from RNA.

*Why not:* GATK's RNA-seq short-variant path (`SplitNCigarReads` →
HaplotypeCaller) is a **germline** workflow. There is no validated somatic
caller for RNA. Allele fractions from RNA reflect expression as much as
genotype — a variant in a silenced allele is absent from RNA and present in
the genome — so a "VAF" from RNA is not the quantity the DNA branch's
thresholds were built for, and reusing them would produce numbers that look
comparable and are not.

*What it would take:* a clear statement of purpose. Confirming that a DNA
variant is expressed is a defensible, modest goal and needs the DNA call set
as input. *Discovering* somatic variants in RNA is not something this
pipeline should claim.

### UMI consensus calling

**Not built**, on either branch. fastp moves a UMI into the read name; no
consensus or duplex reads are constructed.

*Why not:* it needs fgbio and a grouping/consensus/re-alignment loop, and it
changes what every downstream threshold means.

*Consequence to state out loud:* an RNA kit sold on the strength of its
molecular barcodes will not deliver that strength here.

### STAR-Fusion, and a second caller generally

**Not built.** Arriba only.

*Why not:* STAR-Fusion needs the CTAT genome library — roughly 40 GB on top
of the 30 GB STAR index — against a bundle that already asks for 70 GB.
Arriba reads the STAR BAM this pipeline already produces and ships its own
small reference files inside its conda package.

*What it would take:* the CTAT download, a step in `install_pipeline.py`, and
— the actual work — **a reconciliation policy**. Two callers disagreeing is
the normal case, not the exception, and "run both and intersect" silently
halves sensitivity while "run both and union" doubles the false positives.
The profile schema already accepts `fusion_caller`, so the hook is there.

### Isoform and splice-variant analysis beyond same-gene events

**Not built.** The fusion report labels same-gene events (MET exon 14
skipping, EGFRvIII) because they are clinically important and easy to miss
in a table of gene pairs — but it does not *verify* them. It tells you a
same-gene rearrangement was called and that you must check the exon
boundaries yourself.

### Strandedness inference

**Not built.** Strandedness is recorded from the panel profile, not
measured. Arriba does not need it; anything reporting direction does, and
nothing here does.

---

## What must be true before clinical use

The DNA branch's caveats in the README all still apply. These are the
additional ones, and none of them is optional.

1. **A truth set.** No fusion call from this branch has been checked against
   an independent method. Cell lines with known fusions are the cheap
   starting point — NCI-H2228 (EML4–ALK), HCC78 (SLC34A2–ROS1) — and a
   handful of specimens with FISH or RT-PCR confirmation is the real one.

2. **A limit of detection in your hands.** Fusion sensitivity depends on
   library complexity, read length, and RNA quality, all of which are
   properties of your laboratory rather than of this code.

3. **An RNA quality policy.** DV200 is the best single predictor of whether
   fusion detection can work, it is an instrument measurement of the
   extracted RNA, and **nothing in this pipeline can recover it from the
   sequencing data**. The web form asks for it and the reports state its
   absence. Decide your threshold before you need it.

4. **Orthogonal confirmation.** Fusion callers disagree substantially on
   real data, and an RNA fusion call is evidence of a *transcript*, not
   proof of a genomic rearrangement.

5. **For anchored-PCR panels (Archer, Oncomine): the vendor's caller is the
   validated route.** Arriba and STAR-Fusion were developed on capture and
   whole-transcriptome libraries, and their assumptions about read
   distribution and about what a read-through artefact looks like do not
   hold on AMP data. Results here are orthogonal evidence, not a
   replacement. This is the RNA equivalent of the amplicon duplicate-marking
   trap on the DNA side: the chemistry breaks an assumption the tool does
   not know it is making.

6. **For Ion Torrent RNA: single-end reads and a different error profile.**
   STAR handles single-end data; every paired-end assumption elsewhere in
   this bundle does not, and breakpoint support statistics are weaker
   without mate pairs.

---

## The failure mode this branch is built around

> A degraded FFPE RNA extraction produces a clean, well-formed, **empty**
> fusion table, which is indistinguishable — in the table, in the PDF, in
> the summary line, in every artefact — from a specimen that genuinely
> carries no fusion.

This is the RNA counterpart of the DNA branch's coverage problem, and it is
worse, because it is library-wide rather than regional. Everything about
step 5 follows from it:

- the library QC step runs whenever there is a STAR log to read;
- its verdict is a sentence about *reportability*, not a table of numbers;
- it travels in the run manifest, so the PDF cannot present an
  uninterpretable library as a negative;
- library adequacy is printed **before** the fusion table in every report,
  because a reader who sees the empty table first has already concluded
  "negative" by the time they reach the caveat;
- the things it cannot measure are listed by name, so their absence is not
  read as a pass.

---

## Cost to finish the obvious next things

| Addition | Disk | Real work |
|---|---|---|
| STAR-Fusion as a second caller | ~40 GB | the reconciliation policy, not the install |
| Expression quantification | ~5 GB | a reference cohort of 20+ matched samples |
| UMI consensus (both branches) | none | fgbio loop; re-derive every threshold |
| RNA confirmation of DNA calls | none | define the question precisely first |

Every one of these is a day of plumbing and a season of validation. That
ratio is why they are written down here instead of half-built.
