<!-- Created by Brainstorm, 2026. -->
# Anatomy of a cancer genomics pipeline

*A full dissection: what every tool does, how it does it, why it is there,
and the exact command that runs it — so that someone who has never built
one of these could rebuild this from the ground up.*

This is not the user guide ([README.md](README.md) is). This is the
document I would want if someone handed me a pile of bioinformatics tools
and said "make something a hospital can run."

**Who this is for.** A scientist who codes a little, or a programmer who
has landed in a genomics lab. No prior knowledge of the tools is assumed.
Sections 1–4 are the background; 5–8 are the dissection; 9 onwards is the
hard-won judgement.

---

## Contents

| | |
|---|---|
| **1** | [Why this kind of software is unusual](#1-why-this-kind-of-software-is-unusual) |
| **2** | [The organising idea: make silence impossible](#2-the-organising-idea-make-silence-impossible) |
| **3** | [The map](#3-the-map) |
| **4** | [The file formats, in plain English](#4-the-file-formats-in-plain-english) |
| **5** | [The DNA branch, dissected](#5-the-dna-branch-dissected) |
| **6** | [The RNA branch, dissected](#6-the-rna-branch-dissected) |
| **7** | [PCGR, dissected](#7-pcgr-dissected) |
| **8** | [The four reports written here](#8-the-four-reports-written-here) |
| **9** | [The pearls](#9-the-pearls) |
| **10** | [Architecture patterns worth stealing](#10-architecture-patterns-worth-stealing) |
| **11** | [Rebuilding this from scratch](#11-rebuilding-this-from-scratch) |
| **12** | [A checklist for your own pipeline](#12-a-checklist-for-your-own-pipeline) |

---

## 1. Why this kind of software is unusual

Most software tells you when it fails. A web request 500s, a compiler
errors, a test goes red. Bioinformatics pipelines mostly do not, and that
single fact should shape every decision you make.

Three properties conspire:

**Runs are long.** Aligning one exome takes hours. A mistake in step 2 is
discovered at step 9, tomorrow. Errors must therefore be caught *before*
the expensive part, not wherever they happen to surface.

**There is no ground truth.** When your web app returns the wrong user, the
user notices. When your pipeline returns the wrong variant list, nobody
notices — there is nothing to compare against. A broken run and a correct
run produce the *same shape of output*: a well-formed VCF with plausible
numbers in it.

**Every tool is permissive.** Hand a variant caller the wrong reference and
it will not refuse; it will call variants. Hand it an empty target region
and it will call nothing and exit zero. These tools are libraries, not
guardians. The guarding is your job.

> **The first principle.** The dangerous bug here is not the one that
> crashes. It is the one that produces a confident, well-formatted, wrong
> answer — and the surface area for those is enormous.

---

## 2. The organising idea: make silence impossible

Nearly every design decision in this codebase answers one question:

> *If this goes wrong, how would anyone find out?*

When the honest answer is "they wouldn't", that is the thing to fix — not
by making failure impossible (often you can't) but by making it **loud**.
Three tactics recur:

1. **Check early, where it is cheap.** A PCGR reference bundle without a
   VEP cache is refused during argument parsing, not five hours later.
2. **Say what you did *not* do.** Every report names what it could not
   measure. A missing section reads as "nothing found"; a named absence
   reads as "never measured". Clinically those are opposites.
3. **Record provenance.** Every run writes a manifest: settings, reference
   files, tool versions. Not for you — for whoever reads the result next
   year and must decide whether to trust it.

---

## 3. The map

Two branches share a front end, a genome, and a configuration system.

```
                         ┌──────────────────────┐
     FASTQ files  ──────►│  fastq_qc_clean.py   │  STAGE 1 (shared)
                         │  fastp: trim + QC    │
                         └──────────┬───────────┘
                                    │ run manifest (JSON)
                    ┌───────────────┴────────────────┐
                    ▼                                ▼
        ┌───────────────────────┐        ┌────────────────────────┐
        │   align_reads.py      │        │   fusion_calling.py    │
        │   bwa-mem2 | samtools │        │   ├─ align_rna.py      │
        │   STAGE 2 (DNA only)  │        │   │    STAR            │
        └──────────┬────────────┘        │   ├─ arriba            │
                   │ sorted BAM          │   ├─ rna_qc_report.py  │
                   ▼                     │   └─ fusion_report.py  │
   ┌───────────────────────────────┐     └───────────┬────────────┘
   │ comprehensive_variant_calling │                 │
   │  14 steps: dedup → BQSR →     │                 │
   │  Mutect2 → filter → COSMIC →  │                 │
   │  SnpEff → coverage → PCGR     │                 │
   └───────────────┬───────────────┘                 │
                   │ VCF                             │ fusions.tsv
                   └──────────────┬──────────────────┘
                                  ▼
                          ┌───────────────┐
                          │     PCGR      │  clinical interpretation
                          │  (section 7)  │  takes either — or BOTH,
                          └───────────────┘  for one specimen
```

`pipeline_orchestrator.py` chains it; `--assay dna|rna` picks the branch.
`webapp/` is a Flask front end over the same scripts.

**Why two engines rather than one with a flag?** Because after alignment
*nothing* is shared. Duplicate marking, base recalibration and somatic
calling are either meaningless or actively wrong on RNA. One engine with
`if assay == "rna"` threaded through fourteen steps would be unreadable
and untestable. Separate engines, shared front end, shared config.

---

## 4. The file formats, in plain English

You cannot follow the rest without these. Each is a plain-text format
(usually gzipped) that one tool writes and the next reads.

**FASTQ** — raw reads off the sequencer. Four lines per read: an ID, the
bases, a `+`, and one ASCII character per base encoding its *quality*
(the instrument's confidence). A "paired-end" run gives two files, `R1`
and `R2`, holding the two ends of the same physical DNA fragments, in the
same order.

**SAM / BAM** — reads *after alignment*: the same reads plus where in the
genome each one landed, how well it matched (a CIGAR string like
`100M` = 100 matching bases, or `50M200N50M` = 50 bases, a 200-base gap,
50 more), and flags (is it a duplicate? did its mate align?). BAM is the
compressed binary form. Nearly always sorted by coordinate and indexed
(`.bai`) so tools can seek.

**VCF** — variants. One line per position that differs from the reference:
chromosome, position, reference base(s), alternate base(s), a FILTER
column (`PASS` or why not), an `INFO` column of annotations shared by the
site, and per-sample `FORMAT` fields like depth (`DP`) and allele fraction
(`AF`). Annotation tools work by *adding INFO fields*.

**BED** — regions. Three columns: chromosome, start, end. Used here for
the panel's target regions. **Zero-based and half-open**: `chr1 100 200`
means bases 101–200 in one-based counting. Off-by-one errors between BED
and VCF conventions are a classic bug.

**GTF** — gene annotation. Where every gene, transcript and exon is. The
RNA branch needs it to know where junctions *should* be, and to turn a
breakpoint into a gene name.

> **A trap that will cost you a day.** Two conventions exist for naming
> chromosomes: UCSC (`chr1`, `chrM`) and Ensembl (`1`, `MT`). Files in
> different conventions do not error when combined — they simply match
> nothing. Section 9.I covers this properly.

---

## 5. The DNA branch, dissected

Fourteen steps. For each: **what** it does, **how** it works, the **exact
command** this pipeline runs, and **what goes wrong silently**.

Commands below are real, captured from `--dry-run`. `WORK` stands for the
run directory.

---

### Step 1 — QC and trimming (fastp)

**What.** Removes sequence that is not biology: adapters, low-quality tails,
poly-G runs, reads too short to align.

**How.** fastp streams both FASTQ files together. For adapter removal in
paired-end mode it does not need to be told the adapter sequence: it
*overlaps R1 and R2*: where two mates overlap, anything past the overlap is
adapter by definition. It is a single pass in C++, so it is fast.

```bash
fastp -i S_R1.fastq.gz -I S_R2.fastq.gz \
      -o S_R1.clean.fastq.gz -O S_R2.clean.fastq.gz \
      --thread 8 \
      --length_required 50 \
      --qualified_quality_phred 15 \
      --unqualified_percent_limit 40 \
      --n_base_limit 5 \
      --average_qual 20 \
      --detect_adapter_for_pe \
      --trim_poly_g --poly_g_min_len 10 \
      --overrepresentation_analysis \
      --html S.fastp.html --json S.fastp.json
```

| Flag | Why |
|---|---|
| `--detect_adapter_for_pe` | infer adapters from mate overlap rather than a hard-coded list |
| `--trim_poly_g` | **NextSeq/NovaSeq specific.** Two-colour chemistry encodes "no signal" as G, so a dying cluster emits a run of Gs that is not sequence |
| `--length_required 50` | shorter reads align ambiguously. Lowered to 35 on RNA — FFPE RNA is fragmented |
| `--qualified_quality_phred 15` + `--unqualified_percent_limit 40` | discard a read if >40% of bases are below Q15 |
| `--average_qual 20` | drop reads whose mean quality is poor |
| `--json` | machine-readable metrics; the pipeline parses this into its manifest |

**Deliberately NOT used:** `--correction` (overlap-based base correction).
It rewrites bases where mates disagree — exactly the evidence a low-allele-
fraction variant or a UMI workflow depends on.

**Silent failure.** Adapter left in aligns somewhere and looks like
sequence. Poly-G untrimmed produces spurious alignments in G-rich regions.

---

### Step 2 — Alignment (bwa-mem2)

**What.** Finds, for every read, the position in the 3.1-billion-base genome
it most likely came from.

**How.** bwa-mem2 uses an **FM-index** (a compressed full-text index built
on the Burrows–Wheeler Transform) to find exact-match seeds in
near-constant time, extends them with banded Smith–Waterman alignment, and
scores the result. It is a rewrite of BWA-MEM with vectorised (AVX2/512)
inner loops — same output, roughly 2–3× faster, at the cost of a larger
index (~20 GB for hg38, ~1–2 hours to build, **once**).

```bash
bwa-mem2 mem -t 8 \
  -R '@RG\tID:S\tSM:S\tPL:ILLUMINA\tLB:S_lib\tPU:unknown' \
  hg38.fa S_R1.clean.fastq.gz S_R2.clean.fastq.gz \
| samtools sort -@ 7 -o S.sorted.bam -
samtools index S.sorted.bam
```

**The `-R` read group is not optional.** GATK refuses to run without one.
`SM:` (sample) is the field every downstream tool uses to know whose reads
these are; get it wrong and you will label one patient's variants with
another's name.

**Why piped.** The intermediate SAM for an exome is hundreds of gigabytes
of text. Piping straight into `samtools sort` means it never touches disk.

> **Implementation pearl.** When you pipe like this from Python, you must
> **drain the aligner's stderr from a separate thread**. bwa streams
> progress to stderr; if you only read it after the process exits, a long
> alignment fills the 64 KB OS pipe buffer, bwa blocks writing to it, and
> the whole pipeline deadlocks — hours in, with no error. This pipeline
> learned that the hard way; see `align_sample()`.

---

### Step 3 — Duplicate marking (GATK MarkDuplicatesSpark)

**What.** Identifies reads that are PCR copies of one original molecule.

**How.** Library prep amplifies each fragment. Copies of one molecule land
at *identical* start and end coordinates. MarkDuplicates groups reads by
those coordinates (plus strand and mate position), keeps the
highest-quality representative and sets bit `0x400` on the rest.

```bash
gatk MarkDuplicatesSpark \
  -I S.sorted.bam -O S.dedup.bam -M S.dup_metrics.txt \
  --tmp-dir WORK/tmp --spark-master 'local[8]' --verbosity ERROR
samtools index S.dedup.bam
```

**Why it matters.** Twenty reads from one molecule are *one* observation.
Counted as twenty, a single PCR error becomes a confident variant at 100%
allele fraction.

**Marked, not removed** — duplicates still carry usable base-quality
information, and downstream tools honour the flag.

> `--spark-master local[N]` is the knob that actually parallelises this.
> An earlier version passed `--conf spark.local.cores`, which is not a
> real Spark property — Spark ignores unknown keys silently, so `--threads`
> did nothing here for months.

**Silent failure — and the big one:** on **amplicon** data this step is
catastrophic. See section 9.II.

---

### Step 4 — Base Quality Score Recalibration (BQSR)

**What.** The instrument's quality scores are systematically biased —
by machine cycle, by the preceding base, by context. BQSR measures that
bias empirically and rewrites the scores.

**How.** Two passes. `BaseRecalibrator` walks the BAM, and at every
position *not* in a database of known variants assumes any mismatch is an
error. Tabulating those by covariate (reported quality × cycle × dinucleotide)
gives the true error rate per bucket. `ApplyBQSR` then rewrites each base's
quality accordingly.

```bash
gatk BaseRecalibrator -R hg38.fa -I S.dedup.bam \
  -O S.recal --known-sites dbsnp.vcf.gz \
  -L panel.bed --interval-padding 100

gatk ApplyBQSR -R hg38.fa -I S.dedup.bam \
  -O S.bqsr.bam --bqsr-recal-file S.recal
```

**Why `--known-sites` is mandatory.** Without it, real germline variants —
about 3–4 million per genome — are counted as sequencing errors, and BQSR
"corrects" quality scores downward for being right.

> **Note the asymmetry:** `-L panel.bed` restricts the *model-building*
> pass but **not** `ApplyBQSR`. The model is fit on target reads, where the
> data is; but the recalibration is applied to the whole BAM, because
> rewriting only part of it would leave two quality scales in one file.

**Skipped automatically** when no known-sites resource is supplied, and
worth skipping deliberately on panels under ~1 Mb: too few covered sites
to fit a model, so it fits noise.

---

### Step 5 — Somatic variant calling (Mutect2)

**What.** Proposes positions where this tumour differs from the reference,
somatically.

**How.** Mutect2 is a **local haplotype assembler**, not a
position-by-position caller. In each "active region" showing evidence of
variation it assembles the reads into a De Bruijn graph, extracts candidate
haplotypes, realigns every read to each haplotype with a pair-HMM, and
computes the likelihood that the data came from a model containing the
variant versus one without. That is why it handles indels and clustered
variants far better than naive pileup callers.

```bash
gatk Mutect2 -R hg38.fa -I S.bqsr.bam -tumor S \
  -O S.mutect2.vcf.gz \
  --f1r2-tar-gz S.f1r2.tar.gz \
  --germline-resource af-only-gnomad.hg38.vcf.gz \
  -pon 1000g_pon.hg38.vcf.gz \
  -L panel.bed --interval-padding 100 \
  --tmp-dir WORK/tmp --verbosity ERROR
```

| Flag | What it buys |
|---|---|
| `--germline-resource` | gnomAD allele frequencies. In **tumour-only** mode this is the main defence against reporting inherited variants as somatic |
| `-pon` | panel of normals: positions recurrently mutated across normal samples are recurrent *artefacts* of this assay |
| `--f1r2-tar-gz` | raw counts for the read-orientation model (step 6) |
| `-L` + `--interval-padding 100` | restrict to the panel. **Two orders of magnitude of noise** otherwise — see 9.V |

**Mutect2 is deliberately permissive.** It proposes; `FilterMutectCalls`
disposes. Do not read its raw output as a result set.

---

### Step 6 — Read-orientation model (FFPE artefacts)

**What.** Builds a model of orientation-biased artefacts.

**How.** Formalin fixation deaminates cytosine to uracil, read as thymine —
producing **C>T / G>A** artefacts. Crucially these appear on only one
strand of the original duplex, so they show up in reads of one orientation
(F1R2 vs F2R1) and not the other. Real mutations are present on both.
`LearnReadOrientationModel` fits that asymmetry from the counts step 5
collected.

```bash
gatk LearnReadOrientationModel \
  -I S.f1r2.tar.gz -O read_orientation_model.tar.gz
```

**Keep this step on FFPE material.** It is what separates fixation
chemistry from biology. Residual damage still reaches PASS at low allele
fraction, which is why `--min-allele-fraction` earns its place too.

---

### Step 7 — Contamination estimation

**What.** Measures how much of this "tumour" sample is actually someone
else's DNA.

**How.** At common biallelic SNP sites, a pure sample shows allele
fractions near 0, 0.5 or 1. Contamination pulls those toward the
population allele frequency. `GetPileupSummaries` tabulates the counts;
`CalculateContamination` fits the contamination fraction from the
departure, and segments the genome so that copy-number changes are not
mistaken for contamination.

```bash
gatk GetPileupSummaries -I S.bqsr.bam \
  -V small_exac_common_3.hg38.vcf.gz \
  -L small_exac_common_3.hg38.vcf.gz \
  -O S.pileups.table

gatk CalculateContamination -I S.pileups.table \
  -O S.contamination.table --tumor-segmentation S.segments.table
```

**Without it, filtering assumes a contamination fraction of zero** rather
than measuring one. This pipeline also refuses to trust an estimate built
on too few informative sites (panels often have few common SNPs), and says
so rather than reporting a confident number from thin data.

---

### Step 8 — Microsatellite instability (MSIsensor2)

**What.** Scores whether the tumour has lost DNA mismatch repair — a
predictive marker for immunotherapy.

**How.** Microsatellites are short repeats (e.g. `ACACAC…`). Without
mismatch repair their lengths become unstable. MSIsensor2 compares the
observed length distribution at each microsatellite against a **trained
model** of what a stable sample looks like — which is what lets it work
*tumour-only*, with no matched normal.

```bash
msisensor2 msi -M models_hg38 -t S.bqsr.bam -o S.msi -b 8 -e panel.bed
```

> Two traps. The conda package ships the **binary only** — the models are a
> separate download, and without them tumour-only MSI cannot run at all.
> And PCGR will *not* fill this gap: it restricts its own MSI prediction to
> WGS/WES tumour–normal runs and omits the section without comment.
> Upstream calls MSI-H at **≥20%** for tumour-only data; the 3.5% figure
> belongs to the older paired caller and would call almost anything
> unstable.

---

### Step 9 — Filtering, then two local filters

**9a. FilterMutectCalls.** Applies a single unified probabilistic model
across every error mode — contamination, strand-orientation bias, germline
risk, mapping quality, base quality, fragment length, clustered events —
and chooses the threshold that maximises expected utility (the "F-score
optimal" threshold) rather than applying fixed cut-offs.

```bash
gatk FilterMutectCalls -V S.mutect2.vcf.gz -O S.filtered.vcf.gz \
  -R hg38.fa \
  --ob-priors read_orientation_model.tar.gz \
  --contamination-table S.contamination.table \
  --tumor-segmentation S.segments.table \
  --min-allele-fraction 0.05
```

Calls are **tagged in the FILTER column**, not deleted. `PASS` means it
survived everything.

**9b. Depth floor.** GATK has no depth filter of its own and will PASS a
variant supported by two reads.

```bash
bcftools filter -s low_depth -m + -e 'FORMAT/DP[0] < 50' \
  -O z -o S.depthfloor.vcf.gz S.filtered.vcf.gz
```

`-s low_depth` names the tag, `-m +` *adds* to any existing FILTER rather
than replacing it. Tagged, never removed.

**9c. Foldback insertion tagger.** A custom filter for an artefact no
generic tool catches — see 9.V. It reads each insertion, reverse-complements
it, and compares against the flanking reference; the matched fraction is
recorded in `INFO/FBMATCH` for *every* insertion tested, so a reviewer sees
the evidence rather than trusting a threshold.

---

### Step 10 — COSMIC annotation

**What.** Adds known-somatic-mutation identifiers (`COSV…`).

```bash
bcftools annotate -a Cosmic.chr.vcf.gz -c ID,INFO \
  -O z -o S.cosmic.vcf.gz S.filtered.vcf.gz
```

`-c ID,INFO` says which columns to transfer.

> **The single most instructive failure in this whole pipeline.** COSMIC
> ships in Ensembl contig naming (`1`, `MT`); the installed hg38 is UCSC
> (`chr1`, `chrM`). `bcftools annotate` with mismatched contigs matches
> **nothing at all**, exits 0, and the run reports success with an
> annotation step that did nothing. The installer therefore renames
> COSMIC's contigs at install time, and the pipeline compares contig sets
> and warns before annotating.

---

### Step 11 — Gene annotation (SnpEff)

**What.** Turns coordinates into biology: which gene, which transcript,
what consequence (missense? frameshift? splice site?).

**How.** SnpEff holds a pre-built database of transcript models, finds the
transcripts overlapping each variant, and predicts the effect on the
protein, ranking by severity.

```bash
snpEff -Xmx8g hg38 -noStats -noLog S.cosmic.vcf.gz > S.annotated.vcf
```

`-Xmx8g` is **not optional**: the bioconda launcher defaults to 1 GB, which
is too small for hg38 and dies partway through with an `OutOfMemoryError`.
`-noStats` suppresses the summary HTML (PCGR is the reporting layer here).

---

### Step 12 — Target coverage

Custom (`coverage_report.py`) — see section 8. It answers the question no
VCF can: *which parts of the panel were sequenced deeply enough that an
absent variant means something.*

```bash
samtools depth -a -b panel.bed -Q 20 -q 20 S.bqsr.bam
```

`-a` reports **every** position including zero-depth ones — without it,
uncovered bases are simply missing from the output, which is precisely the
silence this step exists to break. `-Q`/`-q` set the mapping- and
base-quality floors so the measurement matches what the caller counted.

---

### Steps 13–14 — Clinical report and summary

Step 13 is PCGR (section 7). Step 14 writes the run manifest: every
parameter, every reference path, every tool version, which steps completed.

---

## 6. The RNA branch, dissected

Six steps. The differences from DNA are the lesson.

---

### Step 0 — Building the STAR index (once)

```bash
STAR --runMode genomeGenerate --runThreadN 16 \
     --genomeDir star_hg38_150 \
     --genomeFastaFiles hg38.fa \
     --sjdbGTFfile gencode.v50.primary_assembly.annotation.gtf \
     --sjdbOverhang 149
```

~1 hour, ~32 GB RAM, ~30 GB disk. **Three things are baked in**: the
genome, the annotation, and the read length.

`--sjdbOverhang` is read length − 1: the amount of sequence STAR stores on
each side of every known junction so a read can be matched across it.
Built for 100 bp and used on 150 bp reads, STAR **runs perfectly well and
quietly loses junction sensitivity**. So `align_rna.py` writes a build
record (genome, GTF, read length, STAR version) into the index directory
and later runs compare against it.

---

### Step 1 — QC (fastp)

Same tool, one important difference: `--length_required 35` rather than 50.
FFPE RNA is fragmented, and a floor inherited from DNA discards much of a
usable library.

---

### Step 2 — Splice-aware alignment (STAR)

**Why bwa-mem2 cannot be used.** A read crossing an exon–exon junction has
**no contiguous match anywhere in the genome** — the intron is missing from
the transcript. bwa-mem2 has no way to represent that, so it soft-clips the
read back to one exon or discards it. Silently. Those are exactly the reads
that could evidence a fusion.

**How STAR works.** Two phases. *Seed searching*: it finds Maximal Mappable
Prefixes using an uncompressed suffix array — the longest prefix of the read
that maps somewhere; then it restarts from the next unmapped base. A read
spanning a junction naturally produces two seeds in two exons. *Stitching*:
seeds are clustered and stitched with a gap, scoring known junctions (from
the GTF) higher than novel ones. The suffix array is why it needs ~30 GB of
RAM and is extremely fast.

```bash
STAR --runThreadN 8 --genomeDir star_hg38_150 \
  --genomeLoad NoSharedMemory \
  --readFilesIn S_R1.clean.fastq.gz S_R2.clean.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix WORK/rna_aligned/S. \
  --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
  --outFilterMultimapNmax 50 \
  --peOverlapNbasesMin 10 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM SoftClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50
```

**The chimeric block is the whole point, and none of it is default.**

| Flag | Why it is there |
|---|---|
| `--chimSegmentMin 10` | **STAR's default is 0, meaning chimeric detection is OFF.** This single value is the difference between a fusion panel and an empty table |
| `--chimOutType WithinBAM SoftClip` | write chimeric alignments into the main BAM so the caller reads one input |
| `--chimScoreJunctionNonGTAG 0` | a fusion breakpoint is a **genomic rearrangement, not a splice site**, so it must not be penalised for lacking the canonical GT/AG intron motif. STAR's default penalty is −1 |
| `--chimScoreDropMax 30` | tolerate the alignment-score penalty a genuine fusion incurs |
| `--chimMultimapNmax 50`, `--outFilterMultimapNmax 50` | fusion partners are often repetitive; discarding multimappers discards the fusion |
| `--peOverlapNbasesMin 10` | use mate overlap — on the short fragments FFPE RNA yields, that is most of the library |
| `--genomeLoad NoSharedMemory` | a shared genome segment surviving a cancelled run is a 30 GB leak needing manual cleanup |
| `--outBAMcompression 0` | the caller reads this once and it is deleted; compression is pure cost |

---

### Step 3 — Sort and index (for QC and IGV only)

```bash
samtools sort -@ 7 -o S.sorted.bam S.Aligned.out.bam
samtools index -@ 7 S.sorted.bam
```

> **The fusion caller needs the UNSORTED BAM.** Arriba reads mate pairs
> together, and coordinate sorting separates them. Handing it the sorted
> file is a supported command that finds far fewer fusions. This pipeline
> gives Arriba STAR's own output and produces the sorted copy separately.

---

### Step 4 — Fusion calling (Arriba)

**How.** Arriba reads STAR's chimeric alignments and applies a cascade of
filters encoding what artefacts look like: read-through transcription
between neighbouring genes, homologous/paralogous partners, low-complexity
regions, PCR duplicates, mismapped reads at repeats, and a **blacklist** of
recurrent artefacts observed across thousands of normal samples. What
survives is scored `high`, `medium` or `low` confidence based on breakpoint
support, whether breakpoints fall at annotated exon boundaries, and whether
the fusion is in frame.

```bash
arriba -x S.Aligned.out.bam \
       -o S.fusions.tsv -O S.fusions.discarded.tsv \
       -a hg38.fa -g gencode.gtf \
       -b blacklist_hg38.tsv.gz \
       -k known_fusions_hg38.tsv.gz \
       -t known_fusions_hg38.tsv.gz \
       -p protein_domains_hg38.gff3
```

| Flag | Consequence if omitted |
|---|---|
| `-b` blacklist | **recurrent read-through artefacts are reported as high-confidence fusions** and the table fills with noise |
| `-k` known fusions | a real low-support `EML4–ALK` loses its sensitivity boost and can drop below the bar |
| `-t` tags | known fusions are not labelled as such |
| `-p` protein domains | cosmetic: retained domains are not annotated |
| `-O` discarded | you lose the ability to ask *why* an expected fusion is absent |

These reference files **ship inside the conda package** (usually
`$CONDA_PREFIX/var/lib/arriba`), which is why there is no separate download
step — and why this pipeline locates them explicitly and warns loudly if it
cannot.

**No duplicate marking, and no flag to enable it.** On RNA, twenty
identical reads from a highly expressed gene are twenty genuine
observations of an abundant transcript. There is no correct value for that
setting, so the setting does not exist.

---

### Steps 5–6 — Library adequacy, then the fusion report

Both custom; section 8. Step 5 exists because of the single most dangerous
property of this branch: **a degraded FFPE RNA library produces a clean,
well-formed, empty fusion table, indistinguishable from a true negative.**

---

## 7. PCGR, dissected

PCGR (Personal Cancer Genome Reporter) is the clinical interpretation
layer, and it is the most complex single component here. It is *not* a
variant caller — it takes finished calls and answers **"which of these
matter, and why?"**

It runs in its **own conda environment** (plus a sibling `pcgrr` for the R
reporting side) and cannot run inside the pipeline's environment. This
pipeline therefore always invokes it as a separate subprocess.

### 7.1 What it accepts

PCGR 2.x requires **at least one molecular input**, and a somatic VCF is
only one of the options:

| Input | Flag |
|---|---|
| somatic SNVs/indels | `--input_vcf` |
| copy-number segments | `--input_cna` |
| **RNA fusions** | `--input_rna_fusion` |
| RNA expression (TPM) | `--input_rna_expression` |

A VEP cache is required **only if a VCF is given** — so a fusion-only
report needs none, there being nothing to annotate. That is what lets this
pipeline's RNA branch produce a real tiered interpretation, and what lets a
DNA+RNA specimen be reported as **one document**.

### 7.2 Its internal stages

Taken from the running code, in order:

1. **Verify arguments** — configure the workflow, check the bundle matches
   the installed version.
2. **Input validation** — for fusions: required columns `FusionGene`,
   `LeftBreakpoint`, `RightBreakpoint`, `SplitReads`; gene pairs separated
   by `--` or `::`, exactly two parts, breakpoints as `<chrom>:<pos>`.
   **An empty fusion file is rejected outright**, so "no fusions found"
   must produce *no file*, not an empty one.
3. **VEP** — the heavy annotation pass (below).
4. **vcfanno** — overlays the reference databases onto the VCF.
5. **pcgr-summarise** — variant and cancer-gene annotation; oncogenicity
   scoring; mutational hotspot matching.
6. **vcf2maf** — conversion for downstream compatibility.
7. **CNA / RNA-fusion / RNA-expression sections** — each annotated against
   the biomarker databases, or announced as *omitted* when absent.
8. **Two-hit analysis** — tumour-suppressor genes with both a somatic
   variant and loss of the other allele.
9. **OncoKB annotation** — only with an API token; otherwise explicitly
   skipped and NA columns added.
10. **Report generation** — a Quarto/R interactive HTML report, plus TSV
    and XLSX.

### 7.3 The VEP pass, in detail

PCGR does not shell out to VEP casually — it pins a long, deliberate flag
set:

```
vep --input_file in.vcf --output_file out.vcf \
    --dir <vep_cache> --assembly GRCh38 --cache_version <pinned> \
    --fasta hg38.fa --offline --cache --species homo_sapiens \
    --hgvs --af_gnomade --af_gnomadg --max_af \
    --variant_class --domains --symbol --protein --ccds --mane \
    --uniprot --appris --biotype --tsl --canonical \
    --numbers --total_length --allele_number \
    --failed 1 --no_stats --no_escape --xref_refseq --vcf \
    --check_ref --dont_skip --flag_pick_allele_gene \
    --pick_order <configured> --buffer_size <n> --fork <n> \
    --plugin NearestExonJB,max_range=50000 \
    --plugin MaxEntScan,<dir> --plugin NMD
```

The ones worth understanding:

- **`--check_ref` / `--dont_skip`** — verify that the REF base in the VCF
  actually matches the reference genome, and *do not silently drop* the
  ones that don't. This is how a wrong-reference VCF gets caught rather
  than quietly half-annotated.
- **`--flag_pick_allele_gene` + `--pick_order`** — a variant hits many
  transcripts. Rather than discard all but one, VEP *flags* one per
  allele-and-gene by a configured priority (MANE, canonical, biotype,
  length…). Downstream, PCGR reads the flagged block.
- **`--mane`, `--canonical`, `--appris`, `--tsl`** — transcript-quality
  annotations, so the choice above is defensible rather than arbitrary.
- **`--af_gnomade --af_gnomadg --max_af`** — population frequencies, the
  main evidence that something is inherited rather than somatic.
- **`--plugin NMD`** — will this premature stop codon actually trigger
  nonsense-mediated decay? A truncating variant in the last exon often
  does not, and that changes its interpretation.
- **`--plugin MaxEntScan`** — maximum-entropy splice-site scoring, for
  variants near exon boundaries.
- **`--offline`** — no network calls; reproducible and fast.

**Why VEP is slow and hungry.** The cache is ~24 GB because it holds every
transcript model, regulatory feature and frequency for the whole genome.
`--fork` parallelises; `--buffer_size` trades memory for throughput.

### 7.4 The knowledge bases it overlays

From the installed bundle:

| Group | Contents |
|---|---|
| **biomarker** | CGI and CiVIC — clinical, literature and variant-level actionability evidence |
| **variant** | ClinVar, dbNSFP (ensemble *in silico* predictors), gnomAD non-cancer, GWAS catalogue, panel of normals, TCGA frequencies, dbMTS |
| **gene** | cancer-gene census, transcript cross-references, virtual gene panels |
| **misc** | mutational **hotspots**, mutational signatures, protein domains, cytobands, TMB reference distributions, clinical trials, Grantham scores, splice databases |
| **drug** | drug–target relationships |

### 7.5 How it decides "oncogenic"

PCGR implements the **ClinGen/CGC/VICC** standard for oncogenicity
classification: a points system. Each satisfied criterion contributes a
score, positive (oncogenic pole) or negative (benign pole) — for example
`ONCG_OVS1_A` (very strong evidence, e.g. a null variant in a
tumour-suppressor) scores **+8**, while `ONCG_SBS2_A` (strong benign
evidence) scores **−4**. Evidence includes population frequency, *in
silico* predictor consensus (a majority of callers must agree), hotspot
membership, functional-domain location and known biomarker matches.

The total is then banded against fixed thresholds:

| Total score | Classification |
|---|---|
| ≥ 9 | Oncogenic |
| 5 – 8 | Likely oncogenic |
| −1 … 4 | Variant of uncertain significance |
| −6 … −2 | Likely benign |
| ≤ −7 | Benign |

Separately, **actionability tiers** (AMP/ASCO/CAP) are assigned from the
biomarker evidence *for the specified tumour site* — which is exactly why
`--tumor_site` matters: the same variant can be predictive in one tissue
and merely oncogenic in another.

### 7.6 The two settings that quietly ruin a report

**TMB is off by default and its denominator is a guess.** PCGR computes no
TMB unless `--estimate_tmb` is passed, and says nothing when it is off —
the report simply lacks the section, which reads as "nothing found". Worse,
when it *is* on, `--effective_target_size_mb` defaults to an exome-sized
**34 Mb**. On a 4.5 Mb panel that is ~7.5× too low; on a 2 Mb panel ~17×.
See 9.III.

**Depth and allele fraction live in FORMAT, not INFO.** PCGR filters on
INFO tags (`--tumor_dp_tag`, `--tumor_af_tag`), but Mutect2 writes DP/AF as
per-sample FORMAT fields. Unless they are lifted across, PCGR has nothing
to filter on and computes TMB over *unfiltered* calls. This pipeline does
that lift explicitly (`--lift-tags`, writing `TDP`/`TAF`) into a *new* file
— the called VCF is never modified.

> **A third trap this pipeline works around:** PCGR aborts if the input
> already carries an INFO tag it intends to write, and step 10's
> `bcftools annotate -c ID,INFO` copies COSMIC's INFO wholesale (`STRAND`
> being the live collision). So reserved tags are stripped into a separate
> file first, and the original is used unchanged when nothing collides.

---

## 8. The four reports written here

Nothing existing answered these questions, so they are written in this
bundle — all four self-contained HTML plus JSON and TSV, with **no external
CSS or JavaScript**, because they get emailed, archived and opened on
machines that are not this one.

### 8.1 `coverage_report.py` — what a negative is entitled to say

Runs `samtools depth -a -b panel.bed -Q 20 -q 20`, walks every target
region, and counts bases below the depth a call would need. Output names
every region that fell short.

The `-a` is load-bearing: without it, zero-depth positions are simply
absent from the output — the exact silence the report exists to break.

### 8.2 `rna_qc_report.py` — is an empty fusion table a negative?

Deliberately cheap: it needs only samtools and files the run already
produced.

| Source | Measurement |
|---|---|
| STAR `Log.final.out` | input reads, uniquely mapped %, **reads lost as "too short"** (the degraded-library signature), splice junctions, chimeric reads |
| `samtools idxstats` | mitochondrial fraction (reads the BAM *index*, so it is free) |
| the GTF + `samtools view -c -L` | rRNA fraction, from intervals derived on the fly by grepping `gene_type "rRNA"` |

It ends not with a table but a **sentence about reportability**:

> *3 of 6 checks did not clear their bar. An empty fusion table from this
> library should NOT be reported as a negative without explaining these.*

And it lists what it **cannot** measure, by name:

- **duplication rate** — not a quality metric on RNA; an abundant
  transcript produces identical reads legitimately, so a good amplicon RNA
  library approaches 100% "duplication". Reporting it invites someone to
  act on it.
- **DV200 / RIN** — an instrument measurement of the extracted RNA. It
  cannot be recovered from FASTQ, and it is the best single predictor of
  whether fusion detection could work at all. So it is a field the operator
  records, not a number this computes.

### 8.3 `fusion_report.py` — ranking by what matters clinically

Arriba ranks by *evidence quality*. A high-confidence read-through between
two housekeeping genes therefore sorts above a medium-confidence
`EML4–ALK`. This re-sorts by **clinical salience**:

1. calls involving a gene with established fusion-directed therapy or
   diagnostic significance,
2. **same-gene splice events**,
3. caller confidence,
4. read support.

Point 2 is easy to miss and clinically expensive: **MET exon 14 skipping**
is one of the most actionable findings a lung panel produces, it is a
splice event *within one gene* rather than a fusion between two, and in a
table sorted by gene pair it reads as `MET--MET`.

Calls below the confidence bar are **marked, never dropped** — a laboratory
investigating a clinically expected fusion needs to see that it was called
weakly, not be told nothing was found.

It also writes `<sample>.pcgr_fusions.tsv` in PCGR's schema, converting
Arriba's columns:

| PCGR column | From Arriba |
|---|---|
| `FusionGene` | `gene1` + `--` + `gene2`, first symbol of each |
| `LeftBreakpoint` / `RightBreakpoint` | `breakpoint1` / `breakpoint2` |
| `SplitReads` | `split_reads1` + `split_reads2` |

Two judgement calls, both recorded in the code: discordant mates are **not**
added (they support the pairing without crossing the junction, and PCGR's
own threshold is written against split reads), and `Score` is **omitted**
rather than fabricated from high/medium/low.

### 8.4 `webapp/pdfwriter.py` — a dependency-free PDF writer

A PDF is a header, numbered objects, a cross-reference table of byte
offsets, and a trailer. This writes that directly using the standard-14
fonts (no embedding), in a few hundred lines with `import datetime` as its
only import.

**Why not reportlab/WeasyPrint/wkhtmltopdf?** The analysis scripts are
standard-library-only by design, and a clinical report generator is a bad
place to inherit cairo/pango or an external binary that is missing on
cluster nodes. The trade is stated: Latin-1 only, no images — which is why
`report.py` flags characters in a patient name it cannot render rather than
mangling them.

---

## 9. The pearls

Each cost somebody a day. Grouped so the *categories* generalise even where
the specific tools do not.

### I. Reference data is where runs die

**Contig naming will get you at least once.** UCSC writes `chr1`, `chrM`;
Ensembl writes `1`, `MT`. Mixing them does not always error:

- COSMIC in Ensembl naming against a UCSC genome annotates **nothing at
  all** — and the run reports success.
- A panel BED in the wrong naming makes GATK abort (loud, fine) but makes
  the *coverage report* find no reads in any region — which reads as "none
  of the panel was covered". A confident, completely wrong clinical
  statement.

> **Pearl.** Detect the convention, translate, and say so. Never edit the
> user's file — write a corrected copy and name it in the log.

**Verify downloads by size, not existence.** A truncated 20 GB reference
passes `os.path.exists()` forever, and every later check that only asks
"is it there?" will keep passing. Compare against the server's
`Content-Length` before renaming into place.

**Version-lock what is baked in.** The STAR index contains a genome, an
annotation and a read length. Record all three at build time and compare
later; an index outliving its GTF is undetectable at run time.

**Dependency versions can fail in the middle.** GATK's Spark tools do not
run on JDK 24+, and the failure lands at duplicate marking — *after*
alignment has burned two hours. `gatk --version` passes either way. So the
installer parses `java -version` and enforces the range up front.

**Sibling-path assumptions.** PCGR locates its R companion as
`dirname($CONDA_PREFIX)/pcgrr`. Install them anywhere else and it fails
late. Worth knowing that tools *do* make assumptions like this.

### II. The chemistry decides the algorithm

The thing that most surprises programmers: **the correct code depends on
how the DNA was prepared in the lab.**

**Amplicon panels must not have duplicates marked.** Amplicon reads begin
and end at primer coordinates *by construction*, so every read from one
amplicon looks like a duplicate of every other. MarkDuplicates flags nearly
the whole library; the caller then works from a handful of surviving reads.
A catastrophic, silent loss of depth whose only symptom is an implausibly
thin VCF.

**BQSR needs data.** It fits an empirical error model from covered sites.
On a 20 kb hotspot panel there are nowhere near enough, so it fits noise
and applies it to real quality scores.

**Padding depends on chemistry.** Capture probes pull in fragments that
extend past the target, so a variant at the first base of an exon needs
~100 bp of flank to be callable. On an amplicon panel the target *is* the
amplicon, and padding walks into primer sequence — bases that come from the
oligo, not the patient.

**Ion Torrent is not Illumina.** Single-end reads, and a different indel
error profile in homopolymers. Every paired-end assumption in a pipeline
tuned on Illumina data silently does not hold.

> **Pearl.** If a setting's correct value depends on the wet lab, do not
> give it a default and hope. Name the assay and derive the settings from
> that — which is what `panel_profiles.py` does, with the precedence rule
> *explicit flag > profile > default*, printed at the start of every run.

### III. Every denominator is a lie until you check it

**TMB is the worked example.** Tumour mutational burden is mutations per
megabase — a count divided by the size of the region sequenced. PCGR
assumes **34 Mb** for a "targeted" assay because that is exome-sized.

Run a 4.5 Mb panel: TMB ~7.5× too low. A 2 Mb panel: ~17× too low. **PCGR
does not warn.** It prints a number.

The fix: measure the footprint from the BED actually in use — **merging
overlapping regions first**, since vendor BEDs routinely overlap their
probe tiles and summing them unmerged inflates the denominator, deflating
TMB in the same silent direction.

**And know when a number is not interpretable at all.** Below ~1 Mb a
single extra artefact moves TMB by several mutations per Mb. This pipeline
switches TMB off below that and says why, rather than printing a figure.

> **Pearl.** Find every place your system divides by something and ask
> where that denominator came from. In science code it is usually a default
> that was right for the author's data.

### IV. Absence of evidence is not evidence of absence

**The most important idea in the whole codebase.**

A variant caller reports variants. It does not distinguish "I looked here
and found nothing" from "I never looked here". Both appear as *nothing at
all*.

So a region the sequencing never reached and a region that is genuinely
wild type are **identical in the VCF, identical in the report, and
identical to the clinician**. One means "no mutation"; the other means "no
information". Opposite conclusions.

The RNA branch has the same hole and worse: a degraded FFPE extraction
produces a clean, well-formed, **empty** fusion table.

> **Pearl.** In any system that reports findings, ask what an empty result
> means. If "we found nothing" and "we couldn't look" render identically,
> you have a serious problem — and it will not show up in testing.

**The corollary shapes the whole interface:** library adequacy is printed
**before** the findings everywhere — in the HTML reports, in the PDF, in
the ordering of links on the run page. A reader who sees an empty table
first has already concluded "negative" before reaching the caveat.

### V. Artefacts that look like biology

**Foldback insertions.** During library prep a fragment end can fold back
on itself and be extended, producing a sequence followed by its own reverse
complement. Aligned, that looks like an *insertion*.

On the panel this was built for:

- **95.6%** of insertions >10 bp were reverse-complement copies of the
  adjacent reference — against **0.2%** in the forward direction
- 1,882 of 2,047 insertions were this artefact
- they accounted for **206 of 270 PASS frameshift/stop calls**, including
  BRCA1, BRCA2, PALB2, PTEN and ATM frameshifts that PCGR tiered as
  oncogenic

No generic filter catches them: they are not low-depth, not strand-biased,
not in tandem repeats, and they reach allele fractions above 0.3.

> **Pearl.** When a filter is doing nothing, ask whether the artefact you
> are worried about actually *has* the property you are filtering on.

**Off-target reads.** Every capture protocol leaves stray reads across the
genome. Without a target restriction the caller walks the whole genome and
calls on them. On that same 4.5 Mb panel:

| | raw candidates | PASS calls |
|---|---|---|
| without `--intervals` | 635,414 | 127,081 |
| with `--intervals` | 20,761 | ~1,400 |

Two orders of magnitude of noise, every one of which looks like a result.

**Tag, don't delete.** The foldback filter *tags*. Genuine foldback
inversions occur in cancer genomes, so the pattern is *evidence* of
artefact, not proof — and `INFO/FBMATCH` records the matched fraction for
every insertion tested.

> **Pearl.** Deleting data because it is *probably* wrong destroys the
> reviewer's ability to disagree with you. Annotate instead, and make the
> threshold visible.

### VI. The interface is part of the instrument

**Two patients in one directory.** Auto-discovery originally took the first
FASTQ pair as tumour and the second as its matched normal. Point it at a
directory holding two *different patients* and it produces somatic calls on
the genetic difference between two unrelated people — their differing
germline variants reported as somatic, their shared real mutations
subtracted. The output looks entirely ordinary, and PCGR will tier it.

**Validate everything before starting anything.** The hybrid DNA+RNA
submission validates *both* halves before queueing either. Before that fix,
a bad RNA half left a DNA run already executing behind a page reading
"Nothing was started".

**Never address things by position.** Run outputs were linked by index into
a list of files. But those reports are meant to be readable *while the run
is still going*, so the list grows under an already-rendered page — and a
link labelled "Read quality — SAMPLE_B" started serving SAMPLE_B's
*coverage* report. HTTP 200, wrong file, no error anywhere. Ids are now a
digest of each file's path.

> **Pearl.** If a collection can change between rendering a link and
> following it, the link must address **content, not position**. This class
> of bug is invisible to "does it return 200?" testing.

**Patient identifiers never touch the command line.** They live in one JSON
file beside the run, read only when a report is produced. The pipeline
subprocess is never told them, so they cannot reach a tool's stdout, a
manifest, or a crash trace.

**Serial queue, deliberately.** One run at a time. BWA-MEM2 wants ~32 GB
for hg38 and STAR ~30 GB; two concurrent runs do not go faster, they
thrash or OOM.

---

## 10. Architecture patterns worth stealing

**Stages hand off through a manifest.** Each stage writes JSON naming its
outputs; the next reads it. Not filename guessing — an explicit,
inspectable contract. It also means stage 3 can re-run without redoing
stage 1.

**Share the expensive thing.** The genome and its indices live in one
directory every run resolves against. It used to be a download and a 1–2
hour index build *per output directory*, for identical bytes.

**Resumption must refuse on drift.** `--resume` reuses finished steps — but
compares the previous run's parameters first and **aborts** if they differ.
Silently reusing a BQSR BAM built with different known-sites would produce
a call set assembled from two analyses, and would finish successfully.

**Separate environments for conflicting dependencies.** GATK pins a JDK
range; STAR and Arriba do not care; PCGR needs R. Three conda environments,
and the job runner activates the right one per job. One environment means a
solver resolving unrelated constraints, and what gives is the pin nobody is
watching.

**Make expensive optional things opt-in.** The RNA branch costs ~32 GB and
an hour of CPU, so it installs only with `--with-rna`.

**Idempotent installers.** Every step checks first and skips what is
complete; an interrupted install resumes by re-running the same command.

**Report-only update checks.** A tool that *tells* you a newer database
exists, and never installs it, because upgrading mid-cohort makes reports
disagree for reasons unrelated to the samples.

**Write down what you did not build.** [CNV_SCOPE.md](CNV_SCOPE.md) and
[RNA_SCOPE.md](RNA_SCOPE.md) describe capabilities *designed and
deliberately not implemented*, with what each would cost. This is a real
deliverable: it stops the half-built version appearing later, and tells the
next person which half was the hard one.

---

## 11. Rebuilding this from scratch

A build order that keeps you honest, each step testable before the next.

**Phase 0 — environment.** One conda env with `fastp fastqc bwa-mem2
samtools gatk4 bcftools snpeff`. Pin `openjdk>=21,<24`. Verify with
`--version` on each.

**Phase 1 — reference data.** hg38 FASTA + `.fai` + `.dict` + bwa-mem2
index. Then dbSNP, Mills/1000G indels, gnomAD, a panel of normals, a
common-SNP VCF. **Check every one's contig naming against the genome
before going further.**

**Phase 2 — QC → BAM.** fastp, then bwa-mem2 piped to samtools sort, then
index. Test on a few thousand reads. Confirm the read group is present
(`samtools view -H`).

**Phase 3 — the GATK chain.** MarkDuplicates → BaseRecalibrator →
ApplyBQSR → Mutect2 → FilterMutectCalls. Add the target BED as soon as you
have one; compare call counts with and without it, and be shocked.

**Phase 4 — make it honest.** Before adding features: the coverage report,
the run manifest, and the parameter-drift check on resume. This is the
phase people skip, and it is the one that makes the output trustworthy.

**Phase 5 — annotation and interpretation.** COSMIC (rename contigs
first), SnpEff, then PCGR in its own environment with its bundle and VEP
cache. Set the TMB denominator explicitly the first time you turn TMB on.

**Phase 6 — configuration.** Only once you have run two different assays
will the profile system make sense. Build it then, not before.

**Phase 7 — RNA, if you need fusions.** Second environment, GENCODE, STAR
index (record what it was built from), STAR with the chimeric block, Arriba
with its blacklist, then the library-adequacy report *before* you trust any
negative.

**Phase 8 — the interface.** Serial queue, path confinement, patient
details off the command line, and outputs linked by content rather than
position.

---

## 12. A checklist for your own pipeline

- [ ] For every step: **if this silently did nothing, how would anyone
      know?** Fix the ones with no answer.
- [ ] For every threshold and denominator: where did that number come from?
      Is it right for *this* assay, or for the author's?
- [ ] Does an empty result look different from an unmeasured one?
- [ ] Are cheap checks done before expensive work?
- [ ] Does every run record what produced it — settings, versions,
      reference files?
- [ ] Do your reports name what they could **not** measure?
- [ ] Does every reference file's contig naming match the genome?
- [ ] Are downloads verified by size, not existence?
- [ ] Does the interface make the dangerous action harder than the safe
      one?
- [ ] Are outputs addressed by content, not by position in a list?
- [ ] Have you written down what you deliberately did not build?
- [ ] Is there an orthogonal method to check against? If not, say so
      loudly in the output itself.

---

## 13. The honest ending

This pipeline is validated on **one sample, one assay, tumour-only**. The
RNA branch has no validation at all. Every report it produces says so.

That is not modesty, it is the same principle as everything else here:
**the system must not claim more than it has earned**, because the person
reading the output cannot tell the difference — and in this domain, they
may act on it.

See *Read before clinical use* in [README.md](README.md), plus
[CNV_SCOPE.md](CNV_SCOPE.md) and [RNA_SCOPE.md](RNA_SCOPE.md), for exactly
what is and is not established.
