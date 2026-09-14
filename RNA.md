<!-- Created by Brainstorm, 2026. -->
# RNA — running the fusion branch

Fusion detection from RNA, for cancer panels. `ALK`, `ROS1`, `RET`,
`NTRK1/2/3`, `FGFR`, `NRG1`, `MET` exon 14 skipping — the findings that
change treatment and that a DNA panel cannot see.

Read [RNA_SCOPE.md](RNA_SCOPE.md) before using any of it clinically. It says
what this branch deliberately does not do.

---

## Install it

The RNA branch is **opt-in**. A laboratory running only DNA panels should
not pay 32 GB and an hour of CPU for an index it will never open.

```bash
python3 install_pipeline.py --with-rna --rna-read-length 150
```

That adds, on top of a normal install:

| | Location | Size |
|---|---|---|
| RNA environment (STAR, Arriba) | `~/miniconda3/envs/cancer_rna` | 1.5 GB |
| GENCODE annotation | `~/data/references/gencode/` | 1.5 GB |
| STAR index | `~/data/references/star_hg38_<len>/` | ~30 GB |

It reuses the hg38 FASTA the DNA branch already downloaded — the shared
`--reference-dir` design pays for itself here.

If you omit `--read-length` and `--star-index` at run time, the pipeline
**discovers** the installed index rather than guessing a path: with exactly
one installed it uses it and says so; with several it names them and asks
you to be explicit.

**Match `--rna-read-length` to your instrument.** The read length is baked
into the index (STAR's `--sjdbOverhang`), an index built for 100 bp reads
works on 150 bp reads, and what it loses is junction sensitivity — silently.
`align_rna.py` records what each index was built from and warns when a run's
read length or annotation does not match.

```bash
python3 install_pipeline.py --check      # says which of the three parts exist
```

---

## Run it

```bash
conda activate cancer_rna

python fusion_calling.py \
    --panel illumina-tso500-rna \
    -i fastqs/ -o rna_results/ --auto-discover \
    --reference ~/data/references/hg38/Homo_sapiens_assembly38.fasta \
    --gtf ~/data/references/gencode/gencode.v44.primary_assembly.annotation.gtf \
    --star-index ~/data/references/star_hg38_150 \
    --read-length 150 --threads 16
```

Six steps: QC → STAR → sort → Arriba → library QC → fusion report.

In the web interface it is the **New RNA run** form, which offers only RNA
profiles and refuses a DNA one.

### Chaining from one QC run

A DNA and an RNA library from one specimen share stage 1:

```bash
python pipeline_orchestrator.py --assay rna \
    --panel illumina-tso500-rna \
    --stage1-args "-i raw/ -o qc/ --threads 16" \
    --stage3-args "-o rna_results/ -r <hg38.fa> --gtf <gencode.gtf> --threads 16"
```

`--panel-bed` is a DNA option and is **ignored with a warning** under
`--assay rna`: an RNA run is scoped by its annotation (`--gtf`), not by a
target BED, and `fusion_calling.py` has no such flag.

`--assay rna` forces stage 2 to be skipped. That is not an optimisation:
bwa-mem2 cannot align across an exon–exon junction, so it would discard
exactly the reads a fusion is evidenced by, and the resulting BAM looks
ordinary.

---

## The profiles

```bash
python fusion_calling.py --list-panels        # RNA ones are marked [RNA]
python fusion_calling.py --describe-panel archer-fusionplex
```

| Profile | Manufacturer | Chemistry |
|---|---|---|
| `generic-rna-capture` | any | probe capture |
| `generic-rna-amplicon` | any | anchored multiplex PCR |
| `rna-total` | any | whole transcriptome |
| `illumina-tso500-rna` | Illumina | capture |
| `archer-fusionplex` | Invitae / ArcherDX | anchored PCR |
| `thermo-oncomine-rna-fusion` | Thermo Fisher | anchored PCR |

They set the library-size and mapping-rate floors the QC verdict is measured
against, and those differ by an **order of magnitude** between chemistries:
2 M reads is a normal, adequate anchored-PCR library and a badly failed whole
transcriptome. Getting that wrong means either passing unusable libraries or
failing good ones.

Profiles live in the same registry as the DNA ones — `--panel-file`,
`~/.config/cancer_pipeline/panels/`, `--new-panel-template` all work exactly
as [PANELS.md](PANELS.md) describes. An RNA profile has
`"chemistry": "rna-capture" | "rna-amplicon" | "rna-total"`, and that string
is what routes a run to this branch. The engines each refuse the other's
profiles by name rather than running with settings that have no meaning.

---

## The two reports

Both are written per sample, as self-contained HTML plus JSON and TSV, and
both are readable **while the run is still going**.

### `rna_qc/<sample>.rna_qc.html` — can a negative be reported?

This is the one that matters, and it exists because of a single failure
mode:

> A degraded FFPE RNA extraction produces a clean, well-formed, **empty**
> fusion table, indistinguishable from a specimen that genuinely carries no
> fusion.

So it measures the library — input reads, unique mapping, reads lost as "too
short", mitochondrial and rRNA fractions, junction and chimeric counts — and
ends with a sentence about *reportability* rather than a table of numbers:

> *3 of 6 checks did not clear their bar. An empty fusion table from this
> library should NOT be reported as a negative without explaining these:
> Library size; Uniquely mapped; Reads lost as 'too short'.*

It also lists, by name, what it **cannot** measure — duplication rate (not a
quality metric on RNA: an abundant transcript produces identical reads
legitimately) and DV200/RIN (an instrument measurement that cannot be
recovered from FASTQ) — so their absence is not read as a pass.

**Record the DV200.** The web form has a field for it. It is the best single
predictor of whether fusion detection could have worked, and no amount of
computation recovers it.

### `fusion_report/<sample>.fusions.html` — what the calls mean

Arriba ranks by its own confidence, which is a statement about *evidence*. A
high-confidence read-through between two housekeeping genes therefore sorts
above a medium-confidence `EML4–ALK`. This report reorders by **clinical
salience**:

1. calls involving a gene with established fusion-directed therapy or
   diagnostic significance,
2. same-gene splice events — **MET exon 14 skipping** is one of the most
   actionable findings a lung panel produces, and in a table of gene pairs it
   reads as `MET--MET` and is trivially missed,
3. caller confidence,
4. read support.

Calls below the confidence bar are **marked, never dropped**: a laboratory
investigating a clinically expected fusion needs to see that it was called
weakly rather than be told nothing was found.

The actionable-gene list is a **highlighting aid, not a knowledge base**. It
does not claim a marked fusion is actionable in this tumour type — the same
fusion is predictive in one tissue and incidental in another — and it does
not claim an unmarked fusion is unimportant.

### `pcgr/` — the clinical interpretation

PCGR 2.x accepts RNA fusions as a molecular input in their own right, so an
RNA run gets a **tiered clinical interpretation**, not only this pipeline's
own ranked table. `fusion_calling.py` writes `<sample>.pcgr_fusions.tsv` in
PCGR's schema every run, and the web interface runs PCGR automatically when
a reference bundle is configured.

```bash
# Fusions alone. No VEP cache needed -- there is no VCF to annotate.
python pcgr_report.py --input-rna-fusion rna_out/fusion_report/S1.pcgr_fusions.tsv \
    --output-dir rna_out/pcgr --sample-id S1 \
    --pcgr-refdata-dir ~/data/pcgr/<release> --tumour-site 15

# Better: the specimen's DNA and RNA in ONE report.
python pcgr_report.py --input-vcf dna_out/annotated/S1.cosmic.vcf.gz \
    --input-rna-fusion rna_out/fusion_report/S1.pcgr_fusions.tsv \
    --vep-dir ~/data/vep_cache --lift-tags --tumour-only \
    --output-dir combined/ --sample-id S1 \
    --pcgr-refdata-dir ~/data/pcgr/<release> --tumour-site 15
```

On the **New RNA run** form, pick the matching DNA run under *Linked DNA
run* and the combined report is produced for you — which is how a specimen
is reported clinically, rather than as two documents somebody reconciles by
eye.

Two behaviours worth knowing:

- **PCGR rejects an empty fusion file.** "No fusions found" is a result, not
  an error, so the pipeline writes no PCGR input in that case and says so.
  The run page then explains why there is no PCGR report.
- Only calls at or above `--fusion-min-confidence` are handed to PCGR.
  Unlike the HTML report, which *marks* weak calls, anything passed here
  becomes an entry in a clinical interpretation with no way to show it was
  weak.

### The PDF

The web interface's **Generate PDF report** produces a patient-attributed
document for RNA runs too, from `webapp/rna_report.py`. Library adequacy
comes **before** the fusion table in it, deliberately: a reader who sees an
empty table first has already concluded "negative" by the time they reach
the caveat.

---

## Keeping it current

```bash
python3 check_db_updates.py --print
```

Reports two RNA-specific things, and only when the branch is installed:

- **GENCODE release.** A newer one changes gene names and transcript sets,
  so do not mix releases within a cohort.
- **STAR index consistency.** Whether each index was built from the
  annotation that is still installed. Nothing else asks this question
  outside of a run, and getting it wrong is silent: STAR runs, the mapping
  rate looks normal, and the junctions it knows about are the old ones.

The web interface's **Databases** page offers both as buttons, and the
GENCODE one runs *both* steps: upgrading the annotation without rebuilding
the index is not a shortcut, it is the inconsistency the check exists to
catch.

Upgrading the annotation means rebuilding the index, because the annotation
is baked into it:

```bash
python3 install_pipeline.py --only gencode star-index --force --rna-read-length 150
```

The installer pins a GENCODE release rather than tracking the newest, so
two machines installed a week apart agree about which transcripts exist.
The checker reports drift; it never upgrades by itself.

## Things that will bite you

- **The STAR index bakes in the annotation and the read length.** Change the
  GENCODE release and the index must be rebuilt. Runs check the index's build
  record and say so; an index built elsewhere carries no record, and that is
  reported too.
- **Contig naming.** GENCODE's primary-assembly GTF uses UCSC names (`chr1`,
  `chrM`) and matches the installed hg38. An Ensembl GTF (`1`, `MT`) does
  not, and builds an index over nothing. The installer checks and warns.
- **Arriba's blacklist.** It ships inside the conda package and is normally
  found automatically. Without it, recurrent read-through artefacts are
  reported as high-confidence fusions and the table fills with noise. The run
  log says which reference files were located; the PDF states their absence.
- **Arriba wants the unsorted BAM.** It reads mates together, and coordinate
  sorting separates them. `fusion_calling.py` hands it STAR's own output and
  produces the sorted copy separately, for QC and IGV. Pointing Arriba at the
  sorted BAM is a supported command that finds far fewer fusions.
- **STAR wants ~32 GB of RAM** to load the hg38 index. The web interface's
  one-run-at-a-time queue is what keeps that from colliding with anything
  else on the machine.
- **One QC run feeds both branches.** `align_rna.py` reads the
  `fastq_qc_clean.py` manifest through the DNA engine's own reader, so a
  sample whose QC was skipped because its output already existed is
  recovered rather than silently dropped.
- **Duplicates are never marked**, and there is no flag to enable it. Two
  reads at one position in a highly expressed gene are two observations of an
  abundant transcript, not one molecule counted twice.
