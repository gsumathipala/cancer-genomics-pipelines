# Installation Guide

This guide installs every external tool required by the three pipeline scripts
(`fastq_qc_clean.py`, `align_reads.py`, `comprehensive_variant_calling.py`),
the orchestrator (`pipeline_orchestrator.py`), and the optional clinical
reporting step (`pcgr_report.py`).

Conda/Bioconda is the recommended route because it installs the exact binaries
(plus Java, for GATK/SnpEff) into one isolated environment, avoiding system
package-manager conflicts.

## Tool inventory

| Tool       | Needed by                                | Purpose                                | Bioconda package |
|------------|------------------------------------------|----------------------------------------|------------------|
| fastp      | fastq_qc_clean.py, comprehensive_variant_calling.py | QC, adapter/poly-G trimming | `fastp` |
| FastQC     | fastq_qc_clean.py, align_reads.py | Quality reports (raw, cleaned, BAM)    | `fastqc` |
| bwa-mem2   | align_reads.py, comprehensive_variant_calling.py | Read alignment to reference | `bwa-mem2` |
| samtools   | align_reads.py, comprehensive_variant_calling.py, coverage_report.py | BAM sort/index/stats; per-base depth for the coverage check | `samtools` |
| GATK (4.x) | comprehensive_variant_calling.py | MarkDuplicatesSpark, BQSR, Mutect2, FilterMutectCalls | `gatk4` |
| bcftools   | comprehensive_variant_calling.py | COSMIC annotation, overlap counts, stats | `bcftools` |
| SnpEff     | comprehensive_variant_calling.py | Gene/consequence annotation      | `snpeff` |
| PCGR       | (optional) pcgr_report.py, comprehensive_variant_calling.py step 13 | Clinical report: actionability tiers; TMB / MSI / signatures only if the matching `--pcgr-estimate-*` flag is passed | `pcgr` |
| coverage_report.py | (optional) comprehensive_variant_calling.py step 12 | Which target regions the sequencing actually reached, so an absent variant can be told from an unsequenced one; needs only samtools | in this bundle |
| multiqc    | (optional, suggested by fastq_qc_clean.py) | Aggregate reports          | `multiqc` |
| Flask      | (optional) webapp/ only            | Web interface                       | `flask` (pip or conda) |
| wget/curl  | system (reference auto-download)   | Download hg38 reference from Broad, once, into the shared `--reference-dir` | system package |

PCGR is optional and deliberately left out of the main environment below: it
pulls a large dependency set and additionally needs two data downloads that
must match the installed PCGR version — its **reference data bundle**
(~7 GB, passed as `--pcgr-refdata-dir`) and an **Ensembl VEP cache**
(~24 GB, `--vep-dir`).

**PCGR is not on bioconda.** It is published on its own `pcgr` Anaconda
channel, and the supported install builds **two** environments from
version-pinned lock files — `pcgr` (Python + VEP workflow) shells out to
`pcgrr` (R + quarto reporting):

```bash
PCGR_VERSION=2.3.2
LOCK="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock"
conda create -y -n pcgr  --file <(curl -sSL "$LOCK/pcgr-linux-64.lock")
conda create -y -n pcgrr --file <(curl -sSL "$LOCK/pcgrr-linux-64.lock")
```

The two environments must be **siblings**: PCGR locates `pcgrr` as
`dirname($CONDA_PREFIX)/pcgrr`. Named envs under `envs/` satisfy this; use
`--pcgrr_conda <name>` if you rename one.

> **PCGR cannot be run from the pipeline environment.** With
> `cancer_pipeline` active, `pcgr` is not on `PATH`, so step 13 skips with a
> warning. Prepending PCGR's `bin/` to `PATH` is *not* a workaround — PCGR
> resolves its VEP plugin directory from `$CONDA_PREFIX/share`, so a real run
> then fails with `FileNotFoundError: No ensembl-vep directories found`.
> Run the pipeline first, then `conda activate pcgr` and invoke
> `pcgr_report.py` separately (it is standard-library-only, so PCGR's own
> Python runs it).

Version mismatch between the PCGR binary and its data bundle is the usual
cause of a failed first run: the software, the reference bundle and the VEP
cache are locked to one another, so upgrade all three together. Check the
release's installation page for the bundle and VEP version that pair with it.
When PCGR or its bundle is absent the pipeline simply skips the step with a
warning — it is never fatal.

The **analysis scripts** (`fastq_qc_clean.py`, `align_reads.py`,
`comprehensive_variant_calling.py`, `pipeline_orchestrator.py`,
`pcgr_report.py`) use only the Python standard library (3.8+) — no `pip`
packages at all. **Flask is required only by the optional web interface**, and
PDF generation adds nothing further because `webapp/pdfwriter.py` writes PDFs
using the standard library alone.

Two machine-readable manifests accompany this guide, and an installer should
consume them rather than parsing this prose so the two cannot drift apart:

| File | Covers |
|------|--------|
| `environment.yml` | The whole conda environment (all tools + Flask) |
| `requirements.txt` | Python packages only (Flask); analysis scripts need none |

## 0. Automated install (recommended)

`install_pipeline.py` performs everything this guide describes, and encodes
the traps below as checks rather than paragraphs a reader has to notice.

```bash
python3 install_pipeline.py --check     # what is missing; changes nothing
python3 install_pipeline.py --dry-run   # print every command, run none
python3 install_pipeline.py             # install what is missing
```

About 70 GB and roughly two hours, most of it waiting on downloads. Re-running
is safe: each step checks first and skips what is already complete, so an
interrupted install is resumed by running it again.

It verifies rather than assumes — the JDK is in the 21–23 window, indices sit
beside their data, every VCF has its `.tbi`, `pcgr` and `pcgrr` are siblings,
and the pipeline script starts. `--check` alone is a useful health check on a
machine you inherited.

COSMIC needs a registered account, so it cannot be fetched unattended.
Download it yourself and pass `--cosmic <file>`; the script renames its
contigs to UCSC style (without which COSMIC annotation silently matches
nothing) and indexes it.

The rest of this document explains what the installer does and why, and is
what you need when installing by hand or debugging a failure.

## 1. Conda prerequisites

If you do not already have conda, install Miniforge (includes `mamba`):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# then log out / back in, or:  source ~/.bashrc
```

Configure the channels once. Order matters (strict channel priority gives
consistent, tested build combinations):

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

> Note: `~/.condarc` applies to all environments. If you already have
> `defaults` channels listed, prefer keeping `conda-forge` + `bioconda`.

## 2. Recommended: single environment for the whole pipeline

One command installs everything (Java and all binaries included; can take a
few minutes):

```bash
conda create -n cancer_pipeline \
    -y \
    python=3.12 \
    fastp \
    fastqc \
    bwa-mem2 \
    samtools \
    gatk4 \
    bcftools \
    snpeff \
    multiqc
```

Activate for every run:

```bash
conda activate cancer_pipeline
```

Run `which bwa-mem2 samtools gatk bcftools snpEff fastp fastqc` to confirm the
activation picked up the environment binaries (the scripts check `PATH`).

## 3. Alternative: per-tool installs

```bash
# Core QC + alignment + utilities
conda install -n cancer_pipeline -y fastp fastqc bwa-mem2 samtools
# Variant calling + annotation
conda install -n cancer_pipeline -y gatk4 bcftools snpeff
# Optional
conda install -n cancer_pipeline -y multiqc
```

If you prefer to avoid snpEff via conda (some users build it from source):

```bash
# from https://pcingola.github.io/SnpEff/
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
export PATH="$PWD/snpEff:$PATH"     # add to ~/.bashrc
```

## 4. SnpEff (and GATK) data setup

SnpEff needs its genome database before `comprehensive_variant_calling.py`
step 11 will run (the script skips keyed on the `snpEff` binary, but annotation
fails without the DB):

```bash
# bioconda builds default to /usr/local/share/snpeff (or $CONDA_PREFIX/share/snpeff on conda-forge builds)
snpEff download hg38   # or:  java -jar snpEff.jar download hg38
```

If the database was installed to a non-default location, point SnpEff at it
with the `SNPEFF_DATABASE`/`SNPEFF_DATA` env var or edit `snpEff.config`
(`data.dir`).

GATK needs a Java runtime, and the JDK version matters more than it looks.
`environment.yml` pins `openjdk>=21,<24`, and BOTH bounds are load-bearing:

* **>= 21** is snpEff 5.4's own requirement; the solve fails below it.
* **< 24** is GATK's. From JDK 24 (JEP 486) `Subject.getSubject()` throws
  unconditionally, and Hadoop calls it while Spark starts up, so every GATK
  Spark tool dies immediately. In this pipeline that is step 3,
  `MarkDuplicatesSpark`, which exits 3 with
  `java.lang.UnsupportedOperationException: getSubject is not supported`
  before it reads a single alignment. The usual JDK 18-23 escape hatch
  (`-Djava.security.manager=allow`) does not help -- from 24 the launcher
  rejects that flag too.

Do **not** rely on `gatk4` to pull a usable JDK in on its own: it declares only
a floor, so a fresh solve happily installs the newest openjdk and leaves you
with an environment that passes `gatk --version` and then fails at step 3.
Check what you actually got:

```bash
java -version          # must report 21.x, 22.x or 23.x -- not 24+
gatk --version         # should print the Java/GATK banner
```

## 5. Reference genome

The pipeline can auto-download hg38 (`--reference hg38`) from the Broad
Institute's public Google Cloud Storage bucket
(`gs://gcp-public-data--broad-references`, served over https, no login);
this requires `wget` or `curl` on your machine. The older S3 path that used
to be documented here is dead.

```bash
# Debian/Ubuntu
sudo apt-get install -y wget curl
# RHEL/Fedora
sudo dnf install -y wget curl
```

Or prepare hg38 yourself and pass the path:

```bash
# Broad public GCS bucket. Served uncompressed -- there is no .fasta.gz here.
BROAD=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
wget -q -O hg38.fa $BROAD/Homo_sapiens_assembly38.fasta
# The .fai and .dict are published alongside it -- downloading them is much
# faster than regenerating, and they must sit beside the FASTA:
wget -q -O hg38.fa.fai $BROAD/Homo_sapiens_assembly38.fasta.fai
wget -q -O hg38.dict   $BROAD/Homo_sapiens_assembly38.dict
# then pass --reference /path/to/hg38.fa
```

> Keeping Broad's original filenames (`Homo_sapiens_assembly38.fasta` plus
> `.fasta.fai` and `.dict`) avoids a renaming trap: GATK finds the dictionary
> by stripping the FASTA extension and appending `.dict`, so if you rename the
> FASTA to `hg38.fa` you must also rename the dictionary to `hg38.dict`.

> Where it ends up: `--reference hg38` is resolved against `--reference-dir`
> (default `~/data/references`, or `$PIPELINE_REFERENCE_DIR`) — the layout
> this installer writes — and an installed genome found there is used as it
> stands, with no download and no re-indexing. Only a genome missing from
> there is fetched, into that same shared directory. It used to land in
> `{output_dir}/reference/`, inside *that run's* output folder, so every new
> `-o` bought another 3 GB download and another 1–2 hour bwa-mem2 index; that
> now happens only if the shared directory cannot be written, and the run says
> so when it falls back. Reference indexes (`.fai`, `.dict`, `.bwt.2bit.64`,
> ...) are generated beside the FASTA and must stay there.

## 6. Databases (VCF resources)

Beyond the reference FASTA and the SnpEff hg38 database, the variant-calling
stage needs several VCF resources. **All of them must be on the hg38/GRCh38
build** to match the reference, and tabix-indexed (`.tbi`) for GATK/bcftools.

| Database | CLI flag | Source (hg38) | Needed when |
|----------|----------|---------------|-------------|
| COSMIC coding mutations | `--cosmic` | https://cancer.sanger.ac.uk/cosmic/download (free registration) | Step 9 annotation + overlap counts; skipped if omitted |
| gnomAD AF-only (germline) | `--germline-resource` | `af-only-gnomad.hg38.vcf.gz` (somatic-hg38 bucket) | **Essential in tumour-only mode**; Mutect2 filters germline calls with it |
| Panel of Normals | `--panel-of-normals` | `1000g_pon.hg38.vcf.gz` (somatic-hg38 bucket) | Optional; Mutect2 recurrent-artefact subtraction |
| Common biallelic SNPs | `--contamination-resource` | `small_exac_common_3.hg38.vcf.gz` (somatic-hg38 bucket) | Step 7 contamination estimate. Without it FilterMutectCalls assumes **zero** contamination rather than measuring it, so real contamination goes undetected |
| dbSNP | `--dbsnp` | `Homo_sapiens_assembly38.dbsnp138.vcf.gz` (Broad bundle) | BQSR known-sites; BQSR is **silently skipped** if no `--dbsnp`/`--known-indels` |
| Mills + 1000G gold indels | `--known-indels` | `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (+ `Homo_sapiens_assembly38.known_indels.vcf.gz`) | BQSR known-sites (with or instead of dbSNP) |

### Target intervals (capture panels and exomes)

`--intervals` is not a database, so there is nothing to download: it is the
BED your panel or exome vendor ships. Supply it whenever the library is
captured rather than whole-genome. Without it Mutect2 walks the entire
genome and calls on the off-target reads every capture protocol leaves
behind — on a 4.5 Mb panel that was 127,081 PASS calls, mostly noise, some
of them supported by a single read.

Pair it with `--interval-padding 50` (reads run past the target edges) and
`--min-depth 20` (GATK applies no depth floor of its own). The BED must use
the same contig naming as the reference — `chr1`, not `1`, for hg38.

Download sources (public, no login):

```bash
# gnomAD germline resource + Panel of Normals (somatic-hg38 bucket)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz{,.tbi}
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz{,.tbi}

# Common biallelic SNPs for the contamination estimate (same bucket, ~1.3 MB).
# Note this one is NOT in gcp-public-data--broad-references -- that path 404s.
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz{,.tbi}

# dbSNP, Mills indels, known indels (Broad hg38 bundle)
# NOTE: the genomics-public-data/references/hg38/v0 path returns 403 to
# anonymous callers. Use the gcp-public-data--broad-references bucket:
BROAD=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
wget $BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.gz{,.tbi}
wget $BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz{,.tbi}
wget $BROAD/Homo_sapiens_assembly38.known_indels.vcf.gz{,.tbi}
```

COSMIC is a manual download: register at https://cancer.sanger.ac.uk/cosmic,
download the **COSMIC Coding Mutation VCF (hg38)** and index it yourself
(also do this for any other resource you download without a `.tbi`):

```bash
# make sure the downloaded file is bgzip-compressed, then index
tabix -p vcf COSMIC_vXX_hg38.vcf.gz   # use `bgzip` first if it is plain .vcf
```

Example invocation wiring the resources together:

```bash
python comprehensive_variant_calling.py \
    --reference /data/references/hg38/hg38.fa \
    --cosmic /data/resources/hg38/COSMIC_vXX_hg38.vcf.gz \
    --dbsnp /data/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-indels /data/resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --germline-resource /data/resources/hg38/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals /data/resources/hg38/1000g_pon.hg38.vcf.gz \
    --contamination-resource /data/resources/hg38/small_exac_common_3.hg38.vcf.gz \
    --intervals /data/panels/mypanel_targets.bed --interval-padding 50 \
    --min-depth 20 \
    -o results/ --threads 16
```

### Where files and installs must live

**Tools / binaries.** The scripts locate every tool with `shutil.which`, so
there is no required directory — the binary just has to be on `$PATH`. With the
conda environment, activating it puts all binaries in `$CONDA_PREFIX/bin` on
PATH automatically. For a manual install (snpEff zip, Homebrew, apt), the
tool's `bin/` must be added to PATH in `~/.bashrc`. Verify with `which gatk
bcftools snpEff bwa-mem2 samtools fastp fastqc`.

**Reference genome.** Auto-download lands in the shared `--reference-dir`
(see §5), not in the run's output directory. Recommended layout — which is
what `install_pipeline.py` creates, and what `--reference hg38` looks for —
indexed files must stay beside the FASTA:

```text
/data/references/hg38/hg38.fa                       # decompressed FASTA
/data/references/hg38/hg38.fa.fai                   # samtools
/data/references/hg38/hg38.fa.dict                  # GATK
/data/references/hg38/hg38.fa.bwt.2bit.64           # bwa-mem2 (plus .amb/.ann/.pac/.sa)
```

**VCF databases.** Each flag takes an explicit path (only `--cosmic` is
existence-checked at startup; point the others at real absolute paths or GATK
errors later). No auto-discovery, so any directory works — keep each VCF's
tabix index adjacent:

```text
/data/resources/hg38/
├── COSMIC_vXX_hg38.vcf.gz            (+ COSMIC_vXX_hg38.vcf.gz.tbi)
├── Homo_sapiens_assembly38.dbsnp138.vcf.gz    (+ .tbi)
├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  (+ .tbi)
├── af-only-gnomad.hg38.vcf.gz                 (+ .tbi)
├── 1000g_pon.hg38.vcf.gz                      (+ .tbi)
└── small_exac_common_3.hg38.vcf.gz            (+ .tbi)
```

**SnpEff database.** `snpEff hg38` resolves the genome named `hg38` under
SnpEff's data directory (`data.dir` in `snpEff.config`):
- conda `snpeff` → `$CONDA_PREFIX/share/snpeff/data/hg38/` (or `/usr/local/share/snpeff/data/hg38/`);
- manual zip → `<snpEff_dir>/data/hg38/`.

Point `data.dir=` at a shared location in `snpEff.config` (or `SNPEFF_DATA`)
if you want it outside the install.

**One reference, shared.** The pipeline resolves `--reference hg38` against
`--reference-dir` (default `~/data/references`, where this installer puts it)
before it considers downloading anything, so every run reuses the one indexed
copy. Only a genome missing from there is fetched, and it is fetched into that
directory — never into the run's output directory, which used to mean a 3 GB
download and a 1–2 hour bwa-mem2 index for each new output folder.

**Adjacency rule.** GATK/bcftools never search directories for indexes — they
are only found next to their files: VCF → `<name>.vcf.gz.tbi`, BAM →
`<name>.bam.bai`, FASTA → `.fai`/`.dict`/bwa-mem2 files, SnpEff genome →
`<data.dir>/hg38/`.

Minimal sets to keep the download volume sane:

- **Tumour-only run (just calling):** reference + gnomAD (+ SnpEff hg38 for
  step 11).
- **+ COSMIC annotation:** add the COSMIC VCF.
- **+ BQSR:** add dbSNP and/or Mills indels.
- **+ contamination estimate (step 7):** add `small_exac_common_3.hg38.vcf.gz`
  (1.3 MB — the cheapest useful addition on this list).
- **+ MSI (step 8):** add the MSIsensor2 `models_hg38` directory. The only
  route to an MSI answer on a panel; PCGR omits it.
- **+ coverage statement (step 12):** no download — pass your panel's BED as
  `--coverage-bed`. It is what separates "sequenced, wild type" from "never
  sequenced" in the report, and nothing else in the run can tell them apart.

The pipeline has **14 steps**. Nothing here is required to reach a call set;
each missing resource disables one step and says so.

## 7. Verification (post-install)

```bash
conda activate cancer_pipeline
python --version                          # 3.8+
fastp --version
fastqc --version 2>&1 | head -1
bwa-mem2 version
samtools --version | head -1
gatk --version
bcftools --version | head -1
snpEff -version

# THE ONE THAT ACTUALLY CATCHES A BAD ENVIRONMENT.
# gatk --version passes on a JDK that GATK's Spark tools cannot use, so it
# proves less than it looks. Check the Java version itself:
java -version                             # must be 21.x-23.x, NOT 24+

# Optional, and only if you installed it (see the tool inventory):
conda activate pcgr && pcgr --version   # needs the sibling pcgrr env too
```

If `java -version` reports 24 or newer, step 3 will fail with
`UnsupportedOperationException: getSubject is not supported` after alignment
has already run. Fix it before starting a real run:

```bash
conda install -n cancer_pipeline -c conda-forge "openjdk=21"
```

Sanity smoke test of the whole chain (dry run, checks command construction but
downloads nothing if `--reference` avoids `hg38`):

```bash
python fastq_qc_clean.py --help >/dev/null && echo "fastq_qc_clean OK"
python align_reads.py --help >/dev/null && echo "align_reads OK"
python comprehensive_variant_calling.py --help >/dev/null && echo "comprehensive OK"
python pipeline_orchestrator.py --help >/dev/null && echo "orchestrator OK"
python pcgr_report.py --help >/dev/null && echo "pcgr_report OK"
python coverage_report.py --help >/dev/null && echo "coverage_report OK"
```

A dry run exercises the real command construction without touching your data,
the network or the disk — it will not download the reference or write any
output:

```bash
python comprehensive_variant_calling.py \
    --input-dir raw_fastqs/ --output-dir results/ --reference /path/to/hg38.fa \
    --tumour-sample T1 --tumour-r1 T1_R1_001.fastq.gz --tumour-r2 T1_R2_001.fastq.gz \
    --dry-run
```

## 7a. Manual page (optional)

The four pipeline scripts, `pcgr_report.py` and `coverage_report.py` are
documented in one section-1 man page:

```bash
# Read it in place
man ./cancer-dna-pipeline.1

# Install for the current user
mkdir -p ~/.local/share/man/man1
cp cancer-dna-pipeline.1 ~/.local/share/man/man1/
man cancer-dna-pipeline

# Or system-wide
sudo install -m644 cancer-dna-pipeline.1 /usr/share/man/man1/
```

It covers every command-line option, the three input modes, the manifest
handoff between stages, what each `--skip-steps` name implies, which failures
are fatal, and the known caveats.

## 7b. Web interface (optional)

`webapp/` is a Flask front end: fill in patient details and run parameters,
watch the pipeline progress, and produce a PDF report from the result.

**It manages the conda environments itself.** Start it from any shell — base
conda, a plain terminal, a systemd unit. The pipeline is launched in
`cancer_pipeline` and PCGR afterwards in `pcgr`, each with its own `bin/` on
PATH and `CONDA_PREFIX` set, because neither can run in the other's
environment. Environments are found by name beside the running one, falling
back to `$CONDA_ROOT` and the usual `~/miniconda3`, `~/anaconda3`,
`~/miniforge3` installs; `PCGR_PYTHON` overrides the PCGR lookup for an
unusual layout.

**A blank form still runs the full flow.** Any resource field left empty is
filled in at submit time from what is installed — dbSNP, known indels, gnomAD,
the panel of normals, the contamination sites, COSMIC, the PCGR bundle and the
VEP cache — with TMB and signature estimation on. Whatever was filled in is
listed on the job page. Tick **minimal run** to opt out; that is the supported
way to run without a resource, rather than clearing a field and hoping someone
notices it was dropped.

The two things it cannot supply for you are the target BED (`--intervals`,
which is your assay vendor's file) and a matched normal.

**Copy number is not implemented.** The pipeline calls SNVs and short indels
only. CNV is designed in `CNV_SCOPE.md` and waits on normals sequenced on your
own assay — every tool it needs is already installed or one conda package away.
Note that the panel of normals used by Mutect2 is a different artefact and does
not serve this purpose.

**Naming the tumour site changes the tiers.** PCGR ranks actionability against
the tissue: the same variant can be predictive in one site and only oncogenic
in another, so the site decides which evidence counts as on-label. The form
offers all 31 of PCGR's sites, and **"Any / not specified" (code 0) is the
default** — a legitimate choice that reports every tier without asserting a
tissue nobody confirmed. On the command line it is `--pcgr-tumour-site <code>`.
The chosen site is recorded in the run record, by name rather than by code.

**MSI needs `--msi-models`, not PCGR.** MSIsensor2 scores microsatellite
instability from the tumour BAM alone, which is the only route to an MSI answer
on a targeted panel. The installer fetches the hg38 models
(`~/data/msisensor2/models_hg38`, 251 MB) and the web form fills the path in.
Clear it to omit MSI. Upstream calls MSI-H at **≥20%** for tumour-only data —
the 3.5% figure belongs to the older paired caller and would call almost
anything unstable. The site count matters as much as the score: a panel carries
only the microsatellites inside its targets, and the pipeline flags a run with
fewer than 200 covered loci.

**On FFPE material, keep the read-orientation model.** Formalin deaminates
cytosine, producing C>T and G>A changes that are damage rather than biology,
and step 6 exists to catch them — on the FFPE panel this was built against it
removed 3,222 calls, 1,824 of them C>T/G>A. Fixation also fragments the DNA,
which is why FFPE runs show short inserts and foldback artefacts. Residual
damage still reaches PASS: 59% of C>T/G>A PASS calls sat below 10% VAF against
37% of everything else, so `--min-allele-fraction` earns its place on FFPE.

**TMB, MSI and signature fitting are each opt-in and each can be switched off.**
PCGR computes none of them unless asked and says nothing when it has not been,
so the form offers all three as checkboxes: TMB and signatures start ticked,
MSI starts unticked. Unticking any of them genuinely omits it.

Two caveats the job page will also tell you:

* **MSI does nothing on a panel.** PCGR restricts it to WGS/WES
  tumour-control runs and logs `MSI status prediction can be applied for
  WGS/WES tumor-control assays only` before skipping. A targeted panel needs a
  dedicated caller such as MSIsensor2, which this pipeline does not yet carry.
* **Tumour-only TMB is unreliable** — PCGR warns about it directly. Residual
  germline variants inflate it, so treat the figure as internal QC rather than
  a reportable result until a matched normal is available.

**Several samples in one directory: choose what they are.** `--auto-discover`
on its own treats the first pair as the tumour and the **second as that
tumour's matched normal**, ignoring the rest — silently. Two patients in one
directory therefore produce somatic calls on the difference between two
people, and the report looks entirely ordinary. The webapp refuses that
outright; it asks you to pick:

- **Batch** — run every sample found as its own patient. One queued run per
  sample, in turn, each writing to its own subdirectory of the output
  directory and producing its own PCGR report.
- **The two samples are a tumour/normal pair from one person** — only valid
  for exactly two pairs.

The command-line script has no such guard: `--auto-discover` there still means
first-is-tumour, second-is-normal. Give it one directory per patient, or name
the samples explicitly.

**Freeing disk after a run.** A finished run leaves roughly 5.5 GB of
intermediates (cleaned FASTQs, and the aligned, duplicate-marked and
recalibrated BAMs) against about 70 MB worth keeping. The job page lists them
with real sizes and offers to delete them. Nothing is removed automatically,
and never while a run is queued or running. "Keep the recalibrated BAM" is on
by default — re-calling variants from it takes minutes, whereas rebuilding it
takes about an hour. The VCFs, the report, the metrics, the logs and the
manifest are never touched.

```bash
pip install flask          # the only dependency the web interface adds

# --allow-root: where your FASTQs/references live (repeatable)
# --runs-dir:   where runs, patient records and PDFs are written
# Both must already exist, or be creatable by YOU -- a path like /data
# usually needs root and will be refused with an explanation.
python webapp/app.py --allow-root ~/data --runs-dir ~/cancer_runs

# Simplest start: runs go to webapp_runs/ inside the repo (gitignored)
python webapp/app.py --allow-root ~/data

# then open http://127.0.0.1:5000
```

PDF generation needs nothing extra — `webapp/pdfwriter.py` writes PDFs using
only the standard library, so no reportlab/WeasyPrint/wkhtmltopdf install is
required.

**Before pointing it at real patient data, understand its limits:**

- **It binds to `127.0.0.1` by default and has no authentication.** Passing
  `--host 0.0.0.0` publishes both the patient details and the ability to start
  pipeline runs to everyone who can reach the port. The app warns if you do.
- **Never pass `--debug` with real data.** Flask's debugger allows arbitrary
  code execution from the browser.
- **Always pass `--allow-root`.** Without it the form accepts any path this
  user can read; with it, submitted paths are resolved (symlinks included) and
  confined to the directories you name.
- **Patient details are written only to the run directory**
  (`<runs-dir>/<job-id>/patient.json`) and are never passed to the pipeline, so
  they cannot appear in tool output, logs or manifests. `.gitignore` excludes
  the run directories — do not commit them.
- **One run executes at a time.** The pipeline is already parallel internally;
  further submissions queue rather than contend for the same CPUs and RAM.
- It is a single-analyst tool, not a multi-user or validated clinical system.
- **The clinical report runs after the pipeline, as a separate process.** If
  it fails, the run still reports `finished` — the call set is already on disk
  and represents the hours — and the job page carries the reason under
  `pcgr_status`. A run that ends without a report is not a run that failed.

FASTQs are referenced by path on the server rather than uploaded, since they
are routinely tens of gigabytes.

## 7c. Keeping an eye on database releases (optional)

`check_db_updates.py` asks whether anything installed has a newer release and
writes the answer to `~/.cache/cancer_pipeline/db_updates.json`. The web
interface reads that file and shows a notice on its front page.

```bash
python3 check_db_updates.py --print     # check now, show the result
```

Weekly, via cron:

```cron
17 6 * * 1  /usr/bin/python3 /path/to/check_db_updates.py
```

**It reports; it does not update.** It never writes to the data directory,
never downloads a database, and never touches a conda environment — every
check is an HTTP HEAD or a read-only API call. That is deliberate: a new VEP
cache changes the transcript set and a new COSMIC changes identifiers, so an
upgrade mid-project can make this week's report disagree with last week's for
reasons that have nothing to do with the sample. Upgrade on purpose, then
re-run the affected samples rather than mixing releases in one cohort.

It exits 0 even when checks fail, so a cron entry stays quiet; `--strict`
exits 1 when something needs attention, for a monitoring system that wants
that. A source it cannot reach is reported as **unknown**, never as up to
date — silence and good news must not look alike.

COSMIC is reported as **manual**: downloads need a registered account, so no
unattended check is honest. The installed version is shown beside the
release-notes link.

## 8. Troubleshooting / notes

- **Channel/solver trouble**: prefer `mamba` over `conda` for the big solve:
  `mamba create -n cancer_pipeline ...` (Miniforge ships mamba built-in).
- **GATK license**: GATK is released under the GATK open-source license
  (BSD-style for academic/commercial use); the bioconda `gatk4` package is the
  official release.
- **fastp version**: the scripts use `--average_qual`, `--trim_poly_g`,
  `--overrepresentation_analysis` — fastp >= 0.23 is required.
- **FastQC on BAMs**: `align_reads.py` (when FastQC is on PATH) runs FastQC on
  the sorted BAM; align_reads requires Java too. If you did not install FastQC,
  pass `--skip-fastqc` to both stages.
- **PCGR bundle mismatch**: the commonest PCGR failure is a reference data
  bundle whose version does not match the installed `pcgr`. Check the PCGR
  release notes for the bundle that pairs with your version. PCGR is optional
  throughout: if it or its bundle is missing, step 13 is skipped with a
  warning and the rest of the run is unaffected.
- **A clean report is not a negative result.** Nothing in a VCF says "this
  region was never sequenced": a capture dropout and a wild-type exon both
  produce no variant, and PCGR reports both as silence. Pass your panel BED as
  `--coverage-bed` and step 12 writes a per-sample HTML report naming every
  region that missed the depth threshold, with the exact stretches a variant
  could have hidden in. Without it the run says so and continues; the
  difference is that you no longer have any way to tell the two silences
  apart.
- **PCGR reports a suspiciously high TMB**: PCGR filters on depth/allele
  fraction taken from *INFO* tags, but Mutect2 writes those as per-sample
  *FORMAT* fields, so by default there is nothing to filter on. Pass
  `--pcgr-lift-tags` to copy them across (a separate VCF is written; the
  called VCF is untouched). Without it, treat TMB and MSI as indicative only.
- **Resuming an interrupted run**: `comprehensive_variant_calling.py --resume`
  reuses steps whose output is already complete. It is opt-in, and it refuses
  to run if the current arguments differ from the previous run's, so it cannot
  quietly build a call set out of two different analyses. Re-run with exactly
  the arguments of the interrupted run.
- **Not on PATH after activation?** Confirm the environment is active
  (`conda activate cancer_pipeline`) — the scripts resolve tools via `which`.

## 9. Installer reference

**This section is now implemented by `install_pipeline.py` (see §0).** It is
kept because it is the specification that script follows, and because it is
what you need when installing by hand, adapting the install to a different
layout, or working out why a step failed.

Everything an installer needs, in one place. The ordering matters: each step
depends on the ones above it, and every step is independently verifiable so a
partial install can be detected rather than discovered later at runtime.

### 9.1 Component matrix

| # | Component | Required? | Approx. size | Verify with |
|---|-----------|-----------|--------------|-------------|
| 1 | Conda/mamba | yes (or install tools by hand) | ~100 MB | `conda --version` |
| 2 | `environment.yml` env | yes | ~3–5 GB | `conda activate cancer_pipeline` |
| 3 | Pipeline scripts | yes | <1 MB | `python comprehensive_variant_calling.py --help` |
| 4 | Reference genome (hg38) + indices | yes | 3.1 GB + 17 GB indices | `ls hg38.fa.fai hg38.dict hg38.fa.bwt.2bit.64` |
| 5 | SnpEff database | for step 11 | ~450 MB | `snpEff databases \| grep -w hg38` |
| 6 | dbSNP / known indels | for BQSR (step 4) | ~1.6 GB | file exists + `.tbi` |
| 7 | gnomAD germline resource | for tumour-only calling | ~3 GB | file exists + `.tbi` |
| 8 | COSMIC VCF | for step 10 | ~1 GB | file exists + index |
| 8a | `small_exac_common_3.hg38.vcf.gz` | for step 7 (contamination) | ~1.3 MB | file exists + `.tbi` |
| 8b | MSIsensor2 models `models_hg38` | for step 8 (MSI) | 251 MB | directory of ~2,800 files |
| 9 | PCGR (2 envs) + bundle + VEP cache | for step 13 | ~6 GB envs + 31 GB data | `conda activate pcgr && pcgr --version` |
| 10 | Flask | for `webapp/` only | ~10 MB | `python -c "import flask"` |
| 11 | Man page | cosmetic | <1 MB | `man cancer-dna-pipeline` |

Only 1–4 are needed to call variants. Everything else degrades gracefully:
missing SnpEff, COSMIC, PCGR or bcftools cause the corresponding step to be
**skipped with a warning**, never to fail the run. An installer can therefore
offer a "minimal" and a "full" profile without risking a broken pipeline.

### 9.2 Ordered steps

```bash
# 1-2. Environment (mamba solves this much faster than conda)
mamba env create -f environment.yml
conda activate cancer_pipeline

# 3. Scripts: no build step; they are run in place. If installing to a
#    prefix, keep all five .py files in ONE directory -- the orchestrator
#    locates its siblings relative to itself, and the variant caller
#    imports pcgr_report from alongside it.

# 4. Reference (the pipeline can fetch hg38 itself, but indexing is what
#    the first real run otherwise spends its time on).
#
#    MEASURED on a 32-core / 93 GB host: bwa-mem2 index took 12 minutes,
#    single-threaded (extra cores do not help), with a PEAK RESIDENT SET OF
#    75 GB. That memory figure is the real constraint -- budget ~80 GB of RAM
#    or the index build is OOM-killed. The five index files are portable:
#    copy them between hosts rather than rebuilding on a smaller machine.
bwa-mem2 index /data/ref/hg38.fa
samtools faidx /data/ref/hg38.fa
gatk CreateSequenceDictionary -R /data/ref/hg38.fa

# 5. SnpEff database
snpEff download hg38

# 10. Web interface (only if wanted)
pip install -r requirements.txt
```

### 9.3 What an installer must not assume

- **Do not require PCGR.** Its bundle must match its version, it is tens of
  GB, and the pipeline is fully functional without it.
- **Do not bake in absolute paths.** Every resource is passed as a CLI flag;
  nothing is read from a fixed location.
- **Do not pre-create output directories.** The scripts create their own
  layout, and `--resume` decides what to reuse by inspecting it.
- **Do not run the web interface as a service on a shared host** without
  putting authentication in front of it. It binds `127.0.0.1` by design and
  has no login; see §7b.
- **Do not skip the index build to save install time.** Without the BWA-MEM2,
  `.fai` and `.dict` indices the first run builds them anyway, inside what the
  user thinks is their analysis.

### 9.4 Post-install self-test

`python3 install_pipeline.py --check` performs all of this and more — it also
checks the JDK version, that indices sit beside their data and that every VCF
has its `.tbi`. The manual equivalent follows.

This exercises the full command construction without touching real data,
downloading anything, or writing output — a good final installer step:

```bash
for s in fastq_qc_clean align_reads comprehensive_variant_calling \
         pipeline_orchestrator pcgr_report; do
  python "$s.py" --help >/dev/null && echo "$s OK" || echo "$s FAILED"
done

mkdir -p /tmp/selftest/fq && cd /tmp/selftest
touch fq/S_R1_001.fastq.gz fq/S_R2_001.fastq.gz ref.fa
python /path/to/comprehensive_variant_calling.py \
    --input-dir fq --output-dir out --reference ref.fa \
    --tumour-sample S --tumour-r1 fq/S_R1_001.fastq.gz \
    --tumour-r2 fq/S_R2_001.fastq.gz --dry-run
```

A successful dry run prints all **14** step banners and creates no BAMs or
VCFs. Add `java -version` to the loop above: it is the check that catches an
environment `gatk --version` will happily pass (see §7).
