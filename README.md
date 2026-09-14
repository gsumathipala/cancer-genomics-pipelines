<!-- Created by Brainstorm, 2026. -->
# Cancer Genomics Pipelines

> ## ⚠️ RESEARCH USE ONLY — NOT A VALIDATED DIAGNOSTIC SYSTEM
>
> This software must **not** be used as the sole basis for a clinical
> decision. It has been validated on **one sample, one assay, tumour-only**.
> No call it produces has been checked against an independent method. The
> RNA branch has **no validation at all**.
>
> Clinical deployment would need a validation set, orthogonal confirmation
> and a matched normal. Read
> **[Read before clinical use](#read-before-clinical-use)** before you do
> anything else — it lists, specifically, what is not established.

Somatic variant and fusion calling for Illumina cancer panels. Two branches
off one shared QC stage: **DNA** for SNVs and indels (GATK/Mutect2), and
**RNA** for gene fusions (STAR + Arriba), with clinical interpretation
through PCGR.

Everything needed to stand it up on a new machine: copy this directory
across, run one command, and you get the analysis environment, the reference
genome, six resource databases, the clinical reporter and its data. The RNA
branch is opt-in (`--with-rna` — see [RNA.md](RNA.md)).

> **New here?** [PIPELINE_ANATOMY.md](PIPELINE_ANATOMY.md) is a full
> dissection — what every tool does, how it works, the exact command that
> runs it, and the failure modes that are silent. It is written so someone
> who has never built one of these could rebuild it.

**Nothing here is patient data.** This bundle is code, configuration and
documentation only — no sequencing data, no reference databases, and no
clinical records have ever been committed to it. The reference data it needs
is downloaded at install time from its original sources, each under its own
licence.

---

## Install

```bash
# 1. Copy this directory to the target machine, then:
cd DNA_pipeline_installer

# 2. See what is missing. Changes nothing.
python3 install_pipeline.py --check

# 3. Install it.
python3 install_pipeline.py
```

About **70 GB** and, on a fast link, **roughly two hours** — mostly waiting on
downloads. The bwa-mem2 index is the one long CPU step, about 15 minutes.

Re-running is safe. Every step checks first and skips what is already complete,
so an interrupted install is resumed by running the same command again.

### Updating a machine that is already installed

```bash
# Copy the newer bundle across, then run it from there:
python3 /path/to/new_bundle/install_pipeline.py --update
```

That replaces the **scripts** of the installation this machine already has,
refreshes the man page and verifies. It downloads nothing and touches no conda
environment, so it takes seconds — the 70 GB of databases is not what changed.

The installer records where the code lives, so `--update` finds it on its own.
The first time you update an installation that predates that record, point at
it once with `--code-dir ~/DNA_pipeline_installer`; it is remembered after
that. `python3 install_pipeline.py --check` says whether the installed code
matches the bundle you are holding, and changes nothing.

Whatever gets replaced is copied into `.bundle-backup-<timestamp>/` inside the
installation first, and nothing is ever deleted. A file that changed since it
was installed — edited in place, or an older copy dropped over the top, which
cannot be told apart — is named individually as it is replaced; `--keep-local`
leaves those alone instead. Restart the web interface afterwards to pick the
new code up.

### Before you commit two hours

```bash
python3 install_pipeline.py --dry-run    # print every command, run none
```

### Requirements on the target machine

| | |
|---|---|
| conda or mamba | [Miniforge](https://conda-forge.org/download/) is the easiest source |
| Python 3.8+ | to run the installer itself; it uses only the standard library |
| wget or curl | for the downloads |
| ~70 GB free | and ~32 GB RAM to align against hg38 |
| Linux x86-64 | the PCGR lock files are platform-specific |

---

## COSMIC

COSMIC needs a registered account, so it cannot be downloaded unattended.
Register at <https://cancer.sanger.ac.uk/cosmic>, download the **GRCh38** VCF,
then:

```bash
python3 install_pipeline.py --cosmic /path/to/Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz
```

The installer renames its contigs from Ensembl style (`1`, `MT`) to UCSC
(`chr1`, `chrM`) and indexes it. That step matters more than it looks: without
it, COSMIC annotation matches **nothing at all** and the run still reports
success.

Everything else installs without an account.

---

## After installing

```bash
conda activate cancer_pipeline

# Command line
python comprehensive_variant_calling.py --help

# Web interface — no conda activation needed, it switches environments itself
python webapp/app.py --allow-root ~/data --allow-root /path/to/fastqs
# then open http://127.0.0.1:5000
```

**New run** asks one question first: *what are you looking for?*

| Pathway | Finds | Needs |
|---|---|---|
| **DNA** | SNVs and indels, tiered | the kit's target BED |
| **RNA** | gene fusions, splice events | GENCODE GTF + STAR index |
| **Both** | everything, in **one** report | both of the above |

The hybrid pathway is for a kit whose specimen yields two libraries
(TSO500, Oncomine Comprehensive, Archer). It queues both runs, links them,
and PCGR produces a single combined report — rather than leaving someone to
remember to link two runs after the fact. Profiles name their partner, so
choosing one half offers the other.

Every output of a run is then a labelled link on its page, grouped by the
question it answers, with *"can this result be trusted?"* deliberately
ahead of the findings.

`python3 install_pipeline.py --check` doubles as a health check afterwards, and
is the fastest way to tell whether an environment has drifted.

### Running more than one panel

```bash
python comprehensive_variant_calling.py --list-panels
python comprehensive_variant_calling.py --describe-panel thermo-oncomine-cav3

python comprehensive_variant_calling.py \
    --panel illumina-tso500 --panel-bed /data/panels/TSO500.bed \
    -i fastqs/ -o results/ -r hg38 --auto-discover
```

Naming the assay configures the run for its chemistry in one step:
interval padding, the allele-fraction and depth floors, whether duplicate
marking and BQSR are meaningful for that library, the UMI layout, and which
report statistics the target is large enough to support. Each of those is
**silent when it is wrong** — duplicate marking left on for an amplicon
panel finishes normally and calls from a fraction of the real depth.

Profiles ship for Illumina, Thermo Fisher, Agilent, Twist, IDT, QIAGEN,
Archer and Roche kits, plus generic shapes (`generic-capture`,
`generic-amplicon`, `generic-hotspot`, `ctdna-capture`, `wes`, `wgs`) for
anything not listed. **No profile ships a BED** — vendor BEDs are licensed,
version-specific content that comes from your kit. `--panel-bed` fills both
`--intervals` and `--coverage-bed`.

Explicit flags always beat the profile, and what the profile set is printed
at the start of the run and recorded in the manifest. Add your own kit with
`--new-panel-template`, or save a tuned run with `--save-panel-as`; drop
the JSON in `~/.config/cancer_pipeline/panels/` and it is found
automatically. See **[PANELS.md](PANELS.md)**.

Two things happen automatically once a target BED is known: the footprint
is measured and used as the TMB denominator (PCGR otherwise assumes 34 Mb,
so a 2 Mb panel reports a TMB 17x too low, silently), and the BED's contig
naming is checked against the reference — an Ensembl-style BED (`1`, `MT`)
is called over a renamed copy rather than making the coverage report state
that none of the panel was covered.

### RNA fusions

```bash
# Opt in at install time -- 32 GB and an hour of CPU more:
python3 install_pipeline.py --with-rna --rna-read-length 150

conda activate cancer_rna
python fusion_calling.py --panel illumina-tso500-rna \
    -i fastqs/ -o rna_results/ --auto-discover \
    --reference ~/data/references/hg38/Homo_sapiens_assembly38.fasta \
    --gtf ~/data/references/gencode/gencode.v44.primary_assembly.annotation.gtf \
    --star-index ~/data/references/star_hg38_150 --threads 16
```

A separate branch, not a mode: bwa-mem2 cannot align a read that crosses an
exon–exon junction, and every step of the DNA engine after alignment is
either meaningless or actively wrong on RNA. So RNA gets STAR, Arriba, its
own engine and its own two reports, sharing stage 1, the genome and the panel
registry. `pipeline_orchestrator.py --assay rna` chains it.

The failure mode it is built around: **a degraded FFPE RNA library produces a
clean, well-formed, empty fusion table that is indistinguishable from a true
negative.** So every run measures the library and states, in a sentence,
whether an empty result may be reported as a negative — and that statement
precedes the fusion table in every report.

PCGR interprets the fusions too — it takes them as a molecular input in
their own right — so an RNA run produces a tiered clinical report, and a run
linked to its DNA partner produces **one** report covering the specimen's
variants and fusions together.

`MET` exon 14 skipping is called out separately, because it is a splice event
within one gene rather than a fusion between two, and in a table sorted by
gene pair it reads as `MET--MET` and is missed.

SNVs and indels are **not** called from RNA, expression is **not** quantified,
and UMIs are extracted but not collapsed into consensus reads. See
[RNA_SCOPE.md](RNA_SCOPE.md) for what each omission would cost to add.

### Staying current

```bash
python3 check_db_updates.py --print       # any newer database releases?
```

Writes its answer to `~/.cache/cancer_pipeline/db_updates.json`. The script
**reports only** — nothing is downloaded or replaced, and no run reads it.

With the RNA branch installed it also checks the **GENCODE release** and
whether each **STAR index is still consistent with it**. That second check
has no other home: an index outliving the annotation it was built from is
completely silent — STAR runs, the mapping rate looks normal, and the
junctions it knows about are the old ones. Both appear on the **Databases**
page with an Update button, and upgrading GENCODE rebuilds the index in the
same action.

The web interface runs it **at startup** and whenever you press *Check now*, so
its front page is about today rather than about whenever someone last
remembered. Its **Databases** page lists every checked source and offers a
per-database **Update**, which runs `install_pipeline.py` for that one step and
re-checks afterwards. An upgrade is refused while a run is queued or in
progress, and asks you to type UPDATE first — because it is a real decision:
a new VEP cache changes the transcript set and a new COSMIC changes
identifiers, so mixing releases within a cohort makes reports disagree for
reasons unrelated to the samples.

```bash
# The same upgrades from a shell, which is what the buttons run:
python3 install_pipeline.py --only pcgr-data --force --vep-release 116
python3 install_pipeline.py --only resources --force
```

---

## What gets installed, and where

| | Location | Size |
|---|---|---|
| Analysis environment | `~/miniconda3/envs/cancer_pipeline` | 2.8 GB |
| Clinical reporter | `~/miniconda3/envs/pcgr` + `pcgrr` | 6.2 GB |
| hg38 + indices | `~/data/references/hg38/` | 20 GB |
| Resource VCFs | `~/data/resources/hg38/` | 6.4 GB |
| PCGR reference bundle | `~/data/pcgr/<release>/` | 7.3 GB |
| Ensembl VEP cache | `~/data/vep_cache/` | 24 GB |
| SnpEff database | inside the conda environment | 448 MB |
| MSIsensor2 models | `~/data/msisensor2/models_hg38/` | 251 MB |
| RNA environment (`--with-rna`) | `~/miniconda3/envs/cancer_rna` | 1.5 GB |
| GENCODE annotation (`--with-rna`) | `~/data/references/gencode/` | 1.5 GB |
| STAR index (`--with-rna`) | `~/data/references/star_hg38_<len>/` | 30 GB |
| The pipeline scripts | wherever you keep this bundle | 1 MB |
| Install record | `~/data/pipeline_install.json` | 20 KB |

Change the data location with `--data-dir`. The conda environments go wherever
conda keeps its environments.

The install record is how `--update` and `--check` find the scripts later: it
holds the directory they were installed to and a hash per file. Delete it and
they fall back to `--code-dir`.

Every run reuses that one indexed hg38: `--reference hg38` is resolved against
`--reference-dir` (default `~/data/references`) before anything is downloaded,
so a new output directory costs nothing. Only a genome missing from there is
fetched, and it is fetched into that shared directory rather than into the
run.

---

## Contents of this bundle

```
install_pipeline.py               the installer, and `--update` for a
                                  machine that already has one
environment.yml                   the analysis environment definition

comprehensive_variant_calling.py  the engine — all 14 steps, self-contained
align_reads.py                    stage 2 standalone
fastq_qc_clean.py                 stage 1 standalone
pipeline_orchestrator.py          chains the stages; --assay picks DNA or RNA
pcgr_report.py                    clinical report wrapper
coverage_report.py                target coverage: which regions a negative
                                  result is actually entitled to speak for
panel_profiles.py                 panel profiles: one name configures a run
                                  for its kit's chemistry
cancer-genomics-pipelines.1             man page

webapp/                           Flask front end (patient details, batch, cleanup)

install.md                        manual install, and the reference the installer follows
DATABASE_SETUP.md                 what each database is for, and where it goes
PANELS.md                         running more than one panel, and adding
                                  a kit of your own
PIPELINE_ANATOMY.md               how the whole thing fits together, and
                                  what it teaches about writing one
CNV_SCOPE.md                      copy-number: designed, not built — read before starting it

check_db_updates.py               tells you when a database has a newer release

align_rna.py                      RNA stage 2: STAR, chimeric detection on
fusion_calling.py                 the RNA engine -- all 6 steps
rna_qc_report.py                  RNA library QC: whether a negative is
                                  interpretable at all
fusion_report.py                  fusion calls ranked by clinical salience
environment-rna.yml               the RNA environment definition
RNA.md                            running the RNA branch
RNA_SCOPE.md                      what the RNA branch does NOT do, and why
```

**Keep the Python files in one directory.** The orchestrator locates its
siblings relative to itself, and the variant caller imports `pcgr_report` and
`coverage_report` from alongside it. The flat layout is deliberate.

---

## What the installer checks that a manual install will not

Each of these was found the hard way, on a real install, and each is silent or
misleading when it goes wrong:

- **The JDK version.** `gatk4` declares only a floor, so a fresh solve installs
  the newest JDK. GATK's Spark tools cannot run on 24+, and the failure lands at
  duplicate marking — *after* alignment has already run. `gatk --version`
  passes either way. The installer parses `java -version` and enforces 21–23.
- **PCGR's two environments must be siblings.** `pcgr` locates `pcgrr` as
  `dirname($CONDA_PREFIX)/pcgrr`. Built from upstream's pinned lock files,
  because PCGR is not on bioconda at all.
- **The Broad files are split across two buckets**, and asking either for the
  other's files returns 403 or 404.
- **Every download is size-verified** against the server before it is renamed
  into place. A truncated reference otherwise gets reused indefinitely by any
  check that only asks whether the file exists.
- **COSMIC contig naming**, as above.
- **MSIsensor2 ships no models.** The conda package is the binary only, so the
  installer fetches the hg38 models from upstream. Without them, tumour-only
  MSI cannot run at all — and PCGR does not cover the gap, because it restricts
  MSI prediction to WGS/WES tumour–control runs.

---

## Verify this bundle transferred intact

```bash
sha256sum -c SHA256SUMS
```

The same command works inside an installation that `--update` has refreshed:
the manifest travels with the code, so it describes what is actually there.

---

## Read before clinical use

This pipeline has been validated on **one sample, one assay, tumour-only**.
Specifically not yet established:

- **No orthogonal confirmation.** No call produced by this pipeline has been
  checked against an independent method.
- **The tumour–normal path has never processed real data.** It exists and is
  exercised by the dry-run tests, but that is not the same thing.
- **The RNA branch has never been validated at all.** It is newer than
  everything above and has no truth set: no fusion call it produces has been
  checked against FISH, RT-PCR or a second caller. For anchored-PCR panels
  (Archer, Oncomine) the vendor's own caller is the validated route, and
  results here are orthogonal evidence rather than a replacement. See
  [RNA_SCOPE.md](RNA_SCOPE.md).
- **SNV and indel only on the DNA branch.** No copy-number or
  structural-variant calling.
  Fusions are the RNA branch's job. Copy number is scoped in
  [CNV_SCOPE.md](CNV_SCOPE.md) but not implemented: it
  waits on 10+ normals sequenced on the same assay, in the same lab, at the
  same fixation state. Every tool it needs is already present or one conda
  package away.
- **Validated on one chemistry.** Everything here — QC, alignment,
  filtering thresholds — was tuned on Illumina paired-end capture data. The
  panel profiles configure what can be configured for amplicon and Ion
  Torrent chemistries, but that is not the same as having been validated on
  them. A profile is a set of sensible defaults, never a validation.
- **No UMI consensus calling.** UMIs are extracted into the read name; no
  consensus or duplex reads are built. A kit sold for 0.1% ctDNA detection
  will not reach 0.1% here.
- **A negative is only as good as its coverage.** A region the sequencing
  never reached produces no variant, exactly like a region that is wild type,
  and a VCF cannot say which happened. Pass your panel BED as `--coverage-bed`
  and the run writes a per-sample report naming every region that missed the
  depth threshold; without one, absence of a finding carries no information
  about whether the region was even looked at.
- **No tumour purity or HRD scoring**, so variant allele fractions are
  uncorrected and clonality is not claimed.
- **Tumour-only TMB is unreliable** — PCGR says so itself — and should not be
  reported without a matched normal.
- **On FFPE material**, keep step 6's read-orientation model: formalin
  deaminates cytosine, and that step is what separates the resulting C>T/G>A
  damage from biology. Residual damage still reaches PASS at low allele
  fraction, so `--min-allele-fraction` earns its place.

Treat the outputs as research-grade. Clinical deployment needs a validation
set, orthogonal confirmation, and a matched normal.

---

## Licence

MIT — see [LICENSE](LICENSE). The pipeline orchestrates third-party tools
(GATK, BWA-MEM2, fastp, samtools, bcftools, SnpEff, PCGR, Ensembl VEP) and the
reference databases it downloads, each of which carries **its own licence and
terms**. COSMIC in particular requires a licence for commercial use. The MIT
grant here covers this code, not those.
