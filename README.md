# Cancer DNA Pipeline — installer bundle

Everything needed to stand this pipeline up on a new machine. Copy this whole
directory across, run one command, and you get the analysis environment, the
reference genome, six resource databases, the clinical reporter and its data.

Nothing here is patient data. This bundle is code, configuration and
documentation only.

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

`python3 install_pipeline.py --check` doubles as a health check afterwards, and
is the fastest way to tell whether an environment has drifted.

### Staying current

```bash
python3 check_db_updates.py --print       # any newer database releases?
```

Writes its answer to `~/.cache/cancer_pipeline/db_updates.json`. The script
**reports only** — nothing is downloaded or replaced, and no run reads it.

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
pipeline_orchestrator.py          chains the three stages
pcgr_report.py                    clinical report wrapper
coverage_report.py                target coverage: which regions a negative
                                  result is actually entitled to speak for
cancer-dna-pipeline.1             man page

webapp/                           Flask front end (patient details, batch, cleanup)

install.md                        manual install, and the reference the installer follows
DATABASE_SETUP.md                 what each database is for, and where it goes
CNV_SCOPE.md                      copy-number: designed, not built — read before starting it

check_db_updates.py               tells you when a database has a newer release
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
- **SNV and indel only.** No copy-number, structural-variant or fusion calling.
  Copy number is scoped in [CNV_SCOPE.md](CNV_SCOPE.md) but not implemented: it
  waits on 10+ normals sequenced on the same assay, in the same lab, at the
  same fixation state. Every tool it needs is already present or one conda
  package away.
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
