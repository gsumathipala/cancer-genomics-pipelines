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

Change the data location with `--data-dir`. The conda environments go wherever
conda keeps its environments.

---

## Contents of this bundle

```
install_pipeline.py               the installer
environment.yml                   the analysis environment definition

comprehensive_variant_calling.py  the engine — all 12 steps, self-contained
align_reads.py                    stage 2 standalone
fastq_qc_clean.py                 stage 1 standalone
pipeline_orchestrator.py          chains the three stages
pcgr_report.py                    clinical report wrapper
cancer-dna-pipeline.1             man page

webapp/                           Flask front end (patient details, batch, cleanup)

install.md                        manual install, and the reference the installer follows
DATABASE_SETUP.md                 what each database is for, and where it goes
```

**Keep the Python files in one directory.** The orchestrator locates its
siblings relative to itself, and the variant caller imports `pcgr_report` from
alongside it. The flat layout is deliberate.

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

---

## Read before clinical use

This pipeline has been validated on **one sample, one assay, tumour-only**.
Specifically not yet established:

- **No orthogonal confirmation.** No call produced by this pipeline has been
  checked against an independent method.
- **The tumour–normal path has never processed real data.** It exists and is
  exercised by the dry-run tests, but that is not the same thing.
- **SNV and indel only.** No copy-number, structural-variant or fusion calling.
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
