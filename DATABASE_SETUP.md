# Database Setup — where to put COSMIC, gnomAD and the rest

Scope: which reference/annotation databases this pipeline needs, where to copy
them on this machine, and what happens when each one is absent. Companion to
`install.md` §6, which covers the download commands in full.

## Installing these automatically

`install_pipeline.py` downloads every database on this page to the layout
below, from the bucket that actually serves each one, and verifies each
download against the size the server reports:

```bash
python3 install_pipeline.py --only resources    # just the VCFs
python3 install_pipeline.py --check             # what is present
```

COSMIC is the exception — it needs a registered account, so download it
yourself and pass `--cosmic <file>`. The script renames its contigs to UCSC
style, which is the step that is silent when skipped.

The rest of this page explains the layout, what each database is for, and what
happens when one is missing.

## The short answer

There is **no hard-coded database location**. Every resource is passed as an
explicit absolute path on the command line (see the `--dbsnp`,
`--germline-resource`, `--panel-of-normals`, `--known-indels` and `--cosmic`
options in `comprehensive_variant_calling.py`), and nothing is auto-discovered.
The directory is therefore a convention you pick.

Recommended layout on this machine — `/data` does not exist and creating it
needs root, so everything lives under the home directory instead:

```text
~/data/references/hg38/     # genome FASTA + its indexes
~/data/resources/hg38/      # all VCF databases
```

Both directories already exist.

## Adding COSMIC

gnomAD was fetched fresh from the public GATK bucket and is already in place.
COSMIC is the one database still outstanding — it needs registration at
<https://cancer.sanger.ac.uk/cosmic>, so it cannot be downloaded unattended.

| File | Copy to | Passed as |
|------|---------|-----------|
| COSMIC coding mutations VCF (+ index) | `~/data/resources/hg38/` | `--cosmic` |

Two checks before copying:

- **Build must be hg38.** A hg19 COSMIC or gnomAD against an hg38 reference
  fails with contig-mismatch errors. This is the most common first-run failure.
- **The `.tbi` index must sit next to the `.vcf.gz`.** GATK and bcftools never
  search directories for indexes — they look only beside the file. COSMIC in
  particular ships plain-gzipped rather than bgzipped, so it usually needs
  re-compressing first:

```bash
zcat COSMIC_vXX_hg38.vcf.gz | bgzip -c > ~/data/resources/hg38/COSMIC_vXX_hg38.vcf.gz
tabix -p vcf ~/data/resources/hg38/COSMIC_vXX_hg38.vcf.gz
```

## Current state of this machine

**On a machine installed by `install_pipeline.py`,** the layout below is what
you get.
The conda environment `cancer_pipeline` holds every tool, and the reference
genome, its bwa-mem2 index, the SnpEff hg38 database and all public resource
VCFs (gnomAD, dbSNP, Mills, known indels, PoN) are in place under `~/data`.

PCGR 2.3.2 was added on 2026-09-07 with its 20260620 reference bundle and the
Ensembl VEP 115 GRCh38 cache. Note it runs from its **own** conda env, not
`cancer_pipeline` — the web interface switches environments for you, and
standalone it is `conda activate pcgr && python pcgr_report.py ...`.

COSMIC v103 was added on 2026-09-07. **Use the `.chr.` file** — COSMIC ships
Ensembl-style contigs (`1`, `2`, `MT`) and this pipeline needs UCSC-style
(`chr1`, `chr2`, `chrM`); the original annotates nothing, silently. See
the COSMIC section above.

`small_exac_common_3.hg38.vcf.gz` (1.3 MB) was added on 2026-09-07 for the
contamination-estimation step. Note it is **not** in
`gcp-public-data--broad-references` — that path 404s; it exists only in
`gatk-best-practices/somatic-hg38`.

Also on 2026-09-07 the environment's JDK was pinned to 21. GATK's Spark tools
cannot run on 24+, and the environment had picked up openjdk 25, which killed
duplicate marking. `java -version` is the check; `gatk --version` passes either
way. See `install.md` §7.

Nothing is outstanding. The whole pipeline was run end to end against a real
tumour-only panel, end to end through PCGR.

Free space on `/`: ~737 GB.

## What is still needed

Sizes are approximate.

### Required for a tumour-only run — all installed

1. **hg38 reference FASTA** — `~/data/references/hg38/Homo_sapiens_assembly38.fasta`
   (3.1 GB) with `.fai`, `.dict` and the five bwa-mem2 index files (16 GB).
   Kept under Broad's original filename so GATK resolves the `.dict` by
   convention. `--reference hg38` also works now — the dead S3 URL was
   repointed at Broad's public GCS bucket and the auto-download path was
   verified end to end on 2026-09-07. Passing this explicit path is still
   preferable: `--reference hg38` caches a *second* 20 GB copy of the genome
   and its indices under your output directory.
2. **SnpEff hg38 database** — 448 MB at
   `~/miniconda3/envs/cancer_pipeline/share/snpeff-5.4.0c-0/data/hg38`.
   Not under `~/data`; it lives in SnpEff's own `data.dir`.
3. **gnomAD** — `~/data/resources/hg38/af-only-gnomad.hg38.vcf.gz` (3.0 GB).
   Essential in tumour-only mode: without it Mutect2 cannot filter germline
   calls.

### Optional — each degrades gracefully if absent

| Resource | Flag | Effect if missing |
|----------|------|-------------------|
| COSMIC coding mutations (~1 GB) | `--cosmic` | Step 9 annotation and overlap counts skipped |
| `small_exac_common_3.hg38.vcf.gz` (1.3 MB) | `--contamination-resource` | Step 7 skipped; FilterMutectCalls then assumes **zero** contamination instead of measuring it |
| MSIsensor2 `models_hg38` (251 MB) | `--msi-models` | Step 8 skipped, and there is **no MSI answer at all** — PCGR restricts its own MSI prediction to WGS/WES tumour-control runs |
| dbSNP `Homo_sapiens_assembly38.dbsnp138.vcf.gz` (~10 GB) | `--dbsnp` | BQSR **silently skipped** |
| Mills + 1000G gold indels (~20 MB) | `--known-indels` | BQSR skipped (same known-sites role) |
| Panel of Normals `1000g_pon.hg38.vcf.gz` (~20 MB) | `--panel-of-normals` | No recurrent-artefact subtraction |
| PCGR data bundle + Ensembl VEP cache (tens of GB each) | `--pcgr-refdata-dir`, `--vep-dir` | PCGR step skipped with a warning |
| Panel target BED (from your assay vendor) | `--intervals` | Mutect2 walks the **whole genome** and calls on off-target reads — on a 4.5 Mb panel that was 127,081 PASS calls instead of ~1,400 |

dbSNP, Mills and the PoN all come from the Broad public buckets; the `wget`
commands are in `install.md` §6.

> **The panel of normals here is not the one copy-number analysis needs.**
> `1000g_pon.hg38.vcf.gz` is sites-only — zero samples, no coverage — a list of
> recurrent artefact positions for Mutect2. CNV needs a *read-count* panel of
> normals built from BAMs on your own assay. See `CNV_SCOPE.md`.

## Actual layout on disk

```text
~/data/references/hg38/                            20 GB
├── Homo_sapiens_assembly38.fasta                  3.1 GB   FASTA
├── Homo_sapiens_assembly38.fasta.fai              158 KB   samtools
├── Homo_sapiens_assembly38.dict                   569 KB   GATK
├── Homo_sapiens_assembly38.fasta.bwt.2bit.64      9.8 GB   bwa-mem2
├── Homo_sapiens_assembly38.fasta.0123             6.0 GB   bwa-mem2
├── Homo_sapiens_assembly38.fasta.pac              768 MB   bwa-mem2
└── Homo_sapiens_assembly38.fasta.amb / .ann       465 KB   bwa-mem2

~/data/resources/hg38/                             6.4 GB
├── af-only-gnomad.hg38.vcf.gz                     3.0 GB  (+ .tbi)
├── Homo_sapiens_assembly38.dbsnp138.vcf.gz        1.5 GB  (+ .tbi)
├── Homo_sapiens_assembly38.known_indels.vcf.gz     59 MB  (+ .tbi)
├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 20 MB (+ .tbi)
├── 1000g_pon.hg38.vcf.gz                           17 MB  (+ .tbi)
├── Cosmic_GenomeScreensMutant_v103_GRCh38.chr.vcf.gz 870 MB (+ .tbi)  <- USE THIS
└── Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz  904 MB  original, Ensembl contigs
    Cosmic_clean.vcf.gz                            114 MB  INFO stripped; unused
└── small_exac_common_3.hg38.vcf.gz                1.3 MB  (+ .tbi)  contamination sites

~/data/msisensor2/models_hg38/                     251 MB   MSIsensor2 models (2,829 files)

~/data/pcgr/20260620/data/grch38/                  7.3 GB   PCGR reference bundle
~/data/vep_cache/homo_sapiens/115_GRCh38/           24 GB   Ensembl VEP cache
```

PCGR's two paths are passed as `--pcgr-refdata-dir ~/data/pcgr/20260620` and
`--vep-dir ~/data/vep_cache` (the **parent** of `homo_sapiens/`), and only from
a shell with the `pcgr` env activated.

## Example invocation

```bash
conda activate cancer_pipeline

python comprehensive_variant_calling.py \
    -o results/ \
    --reference ~/data/references/hg38/Homo_sapiens_assembly38.fasta \
    --germline-resource ~/data/resources/hg38/af-only-gnomad.hg38.vcf.gz \
    --dbsnp ~/data/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-indels ~/data/resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                   ~/data/resources/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --panel-of-normals ~/data/resources/hg38/1000g_pon.hg38.vcf.gz \
    --contamination-resource ~/data/resources/hg38/small_exac_common_3.hg38.vcf.gz \
    --msi-models ~/data/msisensor2/models_hg38 \
    --cosmic ~/data/resources/hg38/Cosmic_GenomeScreensMutant_v103_GRCh38.chr.vcf.gz \
    --threads 32
```

(The previous version of this example was missing a line continuation after
`--threads 32`, which silently split it into two commands and dropped
`--cosmic`.)

**On a capture panel or exome, add the target restriction.** It matters more
than every optional database above combined:

```bash
    --intervals /path/to/panel_targets.bed \
    --interval-padding 50 \
    --min-depth 20
```

`--min-allele-fraction` is deliberately **not** shown with a value here. It is
type-blind, so it buys indel cleanup with SNV sensitivity, and it must come
from the assay's validated limit of detection rather than from whatever tidies
the output. The foldback filter below is the targeted tool for artefactual
insertions and costs no SNVs at all.

Measured on a 4.5 Mb panel: 635,414 raw calls and 127,081 PASS without them,
20,761 raw and ~1,400 PASS with. Use the BED your assay vendor ships.

**Foldback insertions are filtered by default** (`--foldback-min-match`, 0.7).
A fragment end that folds back on itself during library prep is sequenced as a
sequence followed by its own reverse complement, and aligns as an insertion of
the reverse complement of the adjacent reference. On the panel above, 1,882 of
2,047 insertions were exactly that, and they accounted for 206 of 270 PASS
frameshift/stop calls — including BRCA1, BRCA2, PALB2, PTEN and ATM
frameshifts that PCGR tiered as oncogenic. They are tagged, never removed, and
`INFO/FBMATCH` records the matched fraction for every insertion tested.

## Using the web interface

`webapp/app.py` fills in every resource on this page automatically. Anything
left blank on the form is set at submit time from what is installed under
`~/data`, so a bare submission still runs the full flow and ends in a clinical
report. The job page lists what was filled in; a **minimal run** checkbox opts
out.

The resources directory must still sit inside an `--allow-root` path, or the
form rejects those paths — including the ones it filled in itself:

```bash
python webapp/app.py --allow-root ~/data --allow-root /path/to/fastqs
```

Two things it cannot supply: the target BED, which is your assay vendor's
file and must also be under an `--allow-root` path, and a matched normal.

For **several samples**, tick **Batch** and the app runs each as its own
patient in turn, each with its own output subdirectory and report. It refuses
a multi-sample directory that is not marked either Batch or tumour/normal,
because `--auto-discover` left to itself would treat the second sample as the
first one's matched normal.

No conda environment needs activating first — the app runs the pipeline in
`cancer_pipeline` and PCGR in `pcgr` itself.
