<!-- Created by Brainstorm, 2026. -->

# Known limitations

**Research use only.** This is the register of everything this pipeline is
known *not* to do, not to have proven, or to do in a way that can mislead.
It is kept in one place so that nothing on it depends on a reader finding
the right paragraph in the right document. Each entry has a stable ID; the
other documents refer to these IDs where the limitation applies.

**Status** means:

| Status | Meaning |
|---|---|
| **Open** | A gap in this pipeline that could be closed |
| **By design** | Deliberate: closing it would mean guessing on the operator's behalf |
| **Third-party** | Behaviour of a tool this pipeline runs, documented rather than patched |
| **Inherent** | A property of the method itself; no pipeline can remove it |
| **Not built** | Out of scope today, with the reason recorded |

Last reviewed: 23 September 2026, after every pathway was run end to end on
known-answer data (see [validation/](validation/README.md)).

---

## Validation — what has and has not been proven

### L-01 · Proven on simulated reads only · Open

Every pathway (web, command line, orchestrator, hybrid, worksheet; DNA and
RNA) finds exactly the planted variants and fusions in **simulated** reads
from real hg38. Simulated reads have no formalin damage, no PCR artefacts,
no capture bias, no uneven coverage and no real sequencing error profile.
Passing proves the pipeline finds what is there; it does not measure
sensitivity or specificity on real specimens.

*What to do:* validate on real specimens with orthogonally confirmed
results before any clinical use.

### L-02 · No orthogonal confirmation · Open

No call from this pipeline — SNV, indel or fusion — has been checked
against an independent method (ddPCR, Sanger, FISH, RT-PCR, a second
caller) on a real specimen.

### L-03 · Validated on one chemistry · Open

Thresholds were tuned on Illumina paired-end hybrid-capture data. Amplicon
and Ion Torrent profiles configure what can be configured; that is not
validation. Single-end RNA libraries are accepted by the RNA engine but no
single-end run has been made. See [PANELS.md](PANELS.md).

### L-04 · GRCh37 / hg19 never run · Open

The code has an hg19 path (SnpEff database and PCGR build are switched from
the reference name), but every run to date — including all validation — is
hg38.

### L-05 · The validation kit is not automatic · Open

`validation/` needs the installed reference data and about an hour of
machine time for every pathway, so it is run by hand, not by
`run_tests.py`. There is no continuous integration of any kind; the unit
suite is also run by hand.

*What to do:* run `validation/make_truth.sh` and `check_truth.py` after any
change to a tool version, a reference, or a pipeline step.

---

## Third-party tool behaviour

### L-06 · SnpEff's first annotation is not the clinical transcript · Third-party

SnpEff lists every overlapping transcript, ordered by effect severity. On a
known BRAF V600E its **first** `ANN` entry read `p.Val640Glu` — the same
variant numbered on a longer isoform (`NM_001374258`). The clinical
numbering is present further along the list, and the PCGR report (VEP,
MANE transcript) correctly says `p.Val600Glu`.

*What to do:* take protein nomenclature from the PCGR report, never from
the first `ANN` entry of `annotated/<sample>.annotated.vcf`. Nothing in this
pipeline reads `ANN`. Making SnpEff prefer MANE needs `-canonList` and a
MANE transcript list that the installer does not fetch; plain `-canon` picks
the longest CDS, which for BRAF is the wrong one.

### L-07 · PCGR 2.3.2 under-tiers fusions · Third-party

With the tumour site set to lung — confirmed to have reached PCGR, which
recorded `primary_site: Lung` — a known **EML4::ALK was tiered 2**
("potential significance"), with no primary-site match recorded; BCR::ABL1
likewise. The SNVs in the same report were site-matched and tiered 1
normally. EML4::ALK in NSCLC is a tier-1 finding with approved therapy.

*What to do:* read a fusion's tier as a floor. Interpret canonical driver
fusions against current guidelines, not against PCGR's tier alone.

### L-08 · PCGR's MSI and TMB on panels · Third-party

PCGR omits its MSI classifier on a TARGETED assay; MSIsensor2 (step 8) is
the only MSI answer on a panel. PCGR's TMB is computed over unfiltered calls
unless depth/AF tags are lifted (`--pcgr-lift-tags`), and PCGR itself warns
that tumour-only TMB is unreliable.

---

## The command line and the web interface differ

### L-09 · The command line does not find installed resources · Open

The web form pre-fills every installed resource (dbSNP, known indels,
gnomAD, panel of normals, contamination sites, COSMIC, MSI models, PCGR
bundle, VEP cache). The command-line engine does not: without the flags,
**BQSR, contamination, MSI, COSMIC and PCGR are skipped**, each with a
`[SKIP]` line, and the run exits 0. Only the reference is found by name
(`-r hg38`).

*What to do:* pass the resources explicitly (the man page and README list
them), or use the web interface. Check the `[SKIP]` lines of every
command-line run.

### L-10 · QC without `--merge-lanes` treats each lane as a sample · By design

`fastq_qc_clean.py` (stage 1 of the orchestrator) processes each lane as a
separate sample unless `--merge-lanes` is given, and warns when it sees
lane-split input. The engines' own auto-discovery always merges lanes, and
an explicit single lane of a multi-lane sample is refused.

*What to do:* pass `--merge-lanes` to stage 1 for lane-split data.

### L-11 · The worksheet has no matched-normal columns · Not built

Worksheet runs are tumour-only. A tumour/normal pair has to be run through
the single-sample DNA form (or the command line).

---

## Configuration the operator must supply

### L-12 · Sample roles must be named when several samples are present · By design

With more than one sample in a directory or QC manifest, the tumour must be
named (`--tumour-sample`) and the normal is only ever one named
(`--normal-sample`), or the run is declared `--tumour-only`. Nothing is
taken from file order: that rule once analysed the normal as the tumour for
ordinary names such as PT01-N / PT01-T. A run that does not name them stops
with the list of samples it found.

### L-13 · RNA depth floors come from a profile or the operator · By design

Whether a library is deep enough to support a negative fusion result is
assay-specific. Without an RNA panel profile (or `--rna-min-reads-millions`
/ `--rna-min-unique-mapped-pct`), library size and unique mapping are shown
as **not judged**, and the QC verdict refuses to call an empty fusion table
a negative. See [RNA.md](RNA.md).

*What to do:* choose an RNA panel profile for real runs.

### L-14 · Actionability depends on the tumour site · By design

PCGR tiers against the tumour site. Site 0 ("unspecified") gives tiers that
are not specific to the patient's tumour. The worksheet sets it per row for
this reason.

### L-15 · No coverage statement without a target BED · By design

Without a panel BED, step 12 makes no coverage statement and the TMB
denominator cannot be measured. An absent variant then cannot be told apart
from a region that was never sequenced deeply enough.

---

## Installation state

### L-16 · The STAR index can fall behind the GENCODE release · Open

The STAR index bakes in the annotation it was built from. When the
installer's default GENCODE release changes, the new GTF is downloaded but
the existing index is **not** rebuilt (that costs about an hour and ~32 GB
of RAM). Runs stay internally consistent — the web form offers the index's
own annotation — and `install_pipeline.py --check` now names the mismatch.

*What to do:* when `--check` reports it, rebuild deliberately with
`install_pipeline.py --only star-index --force`, or keep the older release.
(On the development machine, as of this review: index built from v44,
installer default v50.)

---

## Inherent to the method

### L-17 · Tumour-only mode cannot remove novel germline variants · Inherent

Without a matched normal, a germline variant absent from gnomAD is
indistinguishable from a somatic one. In validation, a heterozygous germline
SNP was reported at VAF ~0.50 in every tumour-only run, and correctly
removed in every matched-normal run. A VAF near 0.5 or 1.0 is the usual tell.

### L-18 · Panel TMB below ~1 Mb is not interpretable · Inherent

The denominator is measured from the BED, which makes the figure honest,
not meaningful: on a small panel one artefact moves TMB by several
mutations per Mb. The pipeline warns; do not report it.

### L-19 · A negative is only as good as its coverage and its library · Inherent

See L-13 and L-15. An empty result says nothing about regions the sequencing
did not reach, or about an RNA library too shallow or degraded to sample a
fusion.

### L-20 · The chimeric-output guard's failure mode is unproven · Open

The RNA engine compares STAR's chimeric-read count with the chimeric
alignments in the BAM, to catch an empty fusion table caused by STAR not
writing them. The guard is tested against real BAMs; the failure it guards
against has **not** been observed. (An earlier report that it had been came
from a broken detector, since fixed.)

---

## Not built

### L-21 · No copy-number or structural-variant calling · Not built

SNV and indel only on DNA; fusions only on RNA. Copy number waits on 10+
normals sequenced on the same assay — see [CNV_SCOPE.md](CNV_SCOPE.md).

### L-22 · No UMI consensus · Not built

UMIs are extracted into the read name; no consensus or duplex reads are
built. A kit sold for 0.1% ctDNA detection will not reach 0.1% here.

### L-23 · No purity, clonality, HRD, FLT3-ITD or germline predisposition calling · Not built

Allele fractions are uncorrected for purity; clonality is not claimed; HRD
scores, FLT3 internal tandem duplications (which short-read callers
routinely miss) and germline predisposition reporting are not produced.

---

## Fixed defects

The defects found and fixed are not limitations and are not listed here.
Each is recorded, with what a user would have seen, in
[TESTING.md](TESTING.md) and pinned by a named test in
`tests/test_regressions.py`.
