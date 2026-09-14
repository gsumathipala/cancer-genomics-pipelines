<!-- Created by Brainstorm, 2026. -->
# Panels — running more than one assay

This pipeline was written against one panel. Everything that panel needed
is a separate flag with a general-purpose default, which is fine until the
laboratory runs a second kit — at which point the operator has to know,
from memory and per run, that

- an **amplicon** kit must not have its duplicates marked and a **capture**
  kit must;
- a capture kit wants 100 bp of interval padding and an amplicon kit wants
  none;
- **BQSR** has too few covered sites to fit on a 20 kb hotspot panel;
- **TMB** divided by PCGR's assumed 34 Mb is wrong for both of them, by
  different factors;
- the **QIAseq** run in the next directory carries its UMI in the first 12
  bases of read 2.

Every one of those is silent when it is wrong. The run finishes, the report
renders, and the numbers are confidently incorrect.

A **panel profile** turns that into data. One name sets all of it, prints
what it set, and records it in the run manifest.

```bash
python comprehensive_variant_calling.py \
    --panel thermo-oncomine-cav3 \
    --panel-bed /data/panels/OCAv3.designed.bed \
    -i fastqs/ -o results/ -r hg38 --auto-discover
```

---

## The three commands to start with

```bash
python comprehensive_variant_calling.py --list-panels
python comprehensive_variant_calling.py --describe-panel qiaseq
python comprehensive_variant_calling.py --new-panel-template mine.json
```

None of them needs an output directory, a reference or any data.

In the web interface the same thing is a dropdown on the new-run form, and
the **Panels** page lists every profile with what it sets and what it warns
about. The choice is applied **server-side** at submit time, not by the
browser: a value that only ever existed in a page can go missing when
somebody posts from a script or leaves a tab open while a profile changes.
The run page then names the assay and lists which fields the profile filled
in, and the PDF report carries both into its provenance table and its
limitations.

---

## What a profile is, and is not

A profile is **settings, not content**. It carries the numbers, the
chemistry-driven step choices and the UMI layout. It does **not** carry the
target BED, because

- vendor BEDs are licensed content that travels with the kit, not something
  this bundle may redistribute;
- they are specific to a kit version *and* a genome build, and the wrong one
  is worse than none — it looks like a target region and silently is not;
- a laboratory's spike-ins mean the BED in the freezer room is frequently
  not the one on the vendor's website.

So you always supply the BED:

```bash
--panel-bed /data/panels/kit.bed     # sets --intervals AND --coverage-bed
```

They are the same file in all but unusual laboratories, and filling only
one of them is how a run ends up either calling genome-wide or making no
coverage statement at all.

## Precedence

```
explicit command-line flag   >   panel profile   >   script default
```

A profile only ever fills in what the command line did not state, so adding
`--panel` to a command that already works cannot change it. What was
applied — and what the profile *wanted* to apply and was not allowed to —
is printed in the first few lines of the run:

```
Panel   : Oncomine Comprehensive Assay v3 [thermo-oncomine-cav3]
          Thermo Fisher, amplicon, from built-in
          settings applied from this profile:
            --interval-padding 0
            --min-allele-fraction 0.05
            --min-depth 100
            --skip-steps dedup
          NOT applied -- your command line said otherwise:
            --min-allele-fraction (profile suggested 0.05)
          NOTE: Ion Torrent reads are single-end and have a different indel
                error profile from Illumina reads ...
```

---

## What the built-in profiles cover

| Generic shape | Use it for |
|---|---|
| `generic-capture` | any probe-capture panel not listed below |
| `generic-capture-umi` | the same, with molecular barcodes still in the reads |
| `generic-amplicon` | any multiplex-PCR panel |
| `generic-hotspot` | a small amplicon panel (< ~100 kb) |
| `ctdna-capture` | a deep plasma panel |
| `wes` | whole exome, any vendor |
| `wgs` | whole genome (no target restriction) |

| Named kit | Manufacturer | Chemistry |
|---|---|---|
| `illumina-tso500` | Illumina | hybrid capture |
| `illumina-tst170` | Illumina | hybrid capture |
| `illumina-ampliseq-cancer-hotspot` | Illumina | amplicon |
| `thermo-oncomine-cav3` | Thermo Fisher | amplicon |
| `thermo-oncomine-focus` | Thermo Fisher | amplicon |
| `thermo-ion-ampliseq-chpv2` | Thermo Fisher | amplicon |
| `agilent-sureselect-xths2` | Agilent | hybrid capture |
| `agilent-sureselect-exome` | Agilent | hybrid capture |
| `twist-custom-capture` | Twist Bioscience | hybrid capture |
| `twist-exome` | Twist Bioscience | hybrid capture |
| `idt-xgen-pan-cancer` | IDT | hybrid capture |
| `idt-xgen-exome` | IDT | hybrid capture |
| `qiagen-qiaseq-dna` | QIAGEN | amplicon (SPE + UMI) |
| `archer-variantplus` | Invitae / ArcherDX | amplicon (AMP + UMI) |
| `roche-kapa-target` | Roche | hybrid capture |

**Built-ins are a starting point, not an authority.** They encode published,
nominal figures for *one version* of each kit. Kits are revised,
laboratories spike in extra content, and a validated limit of detection is a
property of your validation rather than of a catalogue number. Copy the
closest one, adjust it, keep the copy.

---

## Adding your own kit

```bash
# Start from the closest built-in:
python comprehensive_variant_calling.py \
    --panel generic-capture --new-panel-template our-lung-panel.json

# Edit it, then use it:
python comprehensive_variant_calling.py \
    --panel-file our-lung-panel.json --panel our-lung-panel \
    --panel-bed /data/panels/lung.bed -i fastqs/ -o results/ -r hg38 ...
```

Drop the file in `~/.config/cancer_pipeline/panels/` and `--panel-file`
becomes unnecessary — it is found automatically, as is anything in a
directory named by `$CANCER_PIPELINE_PANEL_DIR`. **A profile whose `id`
matches a built-in replaces it**, which is how a site corrects a built-in
for a local kit revision without editing the bundle and losing the change
at the next update.

The other direction works too: tune a run with ordinary flags under
`--dry-run` until it is right, then save what you arrived at.

```bash
python comprehensive_variant_calling.py ... --dry-run \
    --min-allele-fraction 0.02 --interval-padding 75 \
    --save-panel-as our-lung-panel.json --panel-id our-lung-panel
```

Paths are deliberately not saved: a profile carrying a BED path would
configure a different target, or none, on any other machine.

### The file format

```json
{
  "id": "our-lung-panel",
  "name": "Lung 42-gene capture panel",
  "manufacturer": "Twist Bioscience",
  "chemistry": "hybrid-capture",
  "aliases": ["lung42"],
  "genome_builds": ["hg38"],
  "nominal_target_size_mb": 0.42,
  "summary": "In-house 42-gene lung panel, v3 design.",
  "requires": ["panel BED"],
  "notes": ["Validated LoD 3% on FFPE; see validation report LV-2024-07."],
  "settings": {
    "interval_padding": 100,
    "min_allele_fraction": 0.03,
    "min_depth": 100,
    "pcgr_assay": "TARGETED",
    "pcgr_estimate_tmb": false
  }
}
```

`chemistry` is one of `hybrid-capture`, `amplicon`, `wgs` and is not
cosmetic: it is what decides whether duplicate marking is meaningful.
Setting keys are named after the pipeline's own options with dashes
replaced by underscores, and **an unknown key is an error, not a silently
ignored line** — a typo that quietly does nothing is the failure this
whole mechanism exists to prevent. One file may hold one profile, a list of
them, or `{"panels": [ ... ]}`.

---

## What a profile does for you beyond the settings

Two things happen automatically once a target BED is known, whether or not
a profile was named:

**The footprint is measured and used as the TMB denominator.** PCGR divides
by an assumed 34 Mb for a `TARGETED` assay unless told otherwise, so a 2 Mb
panel reports a TMB about 17 times too low and says nothing about having
done so. The BED is merged (vendor BEDs routinely overlap their probe
tiles) and the real footprint is used instead. If it disagrees with the
profile's nominal figure by more than 25%, the run says so — that almost
always means the wrong kit version or the wrong genome build. Below ~1 Mb
the run warns that a panel TMB is not interpretable at all.

**The contig naming is checked and fixed.** Vendor BEDs arrive in whichever
convention their designer used: Ensembl writes `1` and `MT`, UCSC and the
hg38 this pipeline installs write `chr1` and `chrM`. A mismatched BED makes
GATK abort, which is at least loud — but the *same* file handed to the
coverage check simply finds no reads in any region, and the report then
states that none of the panel was covered. That is a confident, entirely
wrong clinical statement. The run detects the mismatch and writes a renamed
copy into `<output-dir>/panel/`; your own copy of the vendor document is
never touched.

---

## RNA profiles run a different pipeline

A profile whose `chemistry` begins `rna-` describes a **transcriptome**
library and routes to the RNA branch — `fusion_calling.py`, not the somatic
caller. `--list-panels` marks them `[RNA]`.

```bash
python fusion_calling.py --panel archer-fusionplex ...
```

They carry a different vocabulary: no interval padding, no allele-fraction
floor, no duplicate-marking decision (there is no dedup step to switch off —
on RNA, two reads at one position in a highly expressed gene are two
observations of an abundant transcript). What replaces them is a library-size
floor and a mapping-rate floor, because "was there enough usable library?" is
what a negative fusion result depends on.

Each engine **refuses the other's profiles by name** rather than running with
settings that have no meaning, and so do the two web forms. See
[RNA.md](RNA.md).

## Hybrid kits: one specimen, two profiles

A profile may name the other half of its kit with `pairs_with`:

```json
{ "id": "illumina-tso500", "pairs_with": "illumina-tso500-rna", ... }
```

Built in for TSO500, Oncomine Comprehensive v3 and Archer VariantPlus. The
web interface's **hybrid** pathway reads it: choose the kit and both halves
are filled in, both runs are queued, and the RNA half is linked to the DNA
half so PCGR reports the specimen as one result. Add `pairs_with` to your
own profiles and your kit appears there too.

## The standalone tools take it too

`--panel` is not only a pipeline flag. The two tools that can be run on their
own accept it, so a report regenerated by hand is configured the same way the
run was:

```bash
# The coverage threshold is an assay property: a hotspot panel at 2000x and
# an exome at 100x cannot share one.
python coverage_report.py --panel generic-hotspot     --bam sample.bqsr.bam --bed panel.bed --sample-id S1 --output-dir cov/

# The TMB denominator, measured rather than assumed.
python pcgr_report.py --panel illumina-tso500 --panel-bed panel.bed     --input-vcf calls.vcf.gz --output-dir pcgr/ --sample-id S1     --pcgr-refdata-dir ~/data/pcgr/... --vep-dir ~/data/vep_cache
```

Each spells its options without the prefix the pipeline needs to disambiguate
them — the pipeline's `--coverage-min-depth` is `coverage_report.py`'s
`--min-depth`, and its `--pcgr-assay` is `pcgr_report.py`'s `--assay`. Profiles
are written in *one* vocabulary regardless, and the mapping is printed when a
setting is applied, so you can see both the flag and the key you would edit:

```
    --min-depth (profile key: coverage_min_depth) 100
```

`align_reads.py` has no `--panel`, because nothing in alignment depends on the
kit. Stage 2 is the one stage a profile has nothing to say about.

## Chaining all three stages

Stage 1 (`fastq_qc_clean.py`) and stage 3 both understand `--panel`, and
both must be given the same one — a chain whose QC extracted a UMI and
whose caller was configured for a different chemistry is half-configured,
and the half that is wrong is silent. The orchestrator takes it once:

```bash
python pipeline_orchestrator.py \
    --panel qiagen-qiaseq-dna --panel-bed /data/panels/qiaseq.bed \
    --stage1-args "-i raw/ -o qc/ --threads 16" \
    --skip-stage2 \
    --stage3-args "-o results/ --reference hg38 --threads 16"
```

---

## What a profile cannot fix

**There is no UMI consensus step.** fastp moves a UMI into the read name,
which is enough for auditing and for UMI-aware deduplication downstream,
but it is not error-corrected consensus calling — there is no fgbio, no
duplex collapsing. A kit sold for 0.1% ctDNA detection will not reach 0.1%
here, and the `ctdna-capture` profile's floor of 0.5% says so rather than
claiming a sensitivity this code does not have.

**Primer trimming is the vendor's job.** This pipeline does not trim
amplicon primers. An untrimmed primer reports the oligo's sequence as the
patient's.

**Ion Torrent data is not what any of this was tuned on.** The Thermo
profiles configure what can be configured, but the QC, alignment and
filtering defaults throughout were tuned on Illumina paired-end reads, and
the indel error profile in homopolymers differs. Validate before trusting
indel calls from an Ion run.

**A profile is not a validation.** The numbers in a built-in are
conventional defaults, not your assay's measured performance. Nothing here
replaces a validation set, orthogonal confirmation and a matched normal —
see the last section of the README.
