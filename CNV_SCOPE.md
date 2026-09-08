# Somatic copy number: scope for when normals exist

Status: **not implemented.** This is the design to build against, written
2026-09-08 so the decisions are settled before anyone starts.

The blocker is data, not software. Every tool needed is either already
installed or one conda package away; what is missing is a set of normal
samples from this assay. This document says how many, what they must have in
common, what the workflow does with them, and — importantly — what CNV will
and will not deliver on a 4.17 Mb panel.

---

## 1. Why it does not exist yet

The pipeline already passes a **panel of normals** to Mutect2:

```
--panel-of-normals ~/data/resources/hg38/1000g_pon.hg38.vcf.gz
```

That is a different artefact from the one CNV needs, despite the shared name.

| | Mutect2 PoN (have) | CNV PoN (need) |
|---|---|---|
| Contains | variant positions, sites-only | read counts per interval |
| Samples in file | **0** | one column per normal |
| Built by | `CreateSomaticPanelOfNormals` | `CreateReadCountPanelOfNormals` |
| Answers | "is this site a recurrent artefact?" | "is this region's depth unusual?" |
| Format | VCF | HDF5 count matrix |

`1000g_pon.hg38.vcf.gz` carries no coverage information at all, so there is
nothing for a copy-ratio model to normalise against.

A *public* read-count PoN would not help either. Coverage in a capture assay is
dominated by probe efficiency — which probes pull down well, which sit in
GC-rich sequence, how a particular lab's hybridisation performed. Cancelling
that pattern is precisely what the PoN does, and the pattern belongs to the
capture kit and the protocol. Normalise this panel against someone else's
coverage and every well-captured region reads as amplified, every poorly
captured one as deleted.

---

## 2. What has to be collected

**Normal samples sequenced on this assay, in this lab.**

| Requirement | Why |
|---|---|
| **10 minimum, 30–40 preferred** | GATK denoises by PCA over the count matrix; too few normals and the principal components fit noise |
| **Same capture kit and version** | The PoN cancels probe efficiency; a different kit has a different pattern |
| **Same library prep and sequencer** | Fragment-length and GC bias travel with chemistry |
| **Same fixation state — FFPE normals for FFPE tumours** | Fixation changes coverage uniformity. A fresh-frozen PoN will not model an FFPE tumour's noise, and the residual shows up as false segments |
| **Both sexes, or sex-matched** | chrX and chrY copy number is otherwise systematically wrong |
| **Same reference build (GRCh38)** | Intervals must line up |

Germline normals from the *same patients* are ideal but not required for the
PoN itself; the PoN is a coverage baseline, not a per-patient comparator.

---

## 3. The workflow

### 3a. Build the panel of normals — once per assay

```bash
# Bin the capture targets. For a panel, one bin per target: --bin-length 0.
gatk PreprocessIntervals \
    -R $REF -L panel_targets.bed \
    --bin-length 0 --padding 250 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O targets.preprocessed.interval_list

# Per normal
gatk CollectReadCounts \
    -I normal_NN.bqsr.bam -L targets.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O normal_NN.counts.hdf5

# Once, over all of them
gatk CreateReadCountPanelOfNormals \
    -I normal_01.counts.hdf5 -I normal_02.counts.hdf5 ... \
    --minimum-interval-median-percentile 5.0 \
    -O cnv.pon.hdf5
```

Rebuild the PoN whenever the capture kit, chemistry or fixation protocol
changes. It is not portable across assays.

### 3b. Per tumour

```bash
gatk CollectReadCounts -I tumour.bqsr.bam -L targets.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O tumour.counts.hdf5

gatk DenoiseReadCounts -I tumour.counts.hdf5 --count-panel-of-normals cnv.pon.hdf5 \
    --standardized-copy-ratios tumour.standardizedCR.tsv \
    --denoised-copy-ratios tumour.denoisedCR.tsv

# Allelic counts sharpen segmentation and give the minor-allele fraction.
# small_exac_common_3 is ALREADY INSTALLED for the contamination step.
gatk CollectAllelicCounts -I tumour.bqsr.bam -R $REF \
    -L ~/data/resources/hg38/small_exac_common_3.hg38.vcf.gz \
    -O tumour.allelicCounts.tsv

gatk ModelSegments \
    --denoised-copy-ratios tumour.denoisedCR.tsv \
    --allelic-counts tumour.allelicCounts.tsv \
    -O . --output-prefix tumour

gatk CallCopyRatioSegments -I tumour.cr.seg -O tumour.called.seg
```

Runtime on a 4 Mb panel is minutes, not hours. Output is a few MB.

---

## 4. The catch: PCGR wants allele-specific integer copy number

PCGR's `--input_cna` is validated against these columns
(`pcgr/validate.py`):

```
Chromosome   Start   End   nMajor   nMinor
```

`nMajor`/`nMinor` are **integer copies of the major and minor allele**. GATK's
`ModelSegments` does not produce them. It produces a *relative* log2 copy ratio
and a minor-allele fraction, which cannot be converted to absolute integer
copies without knowing **tumour purity and ploidy**.

So GATK alone gets you segments and plots, but not a PCGR-ready file. Three
ways to close that gap:

| Path | Needs | Gives | Fit here |
|---|---|---|---|
| **PureCN** on GATK output | the same PoN, plus a normal DB and the tumour VCF | purity, ploidy, `nMajor`/`nMinor` | **Best fit** — designed for tumour-only hybrid capture |
| **FACETS** or **ASCAT** | a *matched normal per tumour* | purity, ploidy, allele-specific CN | Cleanest results, but needs paired sequencing |
| **GATK only** | nothing further | relative segments, plots, gene-level gain/loss | Useful internally; not PCGR-ready |

Recommendation: **GATK for segmentation, PureCN for purity/ploidy and the
allele-specific conversion.** PureCN reuses the same PoN, and it would also
close the tumour purity gap that currently blocks clonality interpretation and
inflates tumour-only TMB.

---

## 5. Where it slots into the pipeline

CNV is a BAM-level analysis like contamination and MSI, so it belongs beside
them, after BQSR and before filtering.

| | |
|---|---|
| Step number | **9**, pushing filtering→10 through statistics→14 |
| Skip name | `cnv` |
| Failure mode | enrichment — warn and continue, never fatal |
| Enabled by | `--cnv-pon <cnv.pon.hdf5>`; absent, the step skips with a note |
| Also needs | `--intervals` (the same BED), and the existing `--contamination-resource` for allelic counts |
| Outputs | `cnv/<sample>.called.seg`, `.modelFinal.seg`, denoised plots, and `<sample>.cna.tsv` for PCGR |
| Manifest | segment count, ploidy, purity, fraction of genome altered |
| Passes to PCGR | `--input_cna <sample>.cna.tsv` |

Follows the pattern already used by `run_msisensor2()`: a single function
returning a dict or `None`, wired into `main()` behind `should_run("cnv")`,
never raising.

---

## 6. What CNV will and will not give you on this panel

**Will:** focal gene-level amplification and deletion inside the 4.17 Mb
target — *ERBB2*, *MET*, *EGFR*, *MYC*, *CDKN2A* and the rest of the captured
gene set. That is the clinically actionable majority of panel CNV, and PCGR
tiers it.

**Will not:**

- **Genome-wide ploidy or whole-arm events.** 4.17 Mb over 10,461 intervals is
  ~0.13% of the genome, scattered across gene bodies. Arm-level and
  chromosome-level calls need broad, even coverage.
- **HRD scoring.** LOH, telomeric allelic imbalance and large-scale state
  transitions are counted genome-wide over dense SNP coverage. A gene panel
  cannot supply that, whatever CNV caller sits on top. **HRD stays out of reach
  until WES or WGS**, and a "HRD score" derived from this panel would be
  meaningless rather than approximate.
- **Reliable CNV on FFPE without FFPE normals.** Fragmentation raises coverage
  variance; the PoN must carry the same noise to cancel it.

---

## 7. Risks

| Risk | Mitigation |
|---|---|
| Too few normals — PoN fits noise | Refuse to build below 10; state the count in the manifest |
| Fresh-frozen PoN against FFPE tumours | Match fixation state; record it in the PoN metadata |
| Panel too small for confident segmentation | Report segment count and interval support; flag thin segments the way MSI flags thin site counts |
| Purity estimate unstable at low tumour content | PureCN reports a confidence interval — carry it, do not report a bare number |
| Someone reads panel CNV as genome-wide | Say "focal, within target" on the report, not "copy number profile" |

---

## 8. Definition of done

1. PoN builds from ≥10 assay-matched normals and is checked into `~/data`.
2. `--cnv-pon` runs step 9 end to end on a real tumour and writes segments.
3. Purity and ploidy come out of PureCN with a confidence interval.
4. `<sample>.cna.tsv` validates against PCGR and the report shows CNV tiers.
5. A known-positive control — a sample with an established amplification —
   is called correctly. **Without this, the step is not validated**, only
   working.
6. `install_pipeline.py` gains a `cnv-pon` check; `--check` reports it missing
   rather than silently skipping.

Step 5 is the one that matters. Everything before it proves the plumbing.
