#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
panel_profiles.py
=================
One place that knows what a given cancer NGS panel needs, so a run is
configured by NAMING THE ASSAY rather than by remembering nine flags and
getting each of them right.

WHY THIS EXISTS
---------------
  The pipeline was written against one assay, and every setting that assay
  needed is a separate command-line flag with a general-purpose default.
  That is correct and it is also unusable across a laboratory that runs
  several kits: the operator who has just been handed a Thermo Oncomine
  run and an Illumina TSO500 run has to know, from memory, that

    * the amplicon kit must NOT have its duplicates marked and the capture
      kit must,
    * the capture kit wants 100 bp of interval padding and the amplicon
      kit wants none,
    * BQSR has too few sites to fit on a 20 kb hotspot panel,
    * TMB divided by PCGR's default 34 Mb is wrong for both of them, by
      different factors,
    * and the QIAseq run in the next directory carries its UMI in the
      first 12 bases of read 2.

  Every one of those is silent when it is wrong. The run finishes, the
  report renders, and the numbers are confidently incorrect. A profile
  turns that body of knowledge into data: `--panel thermo-oncomine-cav3`
  sets all of it at once, prints exactly what it set, and records it in
  the run manifest so a later reader can see which assay the numbers
  belong to.

WHAT A PROFILE IS, AND IS NOT
-----------------------------
  A profile is SETTINGS, not CONTENT. It carries the numbers, the
  chemistry-driven step choices and the UMI layout. It does NOT carry the
  target BED, because:

    * vendor BEDs are licensed content that travels with the kit, not
      something this bundle may redistribute;
    * they are version- and genome-build-specific, and the wrong version
      is worse than none at all -- it looks like a target region and
      silently is not one;
    * a laboratory's spike-ins and custom content mean the BED in the
      freezer room is frequently not the one on the vendor's website.

  So every panel profile names the BED as something the operator must
  supply (`--panel-bed`), and the pipeline measures the footprint from
  that file rather than trusting any number written here.

PRECEDENCE -- WHAT BEATS WHAT
-----------------------------
      explicit command-line flag   >   panel profile   >   script default

  A profile only ever fills in a setting the command line did not name, so
  adding `--panel` to a working command can never change a value that
  command was already stating. The application is printed line by line,
  with the reason, because a setting that arrived from somewhere the
  operator cannot see is exactly the kind of thing that makes two runs
  disagree for reasons nobody can reconstruct.

BUILT-INS ARE A STARTING POINT, NOT AN AUTHORITY
------------------------------------------------
  The vendor profiles below encode published, nominal figures for ONE
  version of each kit. Kits are revised, laboratories spike in extra
  content, and a validated assay's limit of detection is a property of
  that laboratory's validation, not of the catalogue number. Treat a
  built-in as a sensible first draft: copy it, adjust it to your validated
  assay, and keep the copy.

      python comprehensive_variant_calling.py --describe-panel tso500
      python comprehensive_variant_calling.py --new-panel-template mine.json
      # edit mine.json, then:
      python comprehensive_variant_calling.py --panel-file mine.json \
             --panel our-lung-panel ...

  User profiles are also picked up automatically from, in order:

      $CANCER_PIPELINE_PANEL_DIR   (colon-separated list of directories)
      ~/.config/cancer_pipeline/panels/*.json
      <directory holding these scripts>/panels/*.json

  A user profile whose id matches a built-in REPLACES it, so a site can
  correct a built-in without editing this file.

THE ONE THING A PROFILE CANNOT FIX
----------------------------------
  This pipeline has no UMI consensus step (no fgbio, no duplex
  collapsing). fastp moves a UMI into the read name, which is enough for
  UMI-aware deduplication downstream but is NOT error-corrected consensus
  calling. The UMI profiles below therefore configure extraction and then
  say so: a kit sold for 0.1% ctDNA detection will not reach 0.1% here,
  and a profile that quietly set --min-allele-fraction 0.001 would be
  claiming a sensitivity this code does not have.
"""

import copy
import json
import os
import re
import sys

# Where user-written profiles are looked for, in precedence order (later
# wins, so a site directory overrides a personal one only if it is listed
# later -- see panel_search_path()).
PANEL_DIR_ENV = "CANCER_PIPELINE_PANEL_DIR"
USER_PANEL_DIR = os.path.expanduser("~/.config/cancer_pipeline/panels")

# TMB below roughly this footprint is not interpretable: the estimate is a
# count divided by a very small denominator, so one extra artefact moves it
# by several mutations per Mb. Panels smaller than this get
# pcgr_estimate_tmb off and a note saying why, rather than a number that
# looks like a result. The figure is the commonly cited floor for panel
# TMB; it is not a bright line, which is why it is stated here once instead
# of being implied by each profile.
TMB_MIN_PANEL_MB = 1.0

# Mutational signature fitting wants a few hundred SNVs before it means
# anything. A 2 Mb panel produces tens, so signatures are left off below
# this footprint for the same reason TMB is.
SIGNATURE_MIN_PANEL_MB = 20.0


# =============================================================================
# SECTION 1: THE SCHEMA
# =============================================================================
# Every key a profile's "settings" block may carry, mapped to the value
# type and the argparse dest it fills. The dest names are deliberately
# IDENTICAL to the pipeline's own, so applying a profile is a setattr loop
# rather than a translation table that can drift away from the flags.
#
# A key that is not in here is an ERROR when a user profile is loaded, not
# a silently ignored line. A typo in a hand-written JSON file that quietly
# does nothing is the exact failure this module exists to prevent.

PANEL_SETTING_TYPES = {
    # --- target region --------------------------------------------------
    # Both usually point at the same BED; they are separate because they
    # answer different questions (what to call over vs what to make a
    # coverage statement about) and because a laboratory sometimes reports
    # coverage over a smaller clinically-reportable subset.
    "intervals": "path",
    "coverage_bed": "path",
    "interval_padding": "int",

    # --- call filtering -------------------------------------------------
    "min_allele_fraction": "float",
    "min_depth": "int",
    "foldback_min_match": "float",

    # --- coverage statement ---------------------------------------------
    "coverage_min_depth": "int",
    "coverage_min_mapq": "int",
    "coverage_min_baseq": "int",

    # --- QC / trimming (stage 1, also the pipeline's own fastp step) -----
    "adapter_r1": "str",
    "adapter_r2": "str",
    "min_read_length": "int",
    "cut_right": "bool",
    "correction": "bool",
    "umi_loc": "str",
    "umi_len": "int",
    "umi_skip": "int",

    # --- which steps run at all -----------------------------------------
    # skip_steps is a LIST and is merged with whatever the command line
    # asked for, never replaced: a profile that skips dedup plus a command
    # line that skips qc must end up skipping both.
    "skip_steps": "list",
    "skip_bqsr": "bool",
    # Deliberately NOT here: two-pass alignment. It is a throughput and
    # sensitivity trade-off for a particular machine and a particular
    # hurry, not a property of a kit, and it is the one setting the web
    # form exposes as a checkbox that a profile could also set -- an
    # unticked box is indistinguishable from an absent one, so a profile
    # that set it would silently re-enable what the operator had just
    # turned off.

    # --- caller passthrough ---------------------------------------------
    "mutect2_extra_args": "str",

    # --- RNA: alignment and fusion detection ----------------------------
    # RNA runs share this vocabulary rather than having a registry of their
    # own, because a laboratory's TSO500 DNA and TSO500 RNA libraries come
    # off the same specimen and belong in the same list. Which keys a
    # profile may usefully set is decided by its chemistry, not by a
    # separate file.
    "gtf": "path",
    "star_index": "path",
    "read_length": "int",
    "strandedness": "str",
    "fusion_caller": "str",
    "fusion_min_confidence": "str",
    "min_fusion_reads": "int",
    "rna_min_reads_millions": "float",
    "rna_min_unique_mapped_pct": "float",
    "star_extra_args": "str",
    "arriba_extra_args": "str",

    # --- clinical report -------------------------------------------------
    "pcgr_assay": "str",
    "pcgr_target_size_mb": "float",
    "pcgr_estimate_tmb": "bool",
    "pcgr_estimate_msi": "bool",
    "pcgr_estimate_signatures": "bool",
}

# The subset each script can actually act on. A profile is written once and
# handed to whichever stage is running; fastq_qc_clean.py must apply the
# trimming keys and ignore the caller keys rather than crashing on them.
QC_SETTINGS = frozenset({
    "adapter_r1", "adapter_r2", "min_read_length", "cut_right",
    "correction", "umi_loc", "umi_len", "umi_skip",
})

# The RNA branch: align_rna.py and fusion_calling.py apply these and
# ignore the somatic-caller keys, exactly as fastq_qc_clean.py ignores
# everything but the trimming ones.
RNA_SETTINGS = frozenset({
    "gtf", "star_index", "read_length", "strandedness", "fusion_caller",
    "fusion_min_confidence", "min_fusion_reads", "rna_min_reads_millions",
    "rna_min_unique_mapped_pct", "star_extra_args", "arriba_extra_args",
    # Shared with the DNA branch: trimming happens the same way, and a
    # degraded FFPE RNA library needs the read-length floor lowered.
    "min_read_length", "adapter_r1", "adapter_r2",
    "umi_loc", "umi_len", "umi_skip",
})

# coverage_report.py drops the "coverage_" prefix, because in a tool that
# does nothing but measure coverage there is nothing to distinguish it
# from.
COVERAGE_SETTINGS = frozenset({
    "coverage_min_depth", "coverage_min_mapq", "coverage_min_baseq",
})
COVERAGE_RENAME = {
    "coverage_min_depth": "min_depth",
    "coverage_min_mapq": "min_mapq",
    "coverage_min_baseq": "min_baseq",
}

# pcgr_report.py drops the "pcgr_" prefix for the same reason, and spells
# the TMB denominator the way PCGR itself does.
PCGR_SETTINGS = frozenset({
    "pcgr_assay", "pcgr_target_size_mb", "pcgr_estimate_tmb",
    "pcgr_estimate_msi", "pcgr_estimate_signatures",
})
PCGR_RENAME = {
    "pcgr_assay": "assay",
    "pcgr_target_size_mb": "effective_target_size_mb",
    "pcgr_estimate_tmb": "estimate_tmb",
    "pcgr_estimate_msi": "estimate_msi",
    "pcgr_estimate_signatures": "estimate_signatures",
}

# Chemistry vocabulary. Not free text: the chemistry is what decides
# whether duplicate marking is meaningful, and a typo'd value would make
# that decision silently wrong.
#
# The rna-* values do more than label. They route: a profile whose
# chemistry begins "rna-" describes a TRANSCRIPTOME library, which cannot
# be aligned by bwa-mem2 at all (reads cross exon-exon junctions), must not
# have its duplicates marked (a highly expressed gene legitimately produces
# thousands of identical reads), and has no somatic caller downstream. The
# web form and the orchestrator both read assay_type() rather than asking
# the operator to remember which scripts apply.
CHEMISTRIES = (
    "hybrid-capture", "amplicon", "wgs",
    "rna-capture",     # probe capture of a transcriptome subset (TSO500 RNA)
    "rna-amplicon",    # anchored multiplex PCR (Archer, Oncomine)
    "rna-total",       # total RNA / whole transcriptome, rRNA-depleted
)

# Which of the two pipelines a chemistry belongs to.
DNA_CHEMISTRIES = frozenset({"hybrid-capture", "amplicon", "wgs"})
RNA_CHEMISTRIES = frozenset({"rna-capture", "rna-amplicon", "rna-total"})

REQUIRED_PROFILE_KEYS = ("id", "name", "manufacturer", "chemistry")

# Strandedness vocabulary, matching the way the library was built. Wrong
# here is not fatal for fusion calling -- Arriba does not need it -- but it
# silently inverts any expression or QC statement that does.
STRANDEDNESS = ("unstranded", "forward", "reverse")


def partner_id(profile):
    """
    The other half of a hybrid kit, or None.

    Many cancer panels ship as a DNA library and an RNA library of the same
    specimen -- TSO500, Oncomine Comprehensive, Archer -- reported together
    as one result. Recording the pairing in the profile means the hybrid
    form can offer the matching half automatically instead of relying on an
    operator to know that 'illumina-tso500' and 'illumina-tso500-rna' are
    two halves of one assay.
    """
    return (profile or {}).get("pairs_with") or None


def assay_type(profile):
    """
    'dna' or 'rna' for a profile, derived from its chemistry.

    Derived rather than declared, so the two can never disagree. A profile
    that said assay_type 'dna' and chemistry 'rna-capture' would be a
    profile nobody could debug.
    """
    chemistry = (profile or {}).get("chemistry", "")
    return "rna" if chemistry in RNA_CHEMISTRIES else "dna"


# =============================================================================
# SECTION 2: THE BUILT-IN PROFILES
# =============================================================================
# Ordered most-general first, so --list-panels reads as "here is the shape
# of the thing, and here are the named kits that fit it".
#
# CONVENTIONS USED BELOW
#   nominal_target_size_mb  the vendor's published footprint for one kit
#                           version, used ONLY as a cross-check against the
#                           BED the operator supplies. The BED always wins;
#                           a disagreement is reported, because it usually
#                           means the wrong BED version or the wrong genome
#                           build.
#   notes                   printed by --describe-panel and carried into
#                           the run manifest. These are the caveats that
#                           would otherwise have to live in somebody's
#                           head.
#   requires                what the operator must still supply. A profile
#                           that cannot know the BED says so here instead
#                           of pretending the run is fully configured.

def _capture_defaults():
    """
    Settings shared by every hybrid-capture profile.

    Hybrid capture shears DNA randomly and then pulls target fragments out
    with probes. Two consequences drive everything here: fragment ends are
    random, so duplicate marking is meaningful and must run; and reads run
    past the probe edges, so a variant at the first or last base of a
    target is only callable if the flanking bases are in the interval.
    """
    return {
        # Reads extend beyond the probes, and a variant sitting on a target
        # boundary needs those flanking reads to be callable at all. 100 bp
        # is the usual figure for a 150 bp read length.
        "interval_padding": 100,
        # FFPE tissue somatic panels: 5% is the conventional floor and is
        # roughly where PCR-slippage indels and residual deamination stop
        # dominating. A laboratory with a validated LoD should use that
        # number instead -- this one is a safe default, not a validation.
        "min_allele_fraction": 0.05,
        # GATK has no depth filter of its own and will PASS a call made on
        # two reads. Panels are sequenced deep enough that anything under
        # 50x is a coverage dropout rather than a finding.
        "min_depth": 50,
        "pcgr_assay": "TARGETED",
    }


def _amplicon_defaults():
    """
    Settings shared by every amplicon profile.

    THE ONE THAT MATTERS: duplicate marking is DISABLED.

    Amplicon reads begin and end at primer coordinates by construction, so
    every read from the same amplicon looks like a duplicate of every other
    one. MarkDuplicates flags nearly the whole library, and because Mutect2
    ignores duplicate-flagged reads the call set is then built from a
    handful of surviving reads per amplicon -- a catastrophic, silent loss
    of depth that shows up only as an implausibly thin VCF. This is why the
    step is skipped rather than merely discouraged.

    Padding is 0 for the mirror-image reason: the target IS the amplicon,
    and padding walks the interval straight into primer sequence, where the
    bases come from the primer oligo rather than from the patient.
    """
    return {
        "skip_steps": ["dedup"],
        "interval_padding": 0,
        # Multiplex PCR generates its own low-level chimeras and slippage
        # products; the same 5% floor applies, for stronger reasons.
        "min_allele_fraction": 0.05,
        "min_depth": 50,
        "pcgr_assay": "TARGETED",
    }


def _small_panel_report_defaults():
    """
    Report settings for a panel too small for genome-wide statistics.

    TMB and signature fitting are arithmetic over a denominator. Below
    about 1 Mb the denominator is small enough that a couple of artefacts
    move TMB by several mutations per Mb, and signature fitting has tens of
    SNVs where it wants hundreds. Both are therefore OFF, and the profile's
    notes say so -- otherwise their absence from the report reads as
    "nothing found" rather than "never calculated".
    """
    return {
        "pcgr_estimate_tmb": False,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    }


def _rna_defaults():
    """
    Settings shared by every RNA profile.

    WHAT IS DELIBERATELY ABSENT HERE is the interesting part. There is no
    interval padding, no allele-fraction floor, no depth floor and no
    duplicate-marking decision, because none of those concepts survives the
    move from DNA to RNA:

      * A target BED does not restrict an RNA run. The transcriptome
        annotation (the GTF) does, and it is a different kind of object --
        it carries exon structure, which is the whole point.
      * Duplicate marking is not merely unhelpful, it is wrong. Two reads
        from the same position in a highly expressed gene are two
        observations of an abundant transcript, not one molecule counted
        twice. The RNA pipeline has no dedup step to switch off.
      * There is no allele fraction, because there is no somatic caller.
        This branch detects FUSIONS. Confidence comes from the number and
        the quality of reads crossing a breakpoint.

    What replaces them is a read-count floor and a mapping-rate floor: on
    RNA, "was there enough usable library?" is the question a negative
    result depends on, the way target coverage is on DNA.
    """
    return {
        "fusion_caller": "arriba",
        # Arriba grades every call low/medium/high from breakpoint support,
        # the read-through and blacklist checks, and whether the fusion is
        # in frame. "medium" keeps the clinically plausible calls without
        # the long tail that low admits; the report shows all three and
        # marks which cleared this bar.
        "fusion_min_confidence": "medium",
        "min_fusion_reads": 3,
        "strandedness": "unstranded",
        # FFPE RNA is fragmented, and a 50 bp floor inherited from the DNA
        # side discards a large fraction of a degraded library. 35 keeps
        # reads that still span a junction.
        "min_read_length": 35,
        "rna_min_reads_millions": 10.0,
        "rna_min_unique_mapped_pct": 60.0,
    }


def _builtin_profiles():
    """
    Build the built-in registry.

    A function rather than a module-level literal so the shared defaults
    above are composed once and cannot drift between profiles that are
    meant to agree.
    """
    profiles = []

    # ---------------------------------------------------------------
    # GENERIC PROFILES -- correct for any kit of that shape, and the
    # right starting point for a manufacturer not listed below.
    # ---------------------------------------------------------------
    generic_capture = dict(_capture_defaults())
    generic_capture.update({
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    profiles.append({
        "id": "generic-capture",
        "name": "Generic hybrid-capture panel",
        "manufacturer": "any",
        "chemistry": "hybrid-capture",
        "aliases": ["capture", "hybrid-capture"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Probe-capture solid-tumour panel of any size. Use this "
                   "when the kit is not listed and is not amplicon-based.",
        "requires": ["panel BED"],
        "notes": [
            "TMB is enabled but is only interpretable above roughly 1 Mb of "
            "target; the run reports the footprint it measured from your "
            "BED, so check it.",
            "Set --min-allele-fraction from your assay's validated limit of "
            "detection. The 0.05 here is a conventional floor, not a "
            "validation result.",
        ],
        "settings": generic_capture,
    })

    generic_capture_umi = dict(generic_capture)
    generic_capture_umi.update({
        # Overlap correction rewrites bases where the two mates disagree,
        # which is precisely the evidence a UMI workflow needs to keep.
        "correction": False,
    })
    profiles.append({
        "id": "generic-capture-umi",
        "name": "Generic hybrid-capture panel with UMIs",
        "manufacturer": "any",
        "chemistry": "hybrid-capture",
        "aliases": ["capture-umi"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "As generic-capture, for a kit whose reads still carry "
                   "their molecular barcodes. Supply --umi-loc/--umi-len "
                   "for your kit's read structure.",
        "requires": ["panel BED", "--umi-loc and --umi-len for your kit"],
        "notes": [
            "This pipeline EXTRACTS UMIs (into the read name) but does not "
            "build consensus reads: there is no fgbio/duplex step. The "
            "barcode is preserved for auditing and for UMI-aware "
            "deduplication, and that is all -- do not expect the error "
            "suppression the kit is sold for.",
            "UMI read structure differs between kit versions. Read it off "
            "your kit handbook rather than assuming this one.",
        ],
        "settings": generic_capture_umi,
    })

    generic_amplicon = dict(_amplicon_defaults())
    generic_amplicon.update(_small_panel_report_defaults())
    profiles.append({
        "id": "generic-amplicon",
        "name": "Generic amplicon / multiplex-PCR panel",
        "manufacturer": "any",
        "chemistry": "amplicon",
        "aliases": ["amplicon"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Multiplex-PCR panel of any size: duplicate marking off, "
                   "no interval padding, report statistics that need a "
                   "large footprint switched off.",
        "requires": ["panel BED"],
        "notes": [
            "Duplicate marking is SKIPPED. Amplicon reads share primer "
            "coordinates, so MarkDuplicates would flag nearly the whole "
            "library and Mutect2 would then call on what little survived.",
            "Primer bases must already be trimmed by the instrument or "
            "vendor software. This pipeline does not trim them, and an "
            "untrimmed primer reports the oligo's sequence as the "
            "patient's.",
            "BQSR is left ON here because panel size varies; on a panel "
            "under about 1 Mb add --skip-bqsr (too few covered sites to fit "
            "a recalibration model).",
        ],
        "settings": generic_amplicon,
    })

    hotspot = dict(_amplicon_defaults())
    hotspot.update(_small_panel_report_defaults())
    hotspot.update({
        # A hotspot panel covers tens of kilobases. BQSR fits an empirical
        # error model from covered sites, and there are nowhere near enough
        # of them here; the model it would produce is noise.
        "skip_bqsr": True,
        # Hotspot panels are sequenced very deep, and the whole point is
        # sensitivity at known positions.
        "min_depth": 100,
    })
    profiles.append({
        "id": "generic-hotspot",
        "name": "Generic hotspot panel (< ~100 kb)",
        "manufacturer": "any",
        "chemistry": "amplicon",
        "aliases": ["hotspot"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Small amplicon hotspot panel: as generic-amplicon, plus "
                   "BQSR off and a deeper depth floor.",
        "requires": ["panel BED"],
        "notes": [
            "BQSR is skipped: a few tens of kilobases do not contain enough "
            "covered sites to fit a recalibration model, and fitting one on "
            "too little data distorts the base qualities it is meant to "
            "correct.",
            "A negative on a hotspot panel means only 'not at the positions "
            "this kit interrogates'. The coverage report is the only thing "
            "that says which those were -- pass --panel-bed.",
        ],
        "settings": hotspot,
    })

    ctdna = dict(_capture_defaults())
    ctdna.update({
        # Plasma panels are sequenced to several thousand x, and the whole
        # design intent is low-fraction detection -- so the depth floor
        # rises and the VAF floor falls. It does NOT fall to the 0.1% the
        # kits advertise; see the note.
        "min_depth": 500,
        "min_allele_fraction": 0.005,
        "correction": False,
        "pcgr_estimate_tmb": False,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    profiles.append({
        "id": "ctdna-capture",
        "name": "Circulating tumour DNA capture panel (plasma)",
        "manufacturer": "any",
        "chemistry": "hybrid-capture",
        "aliases": ["ctdna", "liquid-biopsy"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "Deep plasma panel: high depth floor, low VAF floor, "
                   "report statistics that assume tissue switched off.",
        "requires": ["panel BED", "--umi-loc/--umi-len if UMIs are in-read"],
        "notes": [
            "READ THIS BEFORE USING IT CLINICALLY. Commercial ctDNA assays "
            "reach 0.1-0.5% VAF by building consensus reads from UMI "
            "families. This pipeline has no consensus step, so its floor is "
            "set by raw polymerase and sequencing error, around 0.5-1%. The "
            "0.005 here is the lowest defensible setting for that, not the "
            "kit's specification.",
            "Clonal haematopoiesis is the dominant false-positive source in "
            "plasma and is indistinguishable from tumour signal without a "
            "matched buffy-coat normal. Sequence one.",
            "TMB and signatures are off: neither is meaningful from a "
            "low-fraction plasma call set on a small footprint.",
        ],
        "settings": ctdna,
    })

    wes = {
        "interval_padding": 100,
        "min_allele_fraction": 0.05,
        "min_depth": 20,
        "pcgr_assay": "WES",
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": True,
        "pcgr_estimate_msi": False,
    }
    profiles.append({
        "id": "wes",
        "name": "Whole exome",
        "manufacturer": "any",
        "chemistry": "hybrid-capture",
        "aliases": ["exome"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": 34.0,
        "summary": "Exome capture of any vendor. TMB and signature fitting "
                   "are enabled -- an exome is large enough for both.",
        "requires": ["exome BED"],
        "notes": [
            "MSI stays off: PCGR restricts its own MSI prediction to "
            "tumour-CONTROL WGS/WES runs, so on a tumour-only exome the "
            "flag produces nothing. Supply --msi-models for MSIsensor2 "
            "instead.",
            "Depth floor of 20 suits a 100x exome; raise it if yours is "
            "deeper.",
        ],
        "settings": wes,
    })

    wgs = {
        "interval_padding": 0,
        "min_allele_fraction": 0.05,
        "min_depth": 10,
        "pcgr_assay": "WGS",
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": True,
        "pcgr_estimate_msi": False,
    }
    profiles.append({
        "id": "wgs",
        "name": "Whole genome",
        "manufacturer": "any",
        "chemistry": "wgs",
        "aliases": ["genome"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "No target restriction. PCGR computes its own WGS "
                   "denominator, so no panel size is set.",
        "requires": [],
        "notes": [
            "No --intervals, so Mutect2 walks the whole genome: this is the "
            "one configuration where that is correct rather than a mistake.",
            "Expect a long run. Nothing here is tuned for throughput.",
        ],
        "settings": wgs,
    })

    # ---------------------------------------------------------------
    # NAMED VENDOR KITS
    #
    # Figures are the vendor's published, nominal values for ONE version
    # of each kit, carried here as a cross-check against the BED you
    # supply -- never as a substitute for it. Confirm against your kit
    # documentation; where a value could not be stated with confidence it
    # is left null rather than guessed, because a wrong number here would
    # be silently believed.
    # ---------------------------------------------------------------

    tso500 = dict(_capture_defaults())
    tso500.update({
        "correction": False,
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    profiles.append({
        "id": "illumina-tso500",
        "pairs_with": "illumina-tso500-rna",
        "name": "TruSight Oncology 500 (DNA workflow)",
        "manufacturer": "Illumina",
        "chemistry": "hybrid-capture",
        "aliases": ["tso500"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": 1.94,
        "summary": "523-gene pan-cancer capture panel, nominally 1.94 Mb of "
                   "DNA target, UMI-tagged library.",
        "requires": ["panel BED (ships with the kit / DRAGEN TSO500 app)"],
        "notes": [
            "TSO500 libraries carry UMIs. If your demultiplexing leaves "
            "them in the reads, add --umi-loc/--umi-len for the structure "
            "your kit version uses; if the vendor pipeline already "
            "consumed them, leave those unset.",
            "1.94 Mb is close enough to the 1 Mb floor that TMB carries a "
            "wide confidence interval, and tumour-only TMB is unreliable "
            "besides. Treat it as internal QC unless a matched normal was "
            "sequenced.",
            "This profile is the DNA half. The RNA workflow -- fusions "
            "and splice variants -- is --panel illumina-tso500-rna, which "
            "runs the RNA pipeline instead of the somatic caller.",
        ],
        "settings": tso500,
    })

    tst170 = dict(_capture_defaults())
    tst170.update(_small_panel_report_defaults())
    profiles.append({
        "id": "illumina-tst170",
        "name": "TruSight Tumor 170 (DNA component)",
        "manufacturer": "Illumina",
        "chemistry": "hybrid-capture",
        "aliases": ["tst170"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "170-gene capture panel. Footprint is measured from your "
                   "BED; no nominal figure is asserted here.",
        "requires": ["panel BED"],
        "notes": [
            "Report statistics needing a large footprint are off. If your "
            "measured footprint turns out to exceed 1 Mb and you want TMB, "
            "pass --pcgr-estimate-tmb explicitly.",
            "DNA component only. For the RNA fusion workflow use an RNA "
            "profile (generic-rna-capture, or the kit's own).",
        ],
        "settings": tst170,
    })

    ampliseq_illumina = dict(_amplicon_defaults())
    ampliseq_illumina.update(_small_panel_report_defaults())
    ampliseq_illumina.update({"skip_bqsr": True, "min_depth": 100})
    profiles.append({
        "id": "illumina-ampliseq-cancer-hotspot",
        "name": "AmpliSeq for Illumina Cancer HotSpot Panel v2",
        "manufacturer": "Illumina",
        "chemistry": "amplicon",
        "aliases": ["ampliseq-illumina"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "50-gene amplicon hotspot panel: duplicates off, BQSR "
                   "off, deep depth floor.",
        "requires": ["panel BED"],
        "notes": [
            "Duplicate marking and BQSR are both off -- see generic-hotspot "
            "for why each would be actively harmful here.",
            "Hotspot coverage means a negative speaks only for the listed "
            "positions. Pass --panel-bed so the run says which those were.",
        ],
        "settings": ampliseq_illumina,
    })

    oncomine_cav3 = dict(_amplicon_defaults())
    oncomine_cav3.update({
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
        "min_depth": 100,
    })
    profiles.append({
        "id": "thermo-oncomine-cav3",
        "pairs_with": "thermo-oncomine-rna-fusion",
        "name": "Oncomine Comprehensive Assay v3",
        "manufacturer": "Thermo Fisher",
        "chemistry": "amplicon",
        "aliases": ["oncomine-cav3", "ocav3"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "161-gene amplicon panel. Duplicate marking off; "
                   "footprint measured from your BED.",
        "requires": ["panel BED (exported from Ion Reporter / AmpliSeq "
                     "Designer)"],
        "notes": [
            "Ion Torrent reads are single-end and have a different indel "
            "error profile from Illumina reads, particularly in "
            "homopolymers. This pipeline's QC, alignment and filtering "
            "defaults were tuned on Illumina paired-end data -- validate "
            "before trusting indel calls from an Ion run.",
            "TMB is on because the panel may exceed 1 Mb; check the "
            "measured footprint the run prints, and turn it off if it does "
            "not.",
            "Vendor BEDs are frequently in Ensembl contig naming (1, MT). "
            "The run detects that and writes a translated copy rather than "
            "failing at the caller.",
        ],
        "settings": oncomine_cav3,
    })

    ion_chpv2 = dict(_amplicon_defaults())
    ion_chpv2.update(_small_panel_report_defaults())
    ion_chpv2.update({"skip_bqsr": True, "min_depth": 100})
    profiles.append({
        "id": "thermo-ion-ampliseq-chpv2",
        "name": "Ion AmpliSeq Cancer Hotspot Panel v2",
        "manufacturer": "Thermo Fisher",
        "chemistry": "amplicon",
        "aliases": ["chpv2", "ion-chpv2"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "50-gene hotspot amplicon panel: duplicates off, BQSR "
                   "off, report statistics off.",
        "requires": ["panel BED"],
        "notes": [
            "See thermo-oncomine-cav3 on Ion Torrent read characteristics; "
            "the same caveat applies and matters more on a hotspot panel, "
            "where a handful of positions carry the entire result.",
        ],
        "settings": ion_chpv2,
    })

    focus = dict(_amplicon_defaults())
    focus.update(_small_panel_report_defaults())
    focus.update({"skip_bqsr": True, "min_depth": 100})
    profiles.append({
        "id": "thermo-oncomine-focus",
        "name": "Oncomine Focus Assay (DNA)",
        "manufacturer": "Thermo Fisher",
        "chemistry": "amplicon",
        "aliases": ["oncomine-focus"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Small solid-tumour amplicon panel, DNA component only.",
        "requires": ["panel BED"],
        "notes": [
            "DNA component only. Fusion detection is the RNA pipeline: "
            "see --panel thermo-oncomine-rna-fusion.",
            "See thermo-oncomine-cav3 for the Ion Torrent caveat.",
        ],
        "settings": focus,
    })

    sureselect_hs2 = dict(_capture_defaults())
    sureselect_hs2.update({
        "correction": False,
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    profiles.append({
        "id": "agilent-sureselect-xths2",
        "name": "SureSelect XT HS2 DNA (custom or catalogue)",
        "manufacturer": "Agilent",
        "chemistry": "hybrid-capture",
        "aliases": ["sureselect-hs2", "xths2"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Capture library with in-line molecular barcodes. "
                   "Footprint depends on the design, so it is measured "
                   "from your BED.",
        "requires": ["design BED from SureDesign",
                     "--umi-loc/--umi-len for your kit version"],
        "notes": [
            "The molecular barcode length and position depend on the kit "
            "version and on whether AGeNT already trimmed them. Set "
            "--umi-loc/--umi-len from your handbook; leave them unset if "
            "AGeNT consumed the barcodes upstream.",
            "Extraction is not consensus calling -- see generic-capture-umi.",
        ],
        "settings": sureselect_hs2,
    })

    sureselect_exome = dict(wes)
    profiles.append({
        "id": "agilent-sureselect-exome",
        "name": "SureSelect Human All Exon",
        "manufacturer": "Agilent",
        "chemistry": "hybrid-capture",
        "aliases": ["sureselect-exome"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": 35.0,
        "summary": "Exome capture. Nominal footprint varies by version "
                   "(V7/V8 are around 35 Mb); yours is measured from the "
                   "BED.",
        "requires": ["exome BED for your kit version"],
        "notes": [
            "Use the covered-regions BED for your exact version. Mixing a "
            "V7 BED with a V8 library reports coverage for regions the "
            "library never targeted.",
        ],
        "settings": sureselect_exome,
    })

    twist_capture = dict(_capture_defaults())
    twist_capture.update({
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    profiles.append({
        "id": "twist-custom-capture",
        "name": "Twist custom / catalogue capture panel",
        "manufacturer": "Twist Bioscience",
        "chemistry": "hybrid-capture",
        "aliases": ["twist"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Any Twist hybridisation panel. Footprint measured from "
                   "the design BED.",
        "requires": ["design BED"],
        "notes": [
            "Twist designs are commonly distributed in hg38 with UCSC "
            "contig naming, which matches this pipeline's reference "
            "without translation. The run checks and says so either way.",
        ],
        "settings": twist_capture,
    })

    twist_exome = dict(wes)
    profiles.append({
        "id": "twist-exome",
        "name": "Twist Exome",
        "manufacturer": "Twist Bioscience",
        "chemistry": "hybrid-capture",
        "aliases": [],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": 36.0,
        "summary": "Twist exome capture; nominal footprint varies by "
                   "version (Exome 2.0 and Comprehensive Exome differ).",
        "requires": ["exome BED for your kit version"],
        "notes": [
            "Check which Twist exome you have. The versions differ by "
            "several Mb, which is a several-percent error in TMB if the "
            "wrong BED is used.",
        ],
        "settings": twist_exome,
    })

    idt_pan_cancer = dict(_capture_defaults())
    idt_pan_cancer.update(_small_panel_report_defaults())
    profiles.append({
        "id": "idt-xgen-pan-cancer",
        "name": "xGen Pan-Cancer Hyb Panel",
        "manufacturer": "IDT",
        "chemistry": "hybrid-capture",
        "aliases": ["xgen-pan-cancer"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Focused pan-cancer capture panel, well under 1 Mb in "
                   "its catalogue form -- report statistics off.",
        "requires": ["panel BED"],
        "notes": [
            "TMB is off: the catalogue design is too small a denominator "
            "for it to mean anything. If you have spiked the design up past "
            "1 Mb, pass --pcgr-estimate-tmb explicitly.",
        ],
        "settings": idt_pan_cancer,
    })

    idt_exome = dict(wes)
    profiles.append({
        "id": "idt-xgen-exome",
        "name": "xGen Exome Hyb Panel",
        "manufacturer": "IDT",
        "chemistry": "hybrid-capture",
        "aliases": ["xgen-exome"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": 34.0,
        "summary": "Exome capture, nominally around 34 Mb of target.",
        "requires": ["exome BED for your kit version"],
        "notes": [
            "Use the version-matched BED; v1 and v2 differ.",
        ],
        "settings": idt_exome,
    })

    qiaseq = dict(_amplicon_defaults())
    qiaseq.update(_small_panel_report_defaults())
    qiaseq.update({
        # QIAseq Targeted DNA panels use single-primer extension: one end
        # of each fragment is a primer, the other is a random shear point,
        # and a 12-base UMI sits at the 5' end of read 2 ahead of an
        # 11-base common sequence. Verify against your handbook -- the
        # layout has varied between product lines.
        "umi_loc": "read2",
        "umi_len": 12,
        "umi_skip": 11,
        "correction": False,
    })
    profiles.append({
        "id": "qiagen-qiaseq-dna",
        "name": "QIAseq Targeted DNA Panel",
        "manufacturer": "QIAGEN",
        "chemistry": "amplicon",
        "aliases": ["qiaseq"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Single-primer-extension panel with in-read UMIs. "
                   "Duplicates off, UMI extracted from read 2.",
        "requires": ["panel BED", "confirmation of the UMI layout"],
        "notes": [
            "The UMI layout set here (12 bases at the start of read 2, then "
            "an 11-base common sequence) is the published QIAseq Targeted "
            "DNA structure. CONFIRM IT against your kit handbook before "
            "relying on it: a wrong offset silently truncates real sequence "
            "off the front of every read 2.",
            "QIAseq's smCounter workflow derives its sensitivity from UMI "
            "consensus. This pipeline extracts the UMI but does not build "
            "consensus reads, so it will not reach the kit's stated limit "
            "of detection.",
            "Single-primer extension leaves one random end per fragment, so "
            "duplicate marking is less catastrophic here than on a "
            "two-primer amplicon -- but it is still primer-anchored at one "
            "end, so it stays off.",
        ],
        "settings": qiaseq,
    })

    archer = dict(_amplicon_defaults())
    archer.update(_small_panel_report_defaults())
    archer.update({"correction": False})
    profiles.append({
        "id": "archer-variantplus",
        "pairs_with": "archer-fusionplex",
        "name": "Archer VariantPlus (AMP chemistry)",
        "manufacturer": "Invitae / ArcherDX",
        "chemistry": "amplicon",
        "aliases": ["archer"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Anchored multiplex PCR panel with molecular barcodes: "
                   "one gene-specific primer, one universal adapter.",
        "requires": ["panel BED", "--umi-loc/--umi-len for your kit"],
        "notes": [
            "Anchored multiplex PCR is primer-anchored at one end only, "
            "like SPE. Duplicate marking still stays off.",
            "Archer's own analysis software performs the barcode "
            "deduplication this pipeline does not. Expect different numbers "
            "from the two, and validate before substituting one for the "
            "other.",
        ],
        "settings": archer,
    })

    roche = dict(_capture_defaults())
    roche.update({
        "pcgr_estimate_tmb": True,
        "pcgr_estimate_signatures": False,
        "pcgr_estimate_msi": False,
    })
    # ---------------------------------------------------------------
    # RNA PROFILES
    #
    # A separate pipeline downstream (align_rna.py, fusion_calling.py),
    # but the same registry, because the DNA and RNA libraries of one kit
    # come off one specimen and an operator should find them in one list.
    # ---------------------------------------------------------------

    rna_capture = dict(_rna_defaults())
    profiles.append({
        "id": "generic-rna-capture",
        "name": "Generic RNA capture panel",
        "manufacturer": "any",
        "chemistry": "rna-capture",
        "aliases": ["rna-capture"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "Probe-capture RNA panel for fusion detection. Use this "
                   "when the kit is not listed and is not anchored-PCR.",
        "requires": ["GENCODE GTF and a STAR index built from it",
                     "the library's read length, for the index"],
        "notes": [
            "Fusion calling only. This branch does not quantify expression "
            "and does not call variants from RNA -- see RNA_SCOPE.md for "
            "what was deliberately left out and why.",
            "A negative fusion result is only as good as the library. Check "
            "the RNA QC report before reporting one: a degraded FFPE "
            "extraction produces a clean, empty fusion table that looks "
            "exactly like a true negative.",
        ],
        "settings": rna_capture,
    })

    rna_amplicon = dict(_rna_defaults())
    rna_amplicon.update({
        # Anchored multiplex PCR gives one gene-specific primer and one
        # universal adapter, so a fusion is detected from the partner side
        # without needing both partners in the design. It also means
        # coverage is wildly uneven by construction, and a mapping-rate
        # floor tuned for capture is the wrong test.
        "rna_min_reads_millions": 2.0,
        "rna_min_unique_mapped_pct": 40.0,
    })
    profiles.append({
        "id": "generic-rna-amplicon",
        "name": "Generic anchored-PCR RNA fusion panel",
        "manufacturer": "any",
        "chemistry": "rna-amplicon",
        "aliases": ["rna-amplicon"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "Anchored multiplex PCR fusion panel: far fewer reads "
                   "needed, and much lower mapping rates are normal.",
        "requires": ["GENCODE GTF and a STAR index", "read length"],
        "notes": [
            "READ THIS BEFORE USING IT CLINICALLY. Arriba and STAR-Fusion "
            "are built for capture and whole-transcriptome libraries. On "
            "anchored-PCR data their assumptions about read distribution "
            "and about what a read-through artefact looks like do not "
            "hold, and vendors ship their own callers for exactly this "
            "reason. Treat results here as orthogonal evidence, not as a "
            "replacement for the vendor pipeline.",
            "Primers must already be trimmed. This pipeline does not trim "
            "them, and an untrimmed primer is reported as the specimen's "
            "own sequence.",
            "The read-count floor is much lower than a capture panel's "
            "because the design enriches so hard -- 2 million reads is a "
            "normal, adequate anchored-PCR library.",
        ],
        "settings": rna_amplicon,
    })

    rna_total = dict(_rna_defaults())
    rna_total.update({
        # Whole transcriptome: plenty of reads, and the full complement of
        # real junctions, so the confidence bar can afford to be higher and
        # the mapping-rate expectation is a proper one.
        "rna_min_reads_millions": 40.0,
        "rna_min_unique_mapped_pct": 75.0,
        "min_read_length": 50,
    })
    profiles.append({
        "id": "rna-total",
        "name": "Total RNA / whole transcriptome",
        "manufacturer": "any",
        "chemistry": "rna-total",
        "aliases": ["wts", "rnaseq", "whole-transcriptome"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "rRNA-depleted total RNA. The configuration the fusion "
                   "callers were actually developed on.",
        "requires": ["GENCODE GTF and a STAR index", "read length"],
        "notes": [
            "This is the library type Arriba and STAR-Fusion were "
            "developed and benchmarked on, so it is the one where their "
            "published performance is the relevant number.",
            "Expect a long run and a large STAR index. Nothing here is "
            "tuned for throughput.",
        ],
        "settings": rna_total,
    })

    tso500_rna = dict(_rna_defaults())
    tso500_rna.update({"rna_min_reads_millions": 20.0})
    profiles.append({
        "id": "illumina-tso500-rna",
        "pairs_with": "illumina-tso500",
        "name": "TruSight Oncology 500 (RNA workflow)",
        "manufacturer": "Illumina",
        "chemistry": "rna-capture",
        "aliases": ["tso500-rna"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "The RNA half of TSO500: 55 genes for fusions and splice "
                   "variants, capture chemistry.",
        "requires": ["GENCODE GTF and a STAR index", "read length"],
        "notes": [
            "Pair this with --panel illumina-tso500 on the DNA library from "
            "the same specimen. They are one assay reported together, and "
            "the web form records the pairing so the two runs can be read "
            "as one result.",
            "MET exon 14 skipping is a splice event within one gene rather "
            "than a fusion between two. The fusion report labels those "
            "separately -- they are easy to miss in a table sorted by "
            "gene pair.",
        ],
        "settings": tso500_rna,
    })

    archer = dict(rna_amplicon)
    profiles.append({
        "id": "archer-fusionplex",
        "pairs_with": "archer-variantplus",
        "name": "Archer FusionPlex (AMP chemistry)",
        "manufacturer": "Invitae / ArcherDX",
        "chemistry": "rna-amplicon",
        "aliases": ["fusionplex"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "Anchored multiplex PCR fusion panel with molecular "
                   "barcodes; detects fusions with an unknown partner.",
        "requires": ["GENCODE GTF and a STAR index", "read length"],
        "notes": [
            "Archer's own analysis software is the validated route for "
            "this chemistry, and it uses the molecular barcodes this "
            "pipeline only extracts. Expect different calls from the two.",
            "AMP detects a fusion from one known partner, so a novel "
            "partner gene is a real result here rather than an artefact -- "
            "do not filter on 'both partners in the panel'.",
            "See generic-rna-amplicon for why capture-trained callers are "
            "being asked to do something they were not built for.",
        ],
        "settings": archer,
    })

    oncomine_rna = dict(rna_amplicon)
    profiles.append({
        "id": "thermo-oncomine-rna-fusion",
        "pairs_with": "thermo-oncomine-cav3",
        "name": "Oncomine fusion assay (RNA)",
        "manufacturer": "Thermo Fisher",
        "chemistry": "rna-amplicon",
        "aliases": ["oncomine-rna"],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "Amplicon RNA fusion panel, the RNA half of the Oncomine "
                   "solid-tumour assays.",
        "requires": ["GENCODE GTF and a STAR index", "read length"],
        "notes": [
            "Ion Torrent reads are single-end. STAR handles that, but every "
            "paired-end assumption in this bundle's QC does not -- and the "
            "fusion callers' breakpoint support statistics are weaker "
            "without mate pairs.",
            "The vendor's own fusion caller is the validated route; see "
            "generic-rna-amplicon.",
        ],
        "settings": oncomine_rna,
    })

    profiles.append({
        "id": "roche-kapa-target",
        "name": "KAPA HyperCap / HyperExome target enrichment",
        "manufacturer": "Roche",
        "chemistry": "hybrid-capture",
        "aliases": ["kapa", "hypercap"],
        "genome_builds": ["hg38", "hg19"],
        "nominal_target_size_mb": None,
        "summary": "Roche capture panel or exome. Footprint measured from "
                   "the design files that ship with the kit.",
        "requires": ["design/primary-target BED"],
        "notes": [
            "Roche ships both a capture-target and a primary-target BED. "
            "Use the PRIMARY target as --panel-bed: the capture target "
            "includes probe regions you should not be reporting coverage "
            "over.",
            "For a HyperExome, prefer --panel wes and pass the exome BED.",
        ],
        "settings": roche,
    })

    return profiles


# =============================================================================
# SECTION 3: LOADING, VALIDATING AND RESOLVING PROFILES
# =============================================================================

class PanelError(Exception):
    """A profile could not be loaded, validated or resolved by name."""


def _validate_profile(profile, origin):
    """
    Check one profile, returning it, or raise PanelError naming the file.

    Strict on purpose. A hand-written profile is configuration that changes
    variant calls, and the failure mode of leniency here is a key that is
    quietly ignored -- a run configured differently from the file its
    operator is reading.
    """
    if not isinstance(profile, dict):
        raise PanelError(f"{origin}: expected a JSON object per profile, "
                         f"found {type(profile).__name__}")

    for key in REQUIRED_PROFILE_KEYS:
        if not profile.get(key):
            raise PanelError(f"{origin}: profile is missing required key "
                             f"'{key}'")

    if not re.match(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", str(profile["id"])):
        raise PanelError(f"{origin}: id {profile['id']!r} must be letters, "
                         f"digits, dot, dash or underscore -- it is typed on "
                         f"a command line")

    if profile["chemistry"] not in CHEMISTRIES:
        raise PanelError(
            f"{origin}: chemistry {profile['chemistry']!r} is not one of "
            f"{', '.join(CHEMISTRIES)}. This is not cosmetic: chemistry is "
            f"what decides whether duplicate marking is meaningful.")

    settings = profile.get("settings") or {}
    if not isinstance(settings, dict):
        raise PanelError(f"{origin}: 'settings' must be a JSON object")

    for key, value in settings.items():
        kind = PANEL_SETTING_TYPES.get(key)
        if kind is None:
            close = _closest(key, PANEL_SETTING_TYPES)
            hint = f" Did you mean '{close}'?" if close else ""
            raise PanelError(
                f"{origin}: unknown setting '{key}'.{hint} Settings are "
                f"named after the pipeline's own options with dashes "
                f"replaced by underscores; run --describe-panel on a "
                f"built-in to see the spelling.")
        if value is None:
            continue                      # explicit null = "leave default"
        if kind == "list":
            if not isinstance(value, list) or \
                    not all(isinstance(v, str) for v in value):
                raise PanelError(f"{origin}: '{key}' must be a list of "
                                 f"strings")
        elif kind == "bool":
            if not isinstance(value, bool):
                raise PanelError(f"{origin}: '{key}' must be true or false")
        elif kind == "int":
            if isinstance(value, bool) or not isinstance(value, int):
                raise PanelError(f"{origin}: '{key}' must be a whole number")
        elif kind == "float":
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise PanelError(f"{origin}: '{key}' must be a number")
        elif not isinstance(value, str):
            raise PanelError(f"{origin}: '{key}' must be a string")

    profile = copy.deepcopy(profile)
    partner = profile.get("pairs_with")
    if partner is not None and not isinstance(partner, str):
        raise PanelError(f"{origin}: 'pairs_with' must be the id of the "
                         f"other half of the kit, as a string")

    profile.setdefault("pairs_with", None)
    profile.setdefault("aliases", [])
    profile.setdefault("notes", [])
    profile.setdefault("requires", [])
    profile.setdefault("genome_builds", [])
    profile.setdefault("summary", "")
    profile.setdefault("nominal_target_size_mb", None)
    profile.setdefault("settings", {})
    profile["source"] = origin
    return profile


def _closest(word, candidates):
    """Nearest candidate by difflib, or None. Used only for error hints."""
    try:
        import difflib
        match = difflib.get_close_matches(word, list(candidates), n=1,
                                          cutoff=0.6)
        return match[0] if match else None
    except Exception:                     # noqa: BLE001 - hints are optional
        return None


def load_panel_file(path):
    """
    Read one JSON file holding a profile, a list of them, or {"panels": []}.

    All three shapes are accepted because all three are what people
    actually write: one kit per file, a laboratory's whole set in one file,
    or a file with a comment-ish wrapper key around the list.
    """
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except OSError as exc:
        raise PanelError(f"cannot read panel file {path}: {exc}")
    except ValueError as exc:
        raise PanelError(f"{path} is not valid JSON: {exc}")

    if isinstance(data, dict) and "panels" in data:
        data = data["panels"]
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        raise PanelError(f"{path}: expected a profile object or a list of "
                         f"them")
    return [_validate_profile(p, path) for p in data]


def panel_search_path(extra=None):
    """
    Directories scanned for user profiles, lowest precedence first.

    The script's own panels/ directory comes first so that a bundle can
    ship site profiles alongside the code; the user's config directory
    overrides it, and $CANCER_PIPELINE_PANEL_DIR overrides both, because
    the environment variable is the most deliberate of the three.
    """
    dirs = [os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "panels"),
            USER_PANEL_DIR]
    env = os.environ.get(PANEL_DIR_ENV, "")
    dirs += [d for d in env.split(os.pathsep) if d]
    dirs += list(extra or [])
    return dirs


def available_panels(extra_files=None, extra_dirs=None, warn=None):
    """
    The full registry: built-ins, then every user profile found, by id.

    A user profile with the id of a built-in REPLACES it. That is the
    supported way to correct a built-in for a local kit revision without
    editing this file and losing the change at the next update.

    A broken file in a search directory WARNS and is skipped, rather than
    stopping an otherwise valid run: an unrelated stray file in a config
    directory must not be able to block variant calling. A file named
    explicitly with --panel-file is different -- that one raises, because
    the operator asked for it by name.
    """
    registry = {}
    for profile in _builtin_profiles():
        _install(registry, _validate_profile(profile, "built-in"))

    for directory in panel_search_path(extra_dirs):
        if not os.path.isdir(directory):
            continue
        for name in sorted(os.listdir(directory)):
            if not name.endswith(".json"):
                continue
            path = os.path.join(directory, name)
            try:
                for profile in load_panel_file(path):
                    _install(registry, profile)
            except PanelError as exc:
                if warn:
                    warn(f"ignoring panel file: {exc}")

    for path in (extra_files or []):
        # Named explicitly, so a directory is accepted too and an error is
        # fatal rather than a warning.
        if os.path.isdir(path):
            for name in sorted(os.listdir(path)):
                if name.endswith(".json"):
                    for profile in load_panel_file(os.path.join(path, name)):
                        _install(registry, profile)
        else:
            for profile in load_panel_file(path):
                _install(registry, profile)
    return registry


def _install(registry, profile):
    """
    Put a profile into the registry, replacing any profile of that id.

    Aliases are INHERITED from whatever is being replaced when the new
    profile declares none of its own. A site that overrides
    'illumina-tso500' to record its validated limit of detection has not
    thereby decided that '--panel tso500' should stop working -- but that
    is what a plain replacement does, and the error it produces ("unknown
    panel") points at the name rather than at the override that broke it.
    """
    previous = registry.get(profile["id"])
    if previous and not profile.get("aliases"):
        profile["aliases"] = list(previous.get("aliases") or [])
    registry[profile["id"]] = profile
    return registry


def resolve_panel(name, registry):
    """
    Find one profile by id or alias, case-insensitively.

    Raises PanelError naming the near misses. Getting this wrong silently
    -- falling back to "no panel" when a name is misspelled -- would mean a
    run configured as though the operator had asked for nothing.
    """
    if not name:
        return None
    wanted = str(name).strip().lower()
    for profile in registry.values():
        if profile["id"].lower() == wanted:
            return profile
    for profile in registry.values():
        if wanted in [a.lower() for a in profile.get("aliases", [])]:
            return profile

    close = _closest(wanted, [p for p in registry])
    hint = f" Did you mean '{close}'?" if close else ""
    raise PanelError(
        f"unknown panel {name!r}.{hint} Run --list-panels to see what is "
        f"available, or --panel-file to add your own.")


# =============================================================================
# SECTION 4: APPLYING A PROFILE TO PARSED ARGUMENTS
# =============================================================================
# The precedence rule -- explicit flag beats profile beats default -- needs
# to know which flags were actually typed. argparse does not record that:
# an option left out and an option passed its own default value are
# indistinguishable in the namespace afterwards. So the command line is
# scanned for option strings before the profile is applied.

def explicit_dests(parser, argv):
    """
    Which of the parser's dests were named on this command line.

    Reads parser._actions, which is private but stable across every Python
    argparse this code can run on, and is the only place the option-string
    to dest mapping exists. The alternative -- re-declaring every flag's
    dest here -- is a second list to keep in step with the first, which is
    the kind of duplication that goes wrong quietly.
    """
    lookup = {}
    for action in parser._actions:        # noqa: SLF001 - see docstring
        for option in action.option_strings:
            lookup[option] = action.dest

    named = set()
    for token in argv:
        if not isinstance(token, str) or not token.startswith("-") \
                or token == "-" or token == "--":
            continue
        # --flag=value and -f=value both name the flag before the '='.
        head = token.split("=", 1)[0]
        if head in lookup:
            named.add(lookup[head])
    return named


def apply_profile(args, profile, explicit, only=None, rename=None):
    """
    Fill in everything the command line did not state, from the profile.

    ARGS:
        args:     the argparse namespace, modified in place.
        profile:  a validated profile dict.
        explicit: dests the user named on the command line (they win).
        only:     restrict to this set of keys, for a script that
                  implements a subset -- fastq_qc_clean.py applies the
                  trimming keys and has no use for the caller ones.
        rename:   {profile key: argparse dest} for the standalone tools,
                  which spell the same setting differently because they
                  are not prefixed by what they belong to: the pipeline's
                  --pcgr-assay is pcgr_report.py's --assay, and its
                  --coverage-min-depth is coverage_report.py's
                  --min-depth. Mapping here keeps ONE vocabulary in the
                  profiles rather than making a laboratory write the same
                  panel twice.

    RETURNS:
        (applied, overridden) -- two lists of (key, value) pairs. Both are
        printed and both go into the run manifest: a reader needs to see
        not only what the profile set but what it WANTED to set and was
        not allowed to, because the second list is where a surprising
        result usually comes from.
    """
    applied = []
    overridden = []
    rename = rename or {}

    for key, value in sorted((profile.get("settings") or {}).items()):
        if value is None:
            continue
        if only is not None and key not in only:
            continue
        # The dest this setting actually writes to on THIS script. Reported
        # under the profile's own name, because that is what the operator
        # would have to edit to change it.
        dest = rename.get(key, key)
        if not hasattr(args, dest):
            # The profile carries a setting this particular script does not
            # implement. Silence is right here: the same profile is handed
            # to three stages on purpose.
            continue

        if key == "skip_steps":
            # Merged, never replaced. A profile that skips dedup and a
            # command line that skips qc must produce a run that skips
            # both -- and the command line cannot "win" by listing fewer
            # steps, because omitting a step is not a statement about it.
            current = list(getattr(args, dest) or [])
            added = [s for s in value if s not in current]
            if added:
                setattr(args, dest, current + added)
                applied.append((key, added))
            continue

        if dest in explicit:
            overridden.append((key, value))
            continue

        setattr(args, dest, value)
        applied.append((key, value))

    return applied, overridden


def _render(value):
    """
    One setting as it should read on screen.

    Booleans become on/off rather than True/False: every boolean here
    backs a store_true flag, and "--pcgr-estimate-tmb False" invites the
    reader to think that string is passed to something.
    """
    if isinstance(value, bool):
        return "on" if value else "off"
    if isinstance(value, list):
        return " ".join(value)
    return str(value)


def format_application(profile, applied, overridden, rename=None):
    """
    The provenance block printed at the start of a run.

    Printed rather than merely recorded because the whole risk of a preset
    is that it configures something the operator did not see. A run that
    says "min_allele_fraction 0.05 (panel)" in its first ten lines cannot
    surprise anybody later.

    `rename` is the same map apply_profile() was given. Where a tool spells
    a setting differently from the profile, BOTH names are printed: the
    flag, because that is what this tool's reader is looking at, and the
    profile's key, because that is what they would have to edit to change
    it.
    """
    rename = rename or {}

    def flag(key):
        dest = rename.get(key, key)
        shown = "--" + dest.replace("_", "-")
        if dest != key:
            shown += f" (profile key: {key})"
        return shown

    lines = []
    lines.append(f"Panel   : {profile['name']} [{profile['id']}]")
    lines.append(f"          {profile['manufacturer']}, "
                 f"{profile['chemistry']}, from {profile['source']}")
    if applied:
        lines.append("          settings applied from this profile:")
        for key, value in applied:
            lines.append(f"            {flag(key)} {_render(value)}")
    else:
        lines.append("          no settings applied (the command line "
                     "stated them all)")
    if overridden:
        lines.append("          NOT applied -- your command line said "
                     "otherwise:")
        for key, value in overridden:
            lines.append(f"            {flag(key)} "
                         f"-- profile suggested {_render(value)}")
    for note in profile.get("notes", []):
        lines.append(f"          NOTE: {note}")
    return "\n".join(lines)


# =============================================================================
# SECTION 5: MEASURING THE TARGET -- THE TMB DENOMINATOR
# =============================================================================
# PCGR divides the mutation count by an assumed 34 Mb for a TARGETED assay
# unless it is told otherwise, so a 2 Mb panel reports a TMB roughly 17x
# too low and says nothing about it. The vendor's published figure is one
# answer; the BED that is actually being called over is a better one,
# because it reflects the spike-ins, the version and the build. So the
# footprint is measured, and the published figure is demoted to a
# cross-check.

def bed_footprint(path):
    """
    Total unique bases covered by a BED, and the number of regions.

    Overlapping and adjacent regions are merged before summing: vendor
    BEDs routinely contain overlapping probe tiles, and summing them
    unmerged inflates the footprint -- which would deflate TMB by the same
    factor, in the same silent direction as the bug this is here to fix.

    RETURNS:
        (bases, regions, contigs) or None if the file cannot be read.
    """
    by_contig = {}
    regions = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t") if "\t" in line else line.split()
                if len(fields) < 3:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue              # a header row, or not a BED
                if end <= start:
                    continue
                by_contig.setdefault(fields[0], []).append((start, end))
                regions += 1
    except OSError:
        return None

    total = 0
    for spans in by_contig.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:          # overlapping or touching
                cur_end = max(cur_end, end)
            else:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
        total += cur_end - cur_start
    return total, regions, sorted(by_contig)


def footprint_mb(path):
    """The BED's merged footprint in megabases, or None."""
    measured = bed_footprint(path)
    if not measured or not measured[0]:
        return None
    return round(measured[0] / 1e6, 4)


def target_size_note(measured_mb, nominal_mb, tolerance=0.25):
    """
    One line comparing the measured footprint with the profile's figure.

    A large disagreement almost always means the wrong BED: the previous
    kit version, the other genome build, or the capture-target file where
    the primary-target file was wanted. Saying so at the start of the run
    is worth a great deal more than saying nothing and dividing by it.
    """
    if measured_mb is None:
        return None
    if nominal_mb in (None, 0):
        return (f"Target footprint measured from the BED: "
                f"{measured_mb:.3f} Mb (this is the TMB denominator).")
    drift = abs(measured_mb - nominal_mb) / float(nominal_mb)
    if drift <= tolerance:
        return (f"Target footprint {measured_mb:.3f} Mb, consistent with "
                f"the profile's nominal {nominal_mb:.3f} Mb.")
    return (f"Target footprint {measured_mb:.3f} Mb differs from this "
            f"profile's nominal {nominal_mb:.3f} Mb by "
            f"{drift * 100:.0f}%. That usually means a BED from a "
            f"different kit version or genome build. The measured value is "
            f"what will be used -- check it is the right file.")


def tmb_advice(measured_mb, estimate_tmb):
    """
    Whether a TMB request makes sense at this footprint.

    Returns a warning string or None. TMB from a small panel is not a
    weaker version of TMB from an exome: it is a count over a denominator
    small enough that ordinary artefact variation swamps the signal.
    """
    if not estimate_tmb or measured_mb is None:
        return None
    if measured_mb >= TMB_MIN_PANEL_MB:
        return None
    return (f"TMB was requested over {measured_mb:.3f} Mb of target, below "
            f"the ~{TMB_MIN_PANEL_MB:g} Mb at which a panel TMB starts to "
            f"be interpretable. A single extra artefact moves this figure "
            f"by several mutations per Mb; do not report it.")


# =============================================================================
# SECTION 6: CONTIG NAMING -- THE VENDOR BED THAT MATCHES NOTHING
# =============================================================================
# Vendor BEDs arrive in whichever naming their designer used. Ensembl
# writes "1" and "MT"; UCSC -- and the hg38 this pipeline installs --
# writes "chr1" and "chrM". Hand GATK a mismatched BED and it aborts with a
# sequence-dictionary error, which is at least loud. Hand the SAME file to
# the coverage check and every region simply finds no reads, so the report
# says the panel was not covered at all: a confident, entirely wrong
# clinical statement. This is the COSMIC contig problem again, one file
# along, and it is handled the same way -- detect it and translate.

def _read_bed_contigs(path, limit=200000):
    """Distinct contig names in a BED, bounded so a huge file stays cheap."""
    names = []
    seen = set()
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for index, line in enumerate(fh):
                if index > limit:
                    break
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                name = (line.split("\t") if "\t" in line
                        else line.split())[0]
                if name not in seen:
                    seen.add(name)
                    names.append(name)
    except OSError:
        return []
    return names


def reference_contigs(fasta):
    """
    Contig names from the reference's .fai, or [] if it is not indexed yet.

    The .fai is the cheapest authoritative statement of what the reference
    calls its chromosomes; parsing the FASTA itself would mean reading
    three gigabytes to learn twenty-five strings.
    """
    fai = fasta + ".fai"
    names = []
    try:
        with open(fai, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.strip():
                    names.append(line.split("\t")[0])
    except OSError:
        return []
    return names


def _translate_contig(name, to_ucsc):
    """One contig name converted between UCSC and Ensembl conventions."""
    if to_ucsc:
        if name.startswith("chr"):
            return name
        if name in ("MT", "M"):
            return "chrM"
        return "chr" + name
    if not name.startswith("chr"):
        return name
    if name == "chrM":
        return "MT"
    return name[3:]


def harmonise_bed_contigs(bed, fasta, out_dir, dry_run=False):
    """
    Make a BED's contig names match the reference, writing a copy if needed.

    Never edits the input: vendor BEDs are controlled documents and a
    laboratory's copy must stay byte-identical to the one in its records.
    A translated copy is written into the run's own directory instead, so
    the run is self-describing about what it actually called over.

    RETURNS:
        (path_to_use, note). path_to_use is the original when nothing
        needed doing. note is a line for the log, or None when the file
        already matched.
    """
    if not bed or not os.path.isfile(bed):
        return bed, None

    ref_names = reference_contigs(fasta)
    if not ref_names:
        # Not indexed yet (or a dry run before indexing). Silence rather
        # than a guess: an unverifiable claim about contig naming is worse
        # than no claim.
        return bed, None

    bed_names = _read_bed_contigs(bed)
    if not bed_names:
        return bed, f"{os.path.basename(bed)} contains no usable BED records."

    ref_set = set(ref_names)
    matched = [n for n in bed_names if n in ref_set]
    if matched:
        # Some contigs already agree. Partial agreement is normal (a BED
        # may carry alt contigs the reference lacks) and is not a naming
        # problem.
        return bed, None

    to_ucsc = any(n.startswith("chr") for n in ref_names)
    translated = {n: _translate_contig(n, to_ucsc) for n in bed_names}
    usable = [n for n, t in translated.items() if t in ref_set]
    if not usable:
        return bed, (
            f"{os.path.basename(bed)} names contigs the reference does not "
            f"have ({', '.join(bed_names[:4])}...), and renaming them does "
            f"not help either. This is the wrong BED for this genome build, "
            f"or the wrong reference. Nothing will be called over it.")

    style = "UCSC (chr1, chrM)" if to_ucsc else "Ensembl (1, MT)"
    destination = os.path.join(
        out_dir, os.path.basename(bed).replace(".bed", "") + ".renamed.bed")
    note = (f"{os.path.basename(bed)} uses the other contig convention from "
            f"the reference; wrote a copy renamed to {style} at "
            f"{destination}. The original is untouched.")
    if dry_run:
        return destination, note + " (dry run: not written)"

    os.makedirs(out_dir, exist_ok=True)
    dropped = 0
    try:
        with open(bed, encoding="utf-8", errors="replace") as src, \
                open(destination, "w", encoding="utf-8") as dst:
            for line in src:
                if not line.strip() or line.startswith(("#", "track",
                                                        "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    fields = line.split()
                if len(fields) < 3:
                    continue
                new_name = translated.get(fields[0],
                                          _translate_contig(fields[0],
                                                            to_ucsc))
                if new_name not in ref_set:
                    # A contig the reference genuinely does not have -- an
                    # alt or a decoy. Dropping it is right: calling over it
                    # is impossible and leaving it in aborts GATK.
                    dropped += 1
                    continue
                fields[0] = new_name
                dst.write("\t".join(fields) + "\n")
    except OSError as exc:
        return bed, (f"could not write a renamed copy of "
                     f"{os.path.basename(bed)} ({exc}); using it as-is, "
                     f"which will fail at the caller.")
    if dropped:
        note += (f" {dropped} record(s) named contigs absent from the "
                 f"reference and were dropped.")
    return destination, note


# =============================================================================
# SECTION 7: PRESENTATION AND AUTHORING
# =============================================================================

def format_panel_list(registry):
    """--list-panels: every profile, grouped by manufacturer."""
    by_vendor = {}
    for profile in registry.values():
        by_vendor.setdefault(profile["manufacturer"], []).append(profile)

    lines = ["", "AVAILABLE PANEL PROFILES", "=" * 72, ""]
    # "any" holds the generic shapes and belongs first: it is where
    # somebody with an unlisted kit should start.
    order = (["any"] if "any" in by_vendor else []) + \
        sorted(v for v in by_vendor if v != "any")
    for vendor in order:
        lines.append(f"{vendor}")
        for profile in sorted(by_vendor[vendor], key=lambda p: p["id"]):
            size = profile.get("nominal_target_size_mb")
            size_text = f", ~{size:g} Mb" if size else ""
            flag = " *" if profile["source"] != "built-in" else ""
            # The assay type is spelled out rather than left to be inferred
            # from the chemistry string: it decides which pipeline runs,
            # which scripts exist, and which reference data is needed.
            kind = assay_type(profile).upper()
            lines.append(f"  {profile['id']:<34s} [{kind}] "
                         f"{profile['chemistry']}{size_text}{flag}")
            lines.append(f"  {'':<34s} {profile['summary']}")
        lines.append("")
    lines.append("  * = loaded from your own panel file, not built in.")
    lines.append("")
    lines.append("  --describe-panel ID   everything a profile sets, and why")
    lines.append("  --new-panel-template FILE   a starting point for a kit")
    lines.append("                              that is not listed")
    lines.append("")
    lines.append("  DNA profiles never supply the target BED. Pass your "
                 "kit's own with")
    lines.append("  --panel-bed, which sets both --intervals and "
                 "--coverage-bed.")
    lines.append("")
    lines.append("  RNA profiles run a DIFFERENT pipeline -- align_rna.py "
                 "and fusion_calling.py,")
    lines.append("  not the somatic caller. They need a GENCODE GTF and a "
                 "STAR index rather")
    lines.append("  than a BED; see RNA.md.")
    return "\n".join(lines)


def format_panel(profile):
    """--describe-panel: one profile in full, including its caveats."""
    kind = assay_type(profile)
    lines = ["", f"{profile['name']}  [{profile['id']}]", "=" * 72,
             f"Manufacturer : {profile['manufacturer']}",
             f"Assay        : {kind.upper()}"
             + ("  -- runs align_rna.py + fusion_calling.py, NOT the "
                "somatic caller" if kind == "rna" else ""),
             f"Chemistry    : {profile['chemistry']}",
             f"Source       : {profile['source']}"]
    if profile.get("aliases"):
        lines.append(f"Also known as: {', '.join(profile['aliases'])}")
    if profile.get("genome_builds"):
        lines.append(f"Builds       : {', '.join(profile['genome_builds'])}")
    size = profile.get("nominal_target_size_mb")
    if size:
        lines.append(f"Nominal size : {size:g} Mb  (cross-check only -- the "
                     f"BED you supply is measured and used)")
    if profile.get("summary"):
        lines += ["", profile["summary"]]

    if profile.get("requires"):
        lines += ["", "YOU MUST STILL SUPPLY"]
        for item in profile["requires"]:
            lines.append(f"  - {item}")

    settings = profile.get("settings") or {}
    if settings:
        lines += ["", "SETTINGS APPLIED (unless your command line says "
                      "otherwise)"]
        for key in sorted(settings):
            value = settings[key]
            if value is None:
                continue
            lines.append(f"  --{key.replace('_', '-'):<28s} "
                         f"{_render(value)}")

    if profile.get("notes"):
        lines += ["", "NOTES"]
        for note in profile["notes"]:
            lines.append(f"  - {note}")
    lines.append("")
    return "\n".join(lines)


def panel_template(base=None):
    """
    A starting-point profile for a kit that is not built in.

    Deliberately filled in rather than blank: the fastest way to get a new
    manufacturer's panel running is to copy the closest generic shape and
    change the four or five values that differ, and a template full of
    nulls does not show which those are.
    """
    template = {
        "id": "my-panel",
        "name": "My laboratory's panel",
        "manufacturer": "Vendor name",
        "chemistry": "hybrid-capture",
        "aliases": [],
        "genome_builds": ["hg38"],
        "nominal_target_size_mb": None,
        "summary": "One line describing the kit.",
        "requires": ["panel BED"],
        "notes": [
            "Record here what a later reader would otherwise have to "
            "rediscover: the validated limit of detection, which BED "
            "version this is, anything the vendor pipeline already did to "
            "the reads."
        ],
        "settings": dict(_capture_defaults()),
    }
    if base:
        template = copy.deepcopy(base)
        template["id"] = base["id"] + "-local"
        template["name"] = base["name"] + " (local copy)"
        template.pop("source", None)
        template.setdefault("notes", []).insert(
            0, f"Copied from the built-in profile '{base['id']}'. Adjust to "
               f"your validated assay.")
    return template


def write_panel_template(path, base=None):
    """Write a template profile and return the path, for --new-panel-template."""
    data = panel_template(base)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    return path


# Steps that say where THIS RUN's data came from rather than what the kit
# is. They are stripped when a profile is saved: a laboratory that once
# resumed from an existing QC manifest must not end up with a profile that
# silently skips quality control on every future run of that assay.
_RUN_LOCAL_SKIP_STEPS = frozenset({"qc", "align"})


def profile_from_args(args, panel_id, name=None, manufacturer="local",
                      chemistry="hybrid-capture", notes=None,
                      keep_target_size=False):
    """
    Capture a run's effective settings as a reusable profile.

    This is the other half of customisation: tune a run with ordinary
    flags and --dry-run until it is right, then --save-panel-as writes what
    you arrived at as a profile, so the next operator names it instead of
    reconstructing it.

    EVERY effective value is written out, including ones that merely
    happened to be at a script default. That is deliberate: a profile is
    meant to reproduce a configuration, and a value left out would quietly
    follow whatever the default becomes in a later version.

    Two kinds of value are excluded, because both would be WRONG rather
    than merely verbose anywhere but the machine that wrote them:

      * paths -- a profile carrying a BED path configures a different
        target, or none at all, on any other machine;
      * the TMB denominator, when it was measured from this run's BED --
        it belongs to that design, not to the assay. It is written to
        nominal_target_size_mb instead, where it serves as the cross-check
        it actually is.
    """
    settings = {}
    for key, kind in PANEL_SETTING_TYPES.items():
        if kind == "path":
            continue
        if key == "pcgr_target_size_mb" and not keep_target_size:
            continue
        value = getattr(args, key, None)
        # `value in (None, "", [], False)` looks equivalent and is not:
        # 0 == False in Python, so that form silently drops
        # --interval-padding 0, which is exactly the value an amplicon
        # profile needs to state. Test the empty cases by identity.
        if value is None or value is False:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        if isinstance(value, list) and not value:
            continue
        if key == "skip_steps":
            value = [s for s in value if s not in _RUN_LOCAL_SKIP_STEPS]
            if not value:
                continue
        # Note on the booleans: a store_true flag that is off is
        # indistinguishable from one that was never given, so an "off"
        # toggle is not written out. It does not need to be -- off is the
        # script's default for every one of them, so a profile that omits
        # it produces the same run.
        settings[key] = list(value) if isinstance(value, list) else value
    return {
        "id": panel_id,
        "name": name or panel_id,
        "manufacturer": manufacturer,
        "chemistry": chemistry,
        "aliases": [],
        "genome_builds": [],
        "nominal_target_size_mb": (None if keep_target_size
                                   else getattr(args, "pcgr_target_size_mb",
                                                None)),
        "summary": "Saved from a run's effective settings.",
        "requires": ["panel BED"],
        "notes": notes or [
            "Saved by --save-panel-as. Target BEDs are not stored in a "
            "profile: pass --panel-bed when you use it."
        ],
        "settings": settings,
    }


def save_profile(profile, path):
    """Write one profile to disk, validating it first."""
    _validate_profile(profile, path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(profile, fh, indent=2)
        fh.write("\n")
    return path


# =============================================================================
# SECTION 8: THE SHARED FRONT END FOR EVERY SCRIPT
# =============================================================================
# --list-panels and --describe-panel have to work WITHOUT the required
# arguments of the script hosting them: somebody asking "which panels do
# you know about?" has not chosen an output directory yet. argparse cannot
# express that, so these are intercepted from sys.argv before the real
# parser runs. The flags are still declared on the real parser, so they
# appear in --help and are rejected if misspelled.

def handle_panel_queries(argv=None):
    """
    Answer --list-panels / --describe-panel / --new-panel-template and exit.

    Called first thing in main(), before the script's own parser runs.
    Returns normally when none of those flags is present.
    """
    argv = list(sys.argv[1:] if argv is None else argv)

    def value_after(flag):
        for index, token in enumerate(argv):
            if token == flag and index + 1 < len(argv):
                return argv[index + 1]
            if token.startswith(flag + "="):
                return token.split("=", 1)[1]
        return None

    files = []
    for index, token in enumerate(argv):
        if token == "--panel-file" and index + 1 < len(argv):
            files.append(argv[index + 1])
        elif token.startswith("--panel-file="):
            files.append(token.split("=", 1)[1])

    wants_list = "--list-panels" in argv
    describe = value_after("--describe-panel")
    template = value_after("--new-panel-template")
    if not (wants_list or describe or template):
        return

    try:
        registry = available_panels(
            extra_files=files,
            warn=lambda msg: print(f"[WARN] {msg}", file=sys.stderr))
        if wants_list:
            print(format_panel_list(registry))
        if describe:
            print(format_panel(resolve_panel(describe, registry)))
        if template:
            base = None
            source = value_after("--panel")
            if source:
                base = resolve_panel(source, registry)
            written = write_panel_template(template, base)
            print(f"Wrote a panel profile template to {written}.")
            print("Edit it, then use it with:")
            print(f"  --panel-file {written} --panel "
                  f"{panel_template(base)['id']}")
    except PanelError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    # Standalone use: `python panel_profiles.py --list-panels` works, and
    # the installer's --check runs `--help` against every helper module, so
    # this file must start cleanly on its own.
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        print("  --list-panels              list every known profile")
        print("  --describe-panel ID        one profile in full")
        print("  --new-panel-template FILE  write a starting point")
        print("  --panel-file FILE          also load profiles from FILE")
        sys.exit(0)
    handle_panel_queries()
