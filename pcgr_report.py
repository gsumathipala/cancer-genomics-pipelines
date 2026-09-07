#!/usr/bin/env python3
"""
pcgr_report.py
==============
Run PCGR (Personal Cancer Genome Reporter) on a somatic VCF to produce a
clinical-grade interpretation report.

WHAT PCGR ADDS THAT THIS PIPELINE DOES NOT ALREADY HAVE
--------------------------------------------------------
  comprehensive_variant_calling.py already gives you gene/consequence
  (SnpEff) and known-mutation identity (COSMIC). PCGR adds the
  *interpretation* layer on top of a somatic call set:

    - Therapeutic actionability, matched against CIViC / CGI and tiered by
      the AMP/ASCO/CAP guidelines (Tier 1-4).
    - Tumour mutational burden (TMB), with the assay-size correction that
      makes a panel TMB comparable to a WES TMB.
    - Microsatellite instability (MSI) classification.
    - Mutational signature fitting against the COSMIC SBS catalogue --
      this is where HRD (SBS3, a PARP-inhibitor marker), MMR deficiency,
      APOBEC and POLE show up. For an FFPE panel it also surfaces the
      artefact signatures, which doubles as a QC signal.
    - A self-contained interactive HTML report plus machine-readable TSV.

  PCGR runs Ensembl VEP itself as part of its own workflow, so it does not
  need (or read) the SnpEff annotation this pipeline produces. That is why
  the pipeline hands PCGR the *filtered* (or COSMIC-annotated) VCF rather
  than the SnpEff output -- see the note on INPUT SELECTION below.

REQUIREMENTS (none of which this module installs for you)
---------------------------------------------------------
  1. The `pcgr` executable on PATH. PCGR is NOT on bioconda -- it ships
     from its own 'pcgr' Anaconda channel and needs TWO sibling conda
     environments built from version-pinned lock files:
         LOCK=https://raw.githubusercontent.com/sigven/pcgr/v2.3.2/conda/env/lock
         conda create -y -n pcgr  --file <(curl -sSL $LOCK/pcgr-linux-64.lock)
         conda create -y -n pcgrr --file <(curl -sSL $LOCK/pcgrr-linux-64.lock)
     `pcgr` shells out to `pcgrr`, which it locates as a sibling of the
     activated env (dirname($CONDA_PREFIX)/pcgrr).
  2. A PCGR reference data bundle, downloaded separately and matched to the
     installed PCGR version and genome build (~7 GB). Point at it with
     --pcgr-refdata-dir.
  3. An Ensembl VEP cache directory (--vep-dir), matched to the VEP version
     PCGR ships (~24 GB). Pass the PARENT of homo_sapiens/, not the cache
     directory itself.

  RUN THIS MODULE FROM THE PCGR ENVIRONMENT, not the pipeline's. PCGR
  resolves its VEP plugin directory from $CONDA_PREFIX/share, so merely
  putting pcgr's bin/ on PATH is not enough -- a real run then dies with
  "FileNotFoundError: No ensembl-vep directories found". This module is
  standard-library-only, so PCGR's own Python runs it happily:
      conda activate pcgr && python pcgr_report.py ...

  If `pcgr` is not on PATH this module reports "skipped" rather than
  failing, matching how the pipeline treats SnpEff and bcftools: PCGR is
  enrichment, not a load-bearing stage.

INPUT SELECTION (why not the SnpEff VCF?)
------------------------------------------
  PCGR performs its own VEP annotation and validates the INFO fields of the
  VCF it is given. Feeding it an already-annotated VCF adds nothing it will
  use and risks tripping that validation. The pipeline therefore passes the
  best *un-over-annotated* call set -- the FilterMutectCalls output, or the
  COSMIC-annotated one when COSMIC ran -- while still running this step
  after SnpEff in pipeline order. Override with --input-vcf when running
  standalone.

!! DEPTH / ALLELE-FREQUENCY TAGS -- READ THIS BEFORE TRUSTING TMB !!
--------------------------------------------------------------------
  PCGR reads sequencing depth and allele fraction from *INFO* tags, named
  via --tumor-dp-tag / --tumor-af-tag. GATK Mutect2 does NOT write those to
  INFO: it writes AD/AF/DP as per-sample *FORMAT* fields. So out of the box
  there is nothing for PCGR to read.

  Consequences if you ignore this:
    - PCGR still runs and still produces a report.
    - But depth/AF-based filtering is inert, so TMB is computed over
      unfiltered calls and will read HIGH. Do not report that number.

  THE FIX: pass --lift-tags (or --pcgr-lift-tags from the pipeline).
  lift_format_to_info() below rewrites the VCF, copying the tumour's (and
  the normal's) FORMAT DP/AF up into INFO tags TDP/TAF (and NDP/NAF), and
  the tag names are then handed to PCGR automatically. The input VCF is
  never modified -- a new file is written alongside the report.

  Guarantees worth knowing when you are debugging a suspicious TMB:
    - The tumour column is identified from Mutect2's own ##tumor_sample
      header, not by column position, which is not guaranteed.
    - If AF is absent the value is derived from AD as alt/(ref+alt).
    - Multiallelic records contribute their FIRST ALT only, and are
      counted and warned about; split them with `bcftools norm -m-any`
      first if that matters.
    - Records with no usable DP/AF are left untagged rather than
      guessed at, and are counted.
    - Re-running over an already-lifted VCF is refused, so you cannot
      produce duplicate INFO keys.
    - If the transform cannot run at all, PCGR still runs on the
      original VCF and warns that TMB is unfiltered.

  Without --lift-tags, tag names you pass are still forwarded verbatim,
  so a VCF you transformed yourself works too.

VERSION DRIFT
-------------
  PCGR's CLI changed between the 1.x and 2.x series; notably the reference
  bundle argument was --pcgr_dir in 1.x and --refdata_dir in 2.x. This
  module emits the 2.x spelling by default; pass --pcgr-legacy-v1 for the
  1.x spelling. Anything else that moves can be forwarded verbatim with
  --pcgr-extra-args, so a CLI change does not require editing this file.

Usage (standalone)
------------------
  python pcgr_report.py \\
      --input-vcf results/mutect2/TUMOUR_01.mutect2.filtered.vcf.gz \\
      --output-dir results/pcgr \\
      --sample-id TUMOUR_01 \\
      --pcgr-refdata-dir /data/pcgr_bundle \\
      --vep-dir /data/vep_cache \\
      --genome-assembly grch38 \\
      --assay TARGETED \\
      --tumour-only
"""

import argparse
import glob
import gzip
import os
import shlex
import shutil
import subprocess
import sys


# =============================================================================
# CONSTANTS
# =============================================================================
# PCGR's --assay vocabulary. TARGETED is the right choice for a capture
# panel; it changes how TMB is scaled.
ASSAY_CHOICES = ("WGS", "WES", "TARGETED")

# PCGR's --genome_assembly vocabulary.
GENOME_CHOICES = ("grch38", "grch37")


# =============================================================================
# SECTION 1: COMMAND EXECUTION
# =============================================================================
# Duplicated from the other pipeline scripts on purpose: every script in
# this directory is runnable on its own, without importing its siblings.

def run_command(cmd, tag, log_path=None, dry_run=False):
    """
    Execute a command, streaming output to the console and optional log.

    Returns the exit code. Never raises on a non-zero exit -- the caller
    decides whether a failure is fatal.
    """
    printable = shlex.join(cmd)
    print(f"\n[{tag}] $ {printable}", flush=True)

    if dry_run:
        print(f"[{tag}] (dry run - not executed)")
        return 0

    log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if log_fh:
            log_fh.write(f"\n### {tag}\n### {printable}\n")

        # PCGR is long-running (VEP annotation dominates), so stream its
        # output line by line rather than buffering until it exits.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        for line in proc.stdout:
            line = line.rstrip("\n")
            print(f"    | {line}", flush=True)
            if log_fh:
                log_fh.write(line + "\n")
        proc.wait()

        if proc.returncode != 0:
            print(f"[{tag}] ERROR: exit code {proc.returncode}")
        return proc.returncode

    except FileNotFoundError:
        print(f"[{tag}] ERROR: '{cmd[0]}' not found on PATH.")
        return 127
    except OSError as exc:
        print(f"[{tag}] ERROR: {exc}")
        return 1
    finally:
        if log_fh:
            log_fh.close()


# =============================================================================
# SECTION 2: FORMAT -> INFO TAG TRANSFORM
# =============================================================================
# PCGR reads depth and allele fraction from INFO tags. GATK Mutect2 writes
# them as per-sample FORMAT fields (AD/AF/DP), so without this transform
# PCGR has nothing to filter on and its TMB is computed over unfiltered
# calls. This lifts the tumour's (and optionally the normal's) values up
# into INFO.
#
# Implemented in plain Python rather than a bcftools query/annotate recipe
# on purpose: no external dependency, no temporary TSV + bgzip + tabix
# dance, and it can be tested directly against synthetic VCFs.

def _open_maybe_gzip(path, mode="rt"):
    """Open a VCF whether or not it is gzip/bgzip compressed."""
    # Sniff the gzip magic bytes rather than trusting the extension: a
    # bgzipped VCF is often named .vcf.gz but a plain one sometimes is too.
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _parse_sample_value(fmt_keys, sample_field):
    """Zip a record's FORMAT keys with one sample's colon-separated values."""
    values = sample_field.split(":")
    # A trailing FORMAT key may be absent from the sample field; zip stops
    # at the shorter of the two, which is the behaviour we want.
    return dict(zip(fmt_keys, values))


def _depth_and_af(sample_values):
    """
    Extract (depth, allele_fraction) for one sample from its FORMAT values.

    Returns (None, None) for either value that cannot be determined.

    AF is taken from the AF field when present. Mutect2 declares AF with
    Number=A, so it carries one value per ALT allele; we take the FIRST
    ALT. Records with multiple ALTs should be split beforehand
    (bcftools norm -m-any) -- the caller counts them so you can tell.

    When AF is absent we fall back to computing it from AD
    (ref_depth, alt1_depth, ...), which is what AF means anyway.
    """
    depth = None
    raw_dp = sample_values.get("DP")
    if raw_dp and raw_dp != ".":
        try:
            depth = int(raw_dp)
        except ValueError:
            depth = None

    af = None
    raw_af = sample_values.get("AF")
    if raw_af and raw_af != ".":
        first = raw_af.split(",")[0]
        try:
            af = float(first)
        except ValueError:
            af = None

    if af is None:
        raw_ad = sample_values.get("AD")
        if raw_ad and raw_ad != ".":
            try:
                counts = [int(x) for x in raw_ad.split(",") if x != "."]
            except ValueError:
                counts = []
            if len(counts) >= 2 and sum(counts) > 0:
                af = counts[1] / float(sum(counts))
                # AD also gives us a depth if DP was missing.
                if depth is None:
                    depth = sum(counts)

    return depth, af


def lift_format_to_info(vcf_in, vcf_out, tumour_sample=None,
                        normal_sample=None, tumour_dp_tag="TDP",
                        tumour_af_tag="TAF", normal_dp_tag="NDP",
                        normal_af_tag="NAF", dry_run=False):
    """
    Copy per-sample FORMAT depth/allele-fraction into INFO tags.

    The tumour and normal columns are identified from Mutect2's own
    ``##tumor_sample=`` / ``##normal_sample=`` header lines, which is the
    only reliable source -- column order is not guaranteed. Explicit
    arguments override them; failing both, a single-sample VCF is treated
    as the tumour.

    The input is never modified: a new VCF is written to `vcf_out`.

    RETURNS:
        dict of counters on success, or None if the transform could not be
        performed (the caller then falls back to running PCGR untagged).
    """
    if dry_run:
        print(f"[DRY RUN] Would lift FORMAT DP/AF into INFO "
              f"({tumour_dp_tag}/{tumour_af_tag}) -> {vcf_out}")
        return {"dry_run": True}

    if not os.path.exists(vcf_in):
        print(f"[WARN] Cannot lift tags: input VCF not found: {vcf_in}")
        return None

    header_lines = []
    hdr_tumour = None
    hdr_normal = None
    chrom_line = None
    samples = []

    # --- Pass 1: read the header, learn the sample layout ---------------
    try:
        with _open_maybe_gzip(vcf_in) as fh:
            for line in fh:
                if line.startswith("##"):
                    header_lines.append(line.rstrip("\n"))
                    if line.startswith("##tumor_sample="):
                        hdr_tumour = line.strip().split("=", 1)[1]
                    elif line.startswith("##normal_sample="):
                        hdr_normal = line.strip().split("=", 1)[1]
                elif line.startswith("#CHROM"):
                    chrom_line = line.rstrip("\n")
                    samples = chrom_line.split("\t")[9:]
                    break
                else:
                    break
    except OSError as exc:
        print(f"[WARN] Cannot lift tags: {exc}")
        return None

    if chrom_line is None:
        print("[WARN] Cannot lift tags: no #CHROM header line found.")
        return None
    if not samples:
        print("[WARN] Cannot lift tags: VCF has no sample columns.")
        return None

    # Refuse to run twice over the same file rather than emit duplicate
    # INFO keys, which would make the VCF invalid.
    for tag in (tumour_dp_tag, tumour_af_tag, normal_dp_tag, normal_af_tag):
        if any(h.startswith(f"##INFO=<ID={tag},") for h in header_lines):
            print(f"[WARN] Cannot lift tags: INFO/{tag} is already declared "
                  f"in {vcf_in}. Refusing to write duplicate INFO keys.")
            return None

    # --- Resolve which column is the tumour, and which the normal -------
    tumour = tumour_sample or hdr_tumour
    normal = normal_sample or hdr_normal
    if tumour is None and len(samples) == 1:
        tumour = samples[0]  # Unambiguous: the only sample is the tumour.
    if tumour is None:
        print("[WARN] Cannot lift tags: no ##tumor_sample header and no "
              "--tumour-sample given, so the tumour column is ambiguous.")
        return None
    if tumour not in samples:
        print(f"[WARN] Cannot lift tags: tumour sample '{tumour}' is not "
              f"among the VCF's samples {samples}.")
        return None
    t_col = 9 + samples.index(tumour)
    n_col = 9 + samples.index(normal) if normal in samples else None

    # --- Build the INFO header declarations we are about to add ---------
    new_headers = [
        f'##INFO=<ID={tumour_dp_tag},Number=1,Type=Integer,'
        f'Description="Tumour sequencing depth, lifted from FORMAT/DP of '
        f'sample {tumour} by pcgr_report.py">',
        f'##INFO=<ID={tumour_af_tag},Number=1,Type=Float,'
        f'Description="Tumour allele fraction (first ALT), lifted from '
        f'FORMAT/AF of sample {tumour} by pcgr_report.py">',
    ]
    if n_col is not None:
        new_headers += [
            f'##INFO=<ID={normal_dp_tag},Number=1,Type=Integer,'
            f'Description="Normal sequencing depth, lifted from FORMAT/DP of '
            f'sample {normal} by pcgr_report.py">',
            f'##INFO=<ID={normal_af_tag},Number=1,Type=Float,'
            f'Description="Normal allele fraction (first ALT), lifted from '
            f'FORMAT/AF of sample {normal} by pcgr_report.py">',
        ]

    stats = {"records": 0, "tumour_tagged": 0, "normal_tagged": 0,
             "multiallelic": 0, "no_depth_or_af": 0,
             "tumour_sample": tumour, "normal_sample": normal if n_col else None}

    os.makedirs(os.path.dirname(os.path.abspath(vcf_out)) or ".",
                exist_ok=True)

    # --- Pass 2: rewrite ------------------------------------------------
    try:
        with _open_maybe_gzip(vcf_in) as src, \
                open(vcf_out, "w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("##"):
                    dst.write(line)
                    continue
                if line.startswith("#CHROM"):
                    # New INFO declarations go immediately before #CHROM,
                    # which is always valid VCF.
                    for h in new_headers:
                        dst.write(h + "\n")
                    dst.write(line)
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) <= t_col:
                    dst.write(line)  # Malformed/truncated: pass through.
                    continue

                stats["records"] += 1
                if "," in fields[4]:
                    stats["multiallelic"] += 1

                fmt_keys = fields[8].split(":")
                additions = []

                t_vals = _parse_sample_value(fmt_keys, fields[t_col])
                t_dp, t_af = _depth_and_af(t_vals)
                if t_dp is not None:
                    additions.append(f"{tumour_dp_tag}={t_dp}")
                if t_af is not None:
                    additions.append(f"{tumour_af_tag}={t_af:.6g}")
                if t_dp is not None or t_af is not None:
                    stats["tumour_tagged"] += 1
                else:
                    stats["no_depth_or_af"] += 1

                if n_col is not None and len(fields) > n_col:
                    n_vals = _parse_sample_value(fmt_keys, fields[n_col])
                    n_dp, n_af = _depth_and_af(n_vals)
                    if n_dp is not None:
                        additions.append(f"{normal_dp_tag}={n_dp}")
                    if n_af is not None:
                        additions.append(f"{normal_af_tag}={n_af:.6g}")
                    if n_dp is not None or n_af is not None:
                        stats["normal_tagged"] += 1

                if additions:
                    # An INFO column of "." means "no fields"; appending to
                    # it literally would produce ".;TDP=30", which is invalid.
                    fields[7] = (";".join(additions) if fields[7] == "."
                                 else fields[7] + ";" + ";".join(additions))
                dst.write("\t".join(fields) + "\n")
    except OSError as exc:
        print(f"[WARN] Cannot lift tags: {exc}")
        return None

    print(f"[INFO] Lifted FORMAT->INFO for {stats['tumour_tagged']}/"
          f"{stats['records']} records "
          f"(tumour='{tumour}'"
          + (f", normal='{normal}'" if n_col is not None else "")
          + f") -> {vcf_out}")
    if stats["multiallelic"]:
        print(f"[WARN] {stats['multiallelic']} multiallelic record(s): only "
              f"the first ALT's AF was lifted. Split them first with "
              f"'bcftools norm -m-any' if that matters for your TMB.")
    if stats["no_depth_or_af"]:
        print(f"[WARN] {stats['no_depth_or_af']} record(s) had neither DP nor "
              f"AF/AD for the tumour and were left untagged.")
    return stats


# =============================================================================
# SECTION 3: COMMAND CONSTRUCTION
# =============================================================================

def pcgr_reserved_info_tags(refdata_dir, genome_assembly="grch38"):
    """
    Read the INFO tag names PCGR reserves for its own annotations.

    PCGR refuses to start if the input VCF already carries any INFO tag it
    intends to write, with

        INFO tag <X> in the input VCF coincides with a VCF annotation tag
        produced by PCGR/CPSR - please remove or rename this tag

    The authoritative list ships inside the reference bundle rather than
    the code, so it is read from there instead of being hardcoded: it
    tracks the installed bundle across upgrades, and hardcoding it would
    quietly rot.

    Returns an empty set if the files are missing or unreadable, which
    turns the stripping step below into a no-op -- PCGR's own check still
    catches the collision, just with its usual error.

    ARGS:
        refdata_dir:     PCGR reference bundle directory.
        genome_assembly: 'grch38' or 'grch37'.

    RETURNS:
        set: Reserved INFO tag names.
    """
    tags = set()
    base = os.path.join(refdata_dir, "data", genome_assembly)
    for name in ("vcf_infotags_vep.tsv", "vcf_infotags_other.tsv"):
        path = os.path.join(base, name)
        try:
            with open(path, encoding="utf-8") as fh:
                for i, line in enumerate(fh):
                    if i == 0 or not line.strip():
                        continue          # header row
                    tags.add(line.split("\t")[0].strip())
        except OSError:
            continue
    return tags


def strip_reserved_info_tags(vcf_in, vcf_out, reserved, dry_run=False):
    """
    Drop INFO tags that collide with PCGR's own annotation namespace.

    WHY THIS IS NEEDED
      Step 9 annotates with 'bcftools annotate -c ID,INFO', which copies
      every INFO field COSMIC carries -- including STRAND, which PCGR also
      writes. The collision is not COSMIC's fault or PCGR's; it is the
      inevitable result of merging two annotation namespaces. Nothing is
      lost by dropping them here: PCGR re-derives all of this from VEP,
      and this rewritten copy is PCGR's input only, never the pipeline's
      output of record.

    Implemented in plain Python, matching lift_format_to_info() above, so
    the wrapper keeps working wherever bcftools is not on PATH.

    ARGS:
        vcf_in:   Input VCF (plain or gzipped).
        vcf_out:  Output VCF (plain text).
        reserved: Set of INFO tag names to remove.
        dry_run:  If True, report what would happen and write nothing.

    RETURNS:
        list: The tag names actually removed, or None if nothing was done
              (no reserved tags given, or the file could not be read).
    """
    if not reserved:
        return None

    def _keep(entry):
        """True if this INFO entry survives. Handles both KEY=V and KEY."""
        return entry.split("=", 1)[0] not in reserved

    removed = set()
    try:
        with _open_maybe_gzip(vcf_in) as src:
            lines = src.readlines()
    except OSError as exc:
        print(f"[WARN] Could not read {vcf_in} to strip INFO tags: {exc}")
        return None

    out_lines = []
    for line in lines:
        if line.startswith("##INFO=<ID="):
            tag = line.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
            if tag in reserved:
                removed.add(tag)
                continue          # drop the header declaration too
            out_lines.append(line)
        elif line.startswith("#"):
            out_lines.append(line)
        else:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 7:
                kept = [e for e in fields[7].split(";") if _keep(e)]
                # An empty INFO column is spelled ".", never blank.
                fields[7] = ";".join(kept) if kept else "."
                out_lines.append("\t".join(fields) + "\n")
            else:
                out_lines.append(line)

    if not removed:
        return []
    if dry_run:
        print(f"[DRY RUN] Would strip INFO tag(s) {', '.join(sorted(removed))}"
              f" -> {vcf_out}")
        return sorted(removed)

    try:
        with open(vcf_out, "w", encoding="utf-8") as dst:
            dst.writelines(out_lines)
    except OSError as exc:
        print(f"[WARN] Could not write {vcf_out}: {exc}")
        return None
    print(f"[INFO] Stripped INFO tag(s) reserved by PCGR: "
          f"{', '.join(sorted(removed))}")
    return sorted(removed)


def build_pcgr_command(input_vcf, output_dir, sample_id, refdata_dir,
                       vep_dir=None, genome_assembly="grch38",
                       assay="TARGETED", tumour_site=None, tumour_only=False,
                       effective_target_size_mb=None, estimate_tmb=False,
                       estimate_msi=False, estimate_signatures=False,
                       tumour_dp_tag=None, tumour_af_tag=None,
                       control_dp_tag=None, control_af_tag=None,
                       force_overwrite=True, legacy_v1=False,
                       extra_args=None):
    """
    Assemble the PCGR command line.

    Split out from run_pcgr() so it can be unit-tested and so --dry-run can
    show the exact command without any of the surrounding machinery.

    ARGS:
        input_vcf:      Somatic VCF (bgzipped is fine).
        output_dir:     Directory PCGR writes the report into.
        sample_id:      Used for the report title and output filenames.
        refdata_dir:    PCGR reference data bundle.
        vep_dir:        Ensembl VEP cache directory.
        genome_assembly: 'grch38' or 'grch37'. Must match the reference the
                        variants were called against, or coordinates are
                        meaningless.
        assay:          WGS / WES / TARGETED -- drives TMB scaling.
        effective_target_size_mb:
                        Size of the region actually sequenced, in Mb. TMB
                        is mutations divided by this, so leaving it at
                        PCGR's default while running a smaller panel
                        understates TMB in direct proportion.
        tumour_site:    PCGR's integer tumour-site code (0 = unspecified).
                        Site-specific actionability depends on it.
        tumour_only:    Set when there is no matched normal, so PCGR applies
                        its germline-filtering heuristics.
        *_dp_tag/*_af_tag: INFO tag names carrying depth / allele fraction.
                        See the module docstring -- Mutect2 does not provide
                        these in INFO by default.
        legacy_v1:      Emit the PCGR 1.x reference-bundle flag spelling.
        extra_args:     List of raw arguments appended verbatim, for options
                        this wrapper does not model.

    RETURNS:
        list: argv for subprocess.
    """
    cmd = ["pcgr", "--input_vcf", input_vcf, "--output_dir", output_dir,
           "--sample_id", sample_id]

    # The reference bundle flag was renamed between the 1.x and 2.x series.
    cmd += ["--pcgr_dir" if legacy_v1 else "--refdata_dir", refdata_dir]

    if vep_dir:
        cmd += ["--vep_dir", vep_dir]

    cmd += ["--genome_assembly", genome_assembly, "--assay", assay]

    # THESE ARE OFF BY DEFAULT IN PCGR AND MUST BE ASKED FOR.
    # PCGR 2.x computes no TMB, no MSI call and no signature fit unless
    # the corresponding --estimate_* flag is passed. It does not warn; the
    # report simply comes back without those sections, which reads as
    # "nothing to report" rather than "never calculated".
    if estimate_tmb:
        cmd.append("--estimate_tmb")
    if estimate_msi:
        cmd.append("--estimate_msi")
    if estimate_signatures:
        cmd.append("--estimate_signatures")

    # TMB is a rate: mutations per megabase of sequenced territory. PCGR
    # assumes 34 Mb for TARGETED, which is an exome-sized guess -- on a
    # 4.5 Mb panel that divides by roughly 7.5x too much and would report
    # a tumour as far less mutated than it is. Only meaningful alongside
    # --estimate_tmb; there is no safe default, only the number for the
    # assay in hand.
    if effective_target_size_mb:
        cmd += ["--effective_target_size_mb", str(effective_target_size_mb)]

    if tumour_site is not None:
        cmd += ["--tumor_site", str(tumour_site)]

    # Tumour-only mode makes PCGR apply germline-subtraction heuristics
    # (gnomAD frequency filters, panel-of-normals style reasoning). Without
    # a matched normal these are the only thing standing between you and a
    # report full of inherited variants.
    if tumour_only:
        cmd.append("--tumor_only")

    # Depth / allele-fraction INFO tags. Passed through verbatim; this
    # module does not invent or synthesise them.
    for flag, value in (("--tumor_dp_tag", tumour_dp_tag),
                        ("--tumor_af_tag", tumour_af_tag),
                        ("--control_dp_tag", control_dp_tag),
                        ("--control_af_tag", control_af_tag)):
        if value:
            cmd += [flag, value]

    if force_overwrite:
        cmd.append("--force_overwrite")

    if extra_args:
        cmd += list(extra_args)

    return cmd


def find_reports(output_dir, sample_id):
    """
    Locate whatever PCGR actually produced.

    Report filenames have shifted between PCGR releases, so glob for them
    rather than reconstructing an exact name and reporting a false negative
    when the convention changes.

    RETURNS:
        dict: {"html": [...], "tsv": [...], "vcf": [...]} of absolute paths.
    """
    found = {"html": [], "tsv": [], "vcf": []}
    if not os.path.isdir(output_dir):
        return found

    patterns = {
        "html": [f"{sample_id}*.html", "*.html"],
        "tsv": [f"{sample_id}*.tsv*", "*.tsv*"],
        "vcf": [f"{sample_id}*.vcf.gz", "*.vcf.gz"],
    }
    for kind, globs in patterns.items():
        for pattern in globs:
            hits = sorted(glob.glob(os.path.join(output_dir, pattern)))
            if hits:
                found[kind] = [os.path.abspath(h) for h in hits]
                break  # First pattern that matches wins (sample-specific).
    return found


# =============================================================================
# SECTION 3: THE STEP ITSELF
# =============================================================================

def run_pcgr(input_vcf, output_dir, sample_id, refdata_dir, vep_dir=None,
             genome_assembly="grch38", assay="TARGETED", tumour_site=None,
             effective_target_size_mb=None, estimate_tmb=False,
             estimate_msi=False, estimate_signatures=False,
             tumour_only=False, tumour_dp_tag=None, tumour_af_tag=None,
             control_dp_tag=None, control_af_tag=None, legacy_v1=False,
             extra_args=None, log_path=None, dry_run=False):
    """
    Run PCGR over a somatic VCF.

    RETURN CONVENTION (matches run_snpeff() / annotate_cosmic() in
    comprehensive_variant_calling.py, so callers can treat all three the
    same way):
        None  -> skipped: PCGR or its inputs are unavailable. Not an error.
        False -> ran and failed.
        True  -> succeeded.
    """
    # --- Preconditions. Each of these is a "skip", not a failure: PCGR is
    # --- an optional enrichment step and a missing bundle should not sink
    # --- an otherwise good variant-calling run.
    if shutil.which("pcgr") is None:
        print("[WARN] 'pcgr' not found on PATH; skipping PCGR report.")
        print("       PCGR lives in its own conda environment and is not on")
        print("       bioconda. If it is installed, activate that env first:")
        print("           conda activate pcgr && python pcgr_report.py ...")
        print("       Putting pcgr's bin/ on PATH is NOT sufficient.")
        return None

    if not refdata_dir:
        print("[WARN] No PCGR reference bundle given (--pcgr-refdata-dir); "
              "skipping PCGR report.")
        return None

    if not dry_run and not os.path.isdir(refdata_dir):
        print(f"[WARN] PCGR reference bundle not found: {refdata_dir}; "
              f"skipping PCGR report.")
        return None

    if not dry_run and not os.path.exists(input_vcf):
        print(f"[WARN] Input VCF not found: {input_vcf}; "
              f"skipping PCGR report.")
        return None

    os.makedirs(output_dir, exist_ok=True)

    # PCGR aborts if the input already carries an INFO tag it means to
    # write, and step 9's 'bcftools annotate -c ID,INFO' copies COSMIC's
    # INFO wholesale -- STRAND being the one that actually collides today.
    # Strip them into a separate file rather than editing the pipeline's
    # own output; if nothing collides, the original is used unchanged.
    reserved = pcgr_reserved_info_tags(refdata_dir, genome_assembly)
    stripped_vcf = os.path.join(output_dir, f"{sample_id}.pcgr_ready.vcf")
    removed = strip_reserved_info_tags(input_vcf, stripped_vcf, reserved,
                                       dry_run=dry_run)
    if removed:
        input_vcf = stripped_vcf

    cmd = build_pcgr_command(
        input_vcf=input_vcf,
        output_dir=output_dir,
        sample_id=sample_id,
        refdata_dir=refdata_dir,
        vep_dir=vep_dir,
        genome_assembly=genome_assembly,
        assay=assay,
        tumour_site=tumour_site,
        effective_target_size_mb=effective_target_size_mb,
        estimate_tmb=estimate_tmb,
        estimate_msi=estimate_msi,
        estimate_signatures=estimate_signatures,
        tumour_only=tumour_only,
        tumour_dp_tag=tumour_dp_tag,
        tumour_af_tag=tumour_af_tag,
        control_dp_tag=control_dp_tag,
        control_af_tag=control_af_tag,
        legacy_v1=legacy_v1,
        extra_args=extra_args,
    )

    # Loud warning rather than a silent bad number: without depth/AF tags
    # PCGR's TMB is computed over unfiltered calls.
    if not (tumour_dp_tag or tumour_af_tag):
        print("[WARN] No --tumor-dp-tag/--tumor-af-tag supplied. PCGR cannot "
              "apply depth/AF filters, so TMB in the report will be computed "
              "over unfiltered calls and will read high. Treat TMB and MSI "
              "as indicative only. See pcgr_report.py's module docstring.")

    code = run_command(cmd, tag=f"PCGR [{sample_id}]",
                       log_path=log_path, dry_run=dry_run)
    if code != 0:
        print(f"[WARN] PCGR failed (exit {code}).")
        return False

    if dry_run:
        return True

    reports = find_reports(output_dir, sample_id)
    if not any(reports.values()):
        # PCGR exited 0 but produced nothing recognisable -- report it
        # rather than claiming success.
        print("[WARN] PCGR exited cleanly but no report files were found in "
              f"{output_dir}.")
        return False

    if reports["html"]:
        print(f"[INFO] PCGR report: {reports['html'][0]}")
    return True


# =============================================================================
# SECTION 4: STANDALONE ENTRY POINT
# =============================================================================

def build_parser():
    """Build the command-line parser for standalone use."""
    parser = argparse.ArgumentParser(
        description="Run PCGR (Personal Cancer Genome Reporter) over a "
                    "somatic VCF to produce a clinical interpretation "
                    "report (actionability, TMB, MSI, mutational "
                    "signatures).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io = parser.add_argument_group("input / output")
    io.add_argument("--input-vcf", required=True,
                    help="Somatic VCF to report on (the FilterMutectCalls "
                         "or COSMIC-annotated output).")
    io.add_argument("--output-dir", required=True,
                    help="Directory for the PCGR report.")
    io.add_argument("--sample-id", required=True,
                    help="Sample identifier used in the report.")
    io.add_argument("--log", default=None,
                    help="Append command output to this log file.")

    res = parser.add_argument_group("PCGR resources")
    res.add_argument("--pcgr-refdata-dir", required=True,
                     help="PCGR reference data bundle directory (downloaded "
                          "separately; must match the installed PCGR "
                          "version and genome build).")
    res.add_argument("--vep-dir", default=None,
                     help="Ensembl VEP cache directory.")

    opt = parser.add_argument_group("report options")
    opt.add_argument("--genome-assembly", default="grch38",
                     choices=list(GENOME_CHOICES),
                     help="Genome build the variants were called against.")
    opt.add_argument("--assay", default="TARGETED",
                     choices=list(ASSAY_CHOICES),
                     help="Assay type; drives TMB scaling.")
    opt.add_argument("--effective-target-size-mb", type=float, default=None,
                     help="Megabases actually sequenced, used as the TMB "
                          "denominator. SET THIS FOR ANY PANEL: PCGR "
                          "defaults to 34 Mb for --assay TARGETED, which is "
                          "exome-sized, so a 4.5 Mb panel reports a TMB "
                          "roughly 7.5x too low. Use the panel's target "
                          "footprint, not the number of bases covered.")
    opt.add_argument("--estimate-tmb", action="store_true",
                     help="Compute tumour mutational burden. OFF in PCGR by "
                          "default, and its absence is silent: the report "
                          "simply has no TMB section. Pair with "
                          "--effective-target-size-mb on a panel.")
    opt.add_argument("--estimate-msi", action="store_true",
                     help="Predict microsatellite instability. OFF in PCGR "
                          "by default, and PCGR accepts it ONLY for WGS/WES "
                          "tumour-control runs -- on a TARGETED or "
                          "tumour-only query it logs a warning and omits the "
                          "analysis, so passing this on a panel is a no-op.")
    opt.add_argument("--estimate-signatures", action="store_true",
                     help="Fit COSMIC mutational signatures (SBS3/HRD, MMR, "
                          "APOBEC). OFF in PCGR by default. Wants a few "
                          "hundred SNVs at minimum; a small panel will not "
                          "give a trustworthy fit.")
    opt.add_argument("--tumour-site", type=int, default=None,
                     help="PCGR tumour-site code (0 = unspecified). "
                          "Site-specific actionability depends on it; see "
                          "the PCGR documentation for the code list.")
    opt.add_argument("--tumour-only", action="store_true",
                     help="No matched normal: let PCGR apply its "
                          "germline-filtering heuristics.")
    opt.add_argument("--tumor-dp-tag", default=None,
                     help="INFO tag holding tumour depth. NOTE: Mutect2 "
                          "writes DP/AD/AF as FORMAT fields, not INFO, so "
                          "this is usually empty unless you transformed the "
                          "VCF first. Without it TMB is unfiltered.")
    opt.add_argument("--tumor-af-tag", default=None,
                     help="INFO tag holding tumour allele fraction. Same "
                          "caveat as --tumor-dp-tag.")
    opt.add_argument("--control-dp-tag", default=None,
                     help="INFO tag holding control/normal depth.")
    opt.add_argument("--control-af-tag", default=None,
                     help="INFO tag holding control/normal allele fraction.")

    lift = parser.add_argument_group("FORMAT -> INFO transform")
    lift.add_argument("--lift-tags", action="store_true",
                      help="Before running PCGR, copy the tumour's (and "
                           "normal's) FORMAT DP/AF into INFO tags and point "
                           "PCGR at them. Needed for meaningful TMB with a "
                           "Mutect2 VCF, which carries those as FORMAT "
                           "fields only. Writes a new VCF; the input is left "
                           "untouched.")
    lift.add_argument("--lift-tumour-sample", default=None,
                      help="Tumour sample name, if the VCF has no "
                           "##tumor_sample header line to read it from.")
    lift.add_argument("--lift-normal-sample", default=None,
                      help="Normal sample name; likewise.")

    rt = parser.add_argument_group("runtime")
    rt.add_argument("--pcgr-legacy-v1", action="store_true",
                    help="Emit the PCGR 1.x reference-bundle flag "
                         "(--pcgr_dir) instead of the 2.x --refdata_dir.")
    rt.add_argument("--pcgr-extra-args", default=None,
                    help="Extra arguments forwarded verbatim to pcgr, as "
                         "one quoted string. Use this when a PCGR release "
                         "moves an option this wrapper does not model.")
    rt.add_argument("--dry-run", action="store_true",
                    help="Show the command without executing it.")
    return parser


def main():
    """Standalone entry point."""
    parser = build_parser()
    args = parser.parse_args()

    extra = shlex.split(args.pcgr_extra_args) if args.pcgr_extra_args else None

    input_vcf = os.path.abspath(args.input_vcf)
    dp_tag, af_tag = args.tumor_dp_tag, args.tumor_af_tag
    control_dp, control_af = args.control_dp_tag, args.control_af_tag

    # Optional FORMAT -> INFO transform. On success PCGR is pointed at the
    # rewritten VCF and the tag names are filled in automatically, so the
    # user does not have to repeat them.
    if args.lift_tags:
        lifted = os.path.join(os.path.abspath(args.output_dir),
                              f"{args.sample_id}.pcgr_input.vcf")
        stats = lift_format_to_info(
            input_vcf, lifted,
            tumour_sample=args.lift_tumour_sample,
            normal_sample=args.lift_normal_sample,
            dry_run=args.dry_run,
        )
        if stats is None:
            # Transform refused (ambiguous samples, already-lifted input,
            # unreadable file). Carry on with the original VCF rather than
            # aborting -- PCGR still produces a report, just an unfiltered
            # TMB, and run_pcgr() warns about exactly that.
            print("[WARN] FORMAT->INFO transform did not run; continuing "
                  "with the original VCF and no depth/AF tags.")
        else:
            input_vcf = lifted
            dp_tag = dp_tag or "TDP"
            af_tag = af_tag or "TAF"
            if stats.get("normal_sample"):
                control_dp = control_dp or "NDP"
                control_af = control_af or "NAF"

    ok = run_pcgr(
        input_vcf=input_vcf,
        output_dir=os.path.abspath(args.output_dir),
        sample_id=args.sample_id,
        refdata_dir=args.pcgr_refdata_dir,
        vep_dir=args.vep_dir,
        genome_assembly=args.genome_assembly,
        assay=args.assay,
        effective_target_size_mb=args.effective_target_size_mb,
        estimate_tmb=args.estimate_tmb,
        estimate_msi=args.estimate_msi,
        estimate_signatures=args.estimate_signatures,
        tumour_site=args.tumour_site,
        tumour_only=args.tumour_only,
        tumour_dp_tag=dp_tag,
        tumour_af_tag=af_tag,
        control_dp_tag=control_dp,
        control_af_tag=control_af,
        legacy_v1=args.pcgr_legacy_v1,
        extra_args=extra,
        log_path=args.log,
        dry_run=args.dry_run,
    )

    # Exit status mirrors the tri-state return: 0 ok, 2 skipped, 1 failed.
    # "Skipped" is distinguishable so a wrapper can tell "PCGR is not
    # installed" apart from "PCGR broke".
    if ok is None:
        print("[SKIP] PCGR did not run.")
        return 2
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
