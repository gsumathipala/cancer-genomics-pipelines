#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
fusion_calling.py
=================
The RNA engine: QC, splice-aware alignment, fusion detection, and the two
reports that say what the result means.

    STEP 1  QC and trimming (fastp)          -- skipped with --manifest
    STEP 2  Splice-aware alignment (STAR)    -- chimeric detection ON
    STEP 3  Coordinate-sorted BAM for QC and IGV
    STEP 4  Fusion calling (Arriba)
    STEP 5  RNA library QC                   -- is a negative interpretable?
    STEP 6  Fusion report                    -- what the calls mean

WHY THIS IS A SEPARATE ENGINE
-----------------------------
  comprehensive_variant_calling.py is fourteen steps and every one of them
  from step 3 onward is DNA: duplicate marking, base recalibration,
  Mutect2, contamination, read-orientation modelling, MSI. None of those
  has an RNA equivalent -- they have RNA CONTRADICTIONS. Duplicate marking
  on RNA removes real signal, BQSR has nothing to recalibrate against, and
  there is no validated somatic caller for RNA at all.

  So the branch forks after stage 1. Both engines consume the same QC
  manifest, both resolve the same shared genome, and the orchestrator picks
  between them with --assay. Nothing is duplicated except the four lines
  that read a manifest.

WHAT THIS ENGINE DOES NOT DO
----------------------------
  No expression quantification, no RNA variant calling, no allele-specific
  expression, no isoform work. RNA_SCOPE.md says why each was left out, and
  what it would take to add. The short version: fusions are where the
  clinical value is, and everything else on that list needs either a matched
  cohort or a validation this bundle does not have.

A NEGATIVE FUSION RESULT IS THE DANGEROUS OUTPUT
------------------------------------------------
  A degraded FFPE library produces a clean, well-formed, EMPTY fusion
  table, which is indistinguishable from a specimen that genuinely carries
  no fusion. That is why step 5 is not optional in the way the DNA
  coverage check is: it runs whenever there is a STAR log to read, and its
  verdict travels in the run manifest so that no report can quietly present
  an uninterpretable library as a negative.

Usage
-----
    python fusion_calling.py \\
        --panel illumina-tso500-rna \\
        -i fastqs/ -o rna_results/ \\
        --reference ~/data/references/hg38/hg38.fa \\
        --gtf ~/data/references/gencode/gencode.v44.annotation.gtf \\
        --star-index ~/data/references/star_hg38_150 \\
        --sample TUMOUR_RNA --threads 16
"""

import argparse
import glob
import json
import os
import shlex
import shutil
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Alignment is not optional and there is no fallback for it, so this is a
# plain import rather than the tolerant lazy one used for the reports.
# install_pipeline.py --check runs every helper's --help for this reason.
from align_rna import (star_align, sort_and_index, parse_star_log,   # noqa: E402
                       index_is_present, check_index_provenance,
                       read_index_record, load_qc_manifest,
                       default_star_index, DEFAULT_REFERENCE_DIR,
                       run_command, tool_version)

TOTAL_STEPS = 6


def banner(step, title):
    """
    The step banner, in the DNA engine's exact format.

    Not cosmetic: webapp/jobs.py derives its progress bar by matching
    "# STEP n/m: title" in the log. Emitting the same shape means the RNA
    branch gets a working progress bar without a second parser that could
    drift from the first.
    """
    print(f"\n{'#' * 72}")
    print(f"# STEP {step}/{TOTAL_STEPS}: {title}")
    print(f"{'#' * 72}")


# =============================================================================
# SECTION 1: ARRIBA AND ITS REFERENCE FILES
# =============================================================================
# Arriba ships its own small reference files -- a blacklist of recurrent
# artefacts, a list of known fusions, and protein domain annotations. They
# live inside the conda package rather than being downloaded, which is why
# the installer has no arriba-data step: installing the package installs
# them. This locates them, because the path differs between conda layouts
# and because a run that silently omits the BLACKLIST reports every
# read-through artefact in the transcriptome as a fusion.

ARRIBA_FILE_PATTERNS = {
    "blacklist": ("blacklist*hg38*.tsv.gz", "blacklist*GRCh38*.tsv.gz",
                  "blacklist*.tsv.gz"),
    "known_fusions": ("known_fusions*hg38*.tsv.gz",
                      "known_fusions*GRCh38*.tsv.gz",
                      "known_fusions*.tsv.gz"),
    "protein_domains": ("protein_domains*hg38*.gff3",
                        "protein_domains*GRCh38*.gff3",
                        "protein_domains*.gff3"),
}


def arriba_resource_dirs(explicit=None):
    """Where Arriba's database files might be, best guess first."""
    candidates = []
    if explicit:
        candidates.append(explicit)
    env = os.environ.get("ARRIBA_FILES")
    if env:
        candidates.append(env)
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        candidates.append(os.path.join(prefix, "var", "lib", "arriba"))
    executable = shutil.which("arriba")
    if executable:
        root = os.path.dirname(os.path.dirname(os.path.abspath(executable)))
        candidates.append(os.path.join(root, "var", "lib", "arriba"))
    return [d for d in candidates if d and os.path.isdir(d)]


def find_arriba_resources(explicit=None):
    """
    Locate the blacklist, known-fusion and protein-domain files.

    Returns (found, notes). Missing files are NOT fatal -- Arriba runs
    without them -- but each one is named, because the consequences differ
    and are all invisible in the output:

      blacklist        without it, every recurrent read-through artefact in
                       the transcriptome is reported as a fusion, and the
                       table fills with high-confidence noise.
      known_fusions    without it, a known recurrent fusion loses the
                       sensitivity boost Arriba gives it, so a real
                       low-support EML4-ALK can drop below the bar.
      protein_domains  cosmetic; it annotates retained domains.
    """
    found = {}
    notes = []
    directories = arriba_resource_dirs(explicit)
    if not directories:
        notes.append(
            "Arriba's reference files could not be located. Without the "
            "BLACKLIST in particular, recurrent read-through artefacts are "
            "reported as fusions and the table fills with high-confidence "
            "noise. They ship inside the conda package, usually at "
            "$CONDA_PREFIX/var/lib/arriba -- point at them with "
            "--arriba-resources.")
        return found, notes

    for key, patterns in ARRIBA_FILE_PATTERNS.items():
        for directory in directories:
            hits = []
            for pattern in patterns:
                hits = sorted(glob.glob(os.path.join(directory, pattern)))
                if hits:
                    break
            if hits:
                found[key] = hits[0]
                break
        if key not in found:
            notes.append(
                f"Arriba's {key.replace('_', ' ')} file was not found in "
                f"{', '.join(directories)}; the run continues without it.")
    return found, notes


def run_arriba(bam, reference, gtf, out_tsv, discarded_tsv, resources,
               log_path=None, dry_run=False, extra_args=None):
    """
    Call fusions from STAR's chimeric-annotated BAM.

    Arriba takes the UNSORTED BAM on purpose: it needs mates adjacent, and
    coordinate sorting separates them. Handing it the sorted BAM is a
    supported command that finds far fewer fusions, which is exactly the
    kind of quiet wrongness this bundle tries to make impossible -- so the
    caller of this function passes star_align()'s own output and the sorted
    copy is produced separately, for QC and IGV.
    """
    cmd = ["arriba", "-x", bam, "-o", out_tsv, "-O", discarded_tsv,
           "-a", reference, "-g", gtf]
    if resources.get("blacklist"):
        cmd += ["-b", resources["blacklist"]]
    if resources.get("known_fusions"):
        # -k is the known-fusion list (a sensitivity boost for recurrent
        # events); -t tags calls that appear in it. Arriba's own
        # documentation passes the same file to both.
        cmd += ["-k", resources["known_fusions"],
                "-t", resources["known_fusions"]]
    if resources.get("protein_domains"):
        cmd += ["-p", resources["protein_domains"]]
    if extra_args:
        cmd += list(extra_args)

    code = run_command(cmd, tag="Arriba", log_path=log_path, dry_run=dry_run)
    if code != 0:
        print("[ERROR] Arriba failed.")
        return False
    return True


# =============================================================================
# SECTION 2: QC (fastp)
# =============================================================================

def run_fastp(sample, r1, r2, out_r1, out_r2, reports_dir, args,
              log_path=None, dry_run=False):
    """
    Trim one sample, with the RNA-appropriate settings.

    Deliberately thin compared with the DNA engine's fastp step. Two
    differences matter:

      * the minimum read length floor is lower (the panel profiles set 35
        rather than 50), because FFPE RNA is fragmented and a floor
        inherited from DNA discards a large part of a usable library;
      * no overlap correction, for the same reason it is off on the DNA
        side -- it rewrites the bases where mates disagree, and on a fusion
        breakpoint those are the bases carrying the evidence.
    """
    os.makedirs(reports_dir, exist_ok=True)
    cmd = ["fastp", "-i", r1, "-o", out_r1]
    if r2:
        cmd += ["-I", r2, "-O", out_r2]
    cmd += [
        "--thread", str(min(args.threads, 16)),
        "--length_required", str(args.min_read_length),
        "--qualified_quality_phred", str(args.min_base_quality),
        "--n_base_limit", str(args.max_n_bases),
        "--html", os.path.join(reports_dir, f"{sample}.fastp.html"),
        "--json", os.path.join(reports_dir, f"{sample}.fastp.json"),
        "--report_title", f"fastp (RNA): {sample}",
    ]
    if r2 and not args.no_detect_adapter_pe:
        cmd += ["--detect_adapter_for_pe"]
    if args.adapter_r1:
        cmd += ["--adapter_sequence", args.adapter_r1]
    if args.adapter_r2 and r2:
        cmd += ["--adapter_sequence_r2", args.adapter_r2]
    if not args.no_poly_g:
        cmd += ["--trim_poly_g"]
    if args.umi_loc:
        cmd += ["--umi", "--umi_loc", args.umi_loc]
        if args.umi_len:
            cmd += ["--umi_len", str(args.umi_len)]
        if args.umi_skip is not None:
            cmd += ["--umi_skip", str(args.umi_skip)]
    if args.reads_to_process:
        cmd += ["--reads_to_process", str(args.reads_to_process)]

    return run_command(cmd, tag=f"fastp [{sample}]", log_path=log_path,
                       dry_run=dry_run) == 0


# =============================================================================
# SECTION 3: DISCOVERY
# =============================================================================

def find_fastqs(input_dir, r1_pattern, r2_pattern, recursive=False):
    """
    FASTQ pairs in a directory, reusing the DNA engine's own pairing rules.

    Imported rather than reimplemented so the two branches cannot disagree
    about what constitutes a pair -- the same reasoning the webapp uses
    when it resolves samples itself.
    """
    try:
        from comprehensive_variant_calling import find_paired_fastqs
    except ImportError:
        print("[WARN] comprehensive_variant_calling.py not importable; "
              "cannot auto-discover FASTQ pairs.")
        return [], []
    return find_paired_fastqs(input_dir, r1_pattern, r2_pattern,
                              recursive=recursive)


# =============================================================================
# SECTION 4: ARGUMENTS
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="RNA fusion detection for cancer panels: QC, STAR, "
                    "Arriba, and the reports that say what the result "
                    "means.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="See RNA.md for the whole branch and RNA_SCOPE.md for what "
               "it deliberately does not do.")

    io = parser.add_argument_group("input / output")
    io.add_argument("-i", "--input-dir", default=None,
                    help="Directory of FASTQs. Not needed with --manifest.")
    io.add_argument("-o", "--output-dir", required=True,
                    help="Directory for every output of this run.")
    io.add_argument("--manifest", default=None,
                    help="fastq_qc_clean.py run manifest. Reuses that run's "
                         "cleaned reads and skips step 1.")
    io.add_argument("--sample", default=None,
                    help="Sample name for the explicit-input mode.")
    io.add_argument("--r1", default=None, help="Read 1 FASTQ.")
    io.add_argument("--r2", default=None,
                    help="Read 2 FASTQ. Omit for a single-end library.")
    io.add_argument("--auto-discover", action="store_true",
                    help="Find FASTQ pairs in --input-dir. Unlike the DNA "
                         "engine there is no tumour/normal pairing here, so "
                         "EVERY pair found is processed as its own sample.")
    io.add_argument("--r1-pattern", default="*_R1_*.fastq*",
                    help="Glob for read 1 files.")
    io.add_argument("--r2-pattern", default="*_R2_*.fastq*",
                    help="Glob for read 2 files.")
    io.add_argument("--recursive", action="store_true",
                    help="Search subdirectories.")

    ref = parser.add_argument_group("reference")
    ref.add_argument("-r", "--reference", default=None,
                     help="Genome FASTA, the same one the DNA branch uses. "
                          "Arriba needs it to read the sequence either side "
                          "of a breakpoint.")
    ref.add_argument("--gtf", default=None,
                     help="GENCODE annotation GTF. Required: it is how a "
                          "breakpoint becomes a gene name.")
    ref.add_argument("--star-index", default=None,
                     help="STAR index directory. Built once by "
                          "align_rna.py --build-index and shared.")
    ref.add_argument("--reference-dir", default=DEFAULT_REFERENCE_DIR,
                     help="Where shared references live, used to guess "
                          "--star-index.")
    ref.add_argument("--read-length", type=int, default=None,
                     help="Library read length, checked against the index's "
                          "build record.")
    ref.add_argument("--arriba-resources", default=None, metavar="DIR",
                     help="Directory holding Arriba's blacklist, known-"
                          "fusion and protein-domain files. Found "
                          "automatically inside the conda environment; "
                          "supply it only when that fails. WITHOUT THE "
                          "BLACKLIST the table fills with read-through "
                          "artefacts reported as high-confidence fusions.")

    fusion = parser.add_argument_group("fusion calling")
    fusion.add_argument("--fusion-caller", default="arriba",
                        choices=["arriba"],
                        help="Only Arriba is implemented. STAR-Fusion is "
                             "supported by the profile schema but needs the "
                             "~40 GB CTAT library, which this bundle does "
                             "not install -- see RNA_SCOPE.md.")
    fusion.add_argument("--fusion-min-confidence", default="medium",
                        choices=["low", "medium", "high"],
                        help="The bar at which a call counts as reportable "
                             "in the report. Calls below it are MARKED, "
                             "never dropped.")
    fusion.add_argument("--min-fusion-reads", type=int, default=3,
                        help="Advisory floor recorded with the run; the "
                             "report shows the support for every call.")
    fusion.add_argument("--strandedness", default="unstranded",
                        choices=["unstranded", "forward", "reverse"],
                        help="Library strandedness, recorded for the "
                             "report. Arriba does not need it; anything "
                             "that reports direction does.")
    fusion.add_argument("--arriba-extra-args", default=None,
                        help="Extra arguments forwarded verbatim to Arriba.")
    fusion.add_argument("--star-extra-args", default=None,
                        help="Extra arguments forwarded verbatim to STAR, "
                             "appended after the chimeric-detection block "
                             "so they override it.")

    qc = parser.add_argument_group("QC / trimming")
    qc.add_argument("--min-read-length", type=int, default=35,
                    help="Minimum read length after trimming. Lower than "
                         "the DNA default of 50 on purpose: FFPE RNA is "
                         "fragmented and a 50 bp floor discards much of a "
                         "usable library.")
    qc.add_argument("--min-base-quality", type=int, default=15,
                    help="Minimum base Phred quality for qualified bases.")
    qc.add_argument("--max-n-bases", type=int, default=5,
                    help="Maximum N bases per read.")
    qc.add_argument("--adapter-r1", default=None, help="Read 1 adapter.")
    qc.add_argument("--adapter-r2", default=None, help="Read 2 adapter.")
    qc.add_argument("--no-detect-adapter-pe", action="store_true",
                    help="Disable PE adapter auto-detection.")
    qc.add_argument("--no-poly-g", action="store_true",
                    help="Disable poly-G trimming.")
    qc.add_argument("--umi-loc", default=None,
                    choices=["index1", "index2", "read1", "read2",
                             "per_index", "per_read"],
                    help="UMI location, for a barcoded RNA kit.")
    qc.add_argument("--umi-len", type=int, default=None, help="UMI length.")
    qc.add_argument("--umi-skip", type=int, default=None,
                    help="Bases to skip after the UMI.")
    qc.add_argument("--rna-min-reads-millions", type=float, default=None,
                    help="Library-size floor for the QC verdict. Comes from "
                         "the panel profile in a normal run.")
    qc.add_argument("--rna-min-unique-mapped-pct", type=float, default=None,
                    help="Unique-mapping floor for the QC verdict.")

    panel = parser.add_argument_group("panel / assay profile")
    panel.add_argument("--panel", default=None, metavar="ID",
                       help="RNA panel profile configuring this run. Run "
                            "--list-panels; RNA profiles are marked [RNA].")
    panel.add_argument("--panel-file", action="append", default=[],
                       metavar="JSON", help="Additional profile file.")
    panel.add_argument("--list-panels", action="store_true",
                       help="List every known panel profile and exit.")
    panel.add_argument("--describe-panel", default=None, metavar="ID",
                       help="Print everything a profile sets, and exit.")

    run = parser.add_argument_group("runtime")
    run.add_argument("--threads", type=int, default=8, help="CPU threads.")
    run.add_argument("--reads-to-process", type=int, default=None,
                     help="Only process this many reads, for a quick test.")
    run.add_argument("--skip-steps", nargs="*", default=[],
                     help="Skip stages: qc, align, sort, fusions, rnaqc, "
                          "report.")
    run.add_argument("--dry-run", action="store_true",
                     help="Print every command without running any.")
    return parser


def apply_panel(args, parser, argv):
    """Configure the run from an RNA panel profile, if one was named."""
    record = {"panel": None, "applied": [], "overridden": [], "notes": []}
    if not args.panel:
        return record
    try:
        import panel_profiles
    except ImportError as exc:
        parser.error(f"--panel needs panel_profiles.py beside this script "
                     f"({exc}).")
        return record
    try:
        registry = panel_profiles.available_panels(
            extra_files=args.panel_file,
            warn=lambda msg: print(f"[WARN] {msg}"))
        profile = panel_profiles.resolve_panel(args.panel, registry)
        if panel_profiles.assay_type(profile) != "rna":
            parser.error(
                f"'{profile['id']}' is a DNA profile "
                f"({profile['chemistry']}), and this is the RNA engine. Its "
                f"settings -- duplicate marking, allele-fraction floors, "
                f"target BEDs -- have no meaning here. Use an RNA profile "
                f"(--list-panels marks them [RNA]), or run "
                f"comprehensive_variant_calling.py instead.")
        # The profile's rna_* keys map onto this script's flags, which drop
        # the prefix for the same reason the standalone tools do.
        rename = {"rna_min_reads_millions": "rna_min_reads_millions",
                  "rna_min_unique_mapped_pct": "rna_min_unique_mapped_pct"}
        applied, overridden = panel_profiles.apply_profile(
            args, profile, panel_profiles.explicit_dests(parser, argv),
            only=panel_profiles.RNA_SETTINGS, rename=rename)
    except panel_profiles.PanelError as exc:
        parser.error(str(exc))
        return record
    print()
    print(panel_profiles.format_application(profile, applied, overridden))
    record["panel"] = dict(profile)
    record["applied"] = [[k, v] for k, v in applied]
    record["overridden"] = [[k, v] for k, v in overridden]
    return record


# =============================================================================
# SECTION 5: MAIN
# =============================================================================

def main():
    try:
        import panel_profiles
        panel_profiles.handle_panel_queries()
    except ImportError:
        pass

    parser = build_parser()
    args = parser.parse_args()
    panel_record = apply_panel(args, parser, sys.argv[1:])

    started = datetime.now(timezone.utc)
    skip = set(args.skip_steps)

    if not args.gtf:
        parser.error(
            "--gtf is required: without an annotation a breakpoint is a "
            "pair of coordinates and never becomes a gene name. The "
            "installer puts one in ~/data/references/gencode.")
    if not args.reference:
        parser.error("-r/--reference is required; Arriba reads the sequence "
                     "either side of each breakpoint from it.")
    if not args.star_index:
        args.star_index = default_star_index(args.reference_dir,
                                             args.read_length)
        print(f"[INFO] --star-index not given; using {args.star_index}")

    # ---- tools ----
    if not args.dry_run:
        needed = {"STAR": "align", "arriba": "fusions", "samtools": "sort"}
        missing = [tool for tool, step in needed.items()
                   if step not in skip and shutil.which(tool) is None]
        if args.manifest is None and "qc" not in skip and \
                shutil.which("fastp") is None:
            missing.append("fastp")
        if missing:
            print(f"[ERROR] Required tools not found: {', '.join(missing)}")
            print("        They live in the cancer_rna environment: "
                  "conda activate cancer_rna")
            return 1

    output_dir = os.path.abspath(args.output_dir)
    dirs = {
        "cleaned": os.path.join(output_dir, "cleaned_fastq"),
        "align": os.path.join(output_dir, "rna_aligned"),
        "fusions": os.path.join(output_dir, "fusions"),
        "rna_qc": os.path.join(output_dir, "rna_qc"),
        "report": os.path.join(output_dir, "fusion_report"),
        "reports": os.path.join(output_dir, "fastp_reports"),
        "logs": os.path.join(output_dir, "logs"),
    }
    for directory in dirs.values():
        os.makedirs(directory, exist_ok=True)

    # ---- samples ----
    samples = {}
    if args.manifest:
        samples = load_qc_manifest(args.manifest)
        if not samples:
            print("[ERROR] no usable samples in the manifest.")
            return 1
        skip.add("qc")
        print(f"[INFO] Reusing cleaned reads from {args.manifest}; step 1 "
              f"is skipped.")
    elif args.auto_discover:
        if not args.input_dir:
            parser.error("--auto-discover needs --input-dir.")
        pairs, warnings = find_fastqs(args.input_dir, args.r1_pattern,
                                      args.r2_pattern, args.recursive)
        for warning in warnings:
            print(f"[WARN] {warning}")
        if not pairs:
            print(f"[ERROR] no FASTQ pairs found in {args.input_dir}.")
            return 1
        # Every pair is its own sample. There is no tumour/normal pairing
        # on this branch -- a fusion is called from one library -- so the
        # DNA engine's "first pair is the tumour, second is its normal"
        # rule, and the danger that comes with it, does not apply here.
        samples = {p["sample"]: {"r1": p["r1"], "r2": p.get("r2")}
                   for p in pairs}
        print(f"[INFO] {len(samples)} sample(s) discovered; each is "
              f"processed independently.")
    elif args.sample and args.r1:
        samples = {args.sample: {"r1": os.path.abspath(args.r1),
                                 "r2": os.path.abspath(args.r2)
                                 if args.r2 else None}}
    else:
        parser.error("Provide --manifest, --auto-discover, or --sample with "
                     "--r1.")

    # ---- index sanity ----
    if "align" not in skip and not args.dry_run:
        if not index_is_present(args.star_index):
            print(f"[ERROR] no STAR index in {args.star_index}. Build one "
                  f"first:")
            print(f"        python align_rna.py --build-index --reference "
                  f"{args.reference} \\")
            print(f"            --gtf {args.gtf} --read-length "
                  f"{args.read_length or '<len>'} --star-index "
                  f"{args.star_index}")
            return 1
        for warning in check_index_provenance(args.star_index, args.gtf,
                                              args.read_length):
            print(f"[WARN] {warning}")

    # ---- Arriba's reference files ----
    resources, resource_notes = ({}, [])
    if "fusions" not in skip:
        resources, resource_notes = find_arriba_resources(
            args.arriba_resources)
        for note in resource_notes:
            print(f"[WARN] {note}")
        if resources:
            print(f"[INFO] Arriba reference files: "
                  f"{', '.join(sorted(resources))}")

    panel_label = None
    if panel_record.get("panel"):
        panel_label = (f"{panel_record['panel']['name']} "
                       f"[{panel_record['panel']['id']}]")

    print(f"\n{'=' * 72}")
    print("RNA FUSION DETECTION PIPELINE")
    print(f"{'=' * 72}")
    print(f"Samples : {', '.join(sorted(samples))}")
    print(f"Panel   : {panel_label or 'none named'}")
    print(f"Index   : {args.star_index}")
    print(f"GTF     : {os.path.basename(args.gtf)}")
    print(f"Threads : {args.threads}")
    print(f"Output  : {output_dir}")
    print(f"Skip    : {', '.join(sorted(skip)) or 'none'}")
    print(f"{'=' * 72}")

    manifest = {
        "script": os.path.basename(__file__),
        "assay_type": "rna",
        "started_utc": started.isoformat(),
        "command_line": shlex.join(sys.argv),
        "parameters": vars(args),
        "panel_profile": panel_record,
        "reference": os.path.abspath(args.reference),
        "gtf": os.path.abspath(args.gtf),
        "star_index": os.path.abspath(args.star_index),
        "index_record": read_index_record(args.star_index),
        "arriba_resources": resources,
        "arriba_resource_notes": resource_notes,
        "tool_versions": {name: tool_version(name)
                          for name in ("STAR", "arriba", "samtools",
                                       "fastp")},
        "samples": {},
        "steps_completed": [],
    }

    results = {}
    for name in sorted(samples):
        reads = samples[name]
        record = {"sample": name}
        log_path = os.path.join(dirs["logs"], f"{name}.log")

        # ---- STEP 1: QC ----
        banner(1, "QC and trimming (fastp)")
        r1, r2 = reads["r1"], reads.get("r2")
        if "qc" in skip:
            print("[SKIP] QC skipped (reads are already cleaned).")
        else:
            clean_r1 = os.path.join(dirs["cleaned"],
                                    f"{name}_R1.clean.fastq.gz")
            clean_r2 = (os.path.join(dirs["cleaned"],
                                     f"{name}_R2.clean.fastq.gz")
                        if r2 else None)
            if not run_fastp(name, r1, r2, clean_r1, clean_r2,
                             dirs["reports"], args, log_path=log_path,
                             dry_run=args.dry_run):
                print(f"[FATAL] fastp failed for {name}.")
                return 1
            r1, r2 = clean_r1, clean_r2
            manifest["steps_completed"].append(f"qc:{name}")
        record["r1"], record["r2"] = r1, r2

        # ---- STEP 2: ALIGN ----
        banner(2, "Splice-aware alignment (STAR, chimeric detection on)")
        aligned = None
        if "align" in skip:
            print("[SKIP] Alignment skipped.")
            guess = os.path.join(dirs["align"], f"{name}.Aligned.out.bam")
            aligned = {"bam": guess,
                       "star_log": os.path.join(dirs["align"],
                                                f"{name}.Log.final.out")}
        else:
            aligned = star_align(
                name, r1, r2, args.star_index, dirs["align"], args.threads,
                log_path=log_path, dry_run=args.dry_run,
                extra_args=(shlex.split(args.star_extra_args)
                            if args.star_extra_args else None))
            if aligned is None:
                print(f"[FATAL] STAR failed for {name}.")
                return 1
            manifest["steps_completed"].append(f"align:{name}")
        record.update({k: v for k, v in aligned.items() if k != "sample"})

        stats = ({} if args.dry_run
                 else parse_star_log(aligned.get("star_log", "")))
        record["star_stats"] = stats
        if stats:
            print(f"[OK] {stats.get('Uniquely mapped reads %', '?')}% "
                  f"uniquely mapped, "
                  f"{stats.get('Number of chimeric reads', 0):,} chimeric "
                  f"reads.")

        # ---- STEP 3: SORT ----
        banner(3, "Coordinate-sorted BAM (for QC and IGV)")
        sorted_bam = os.path.join(dirs["align"], f"{name}.sorted.bam")
        if "sort" in skip:
            print("[SKIP] Sorting skipped.")
        elif sort_and_index(aligned["bam"], sorted_bam, args.threads,
                            log_path=log_path, dry_run=args.dry_run):
            record["sorted_bam"] = sorted_bam
            manifest["steps_completed"].append(f"sort:{name}")
            print("[INFO] The fusion caller uses the UNSORTED BAM (it needs "
                  "mates adjacent); this copy is for everything else.")
        else:
            print("[WARN] sorting failed; the QC step will have less to "
                  "measure, and the calls are unaffected.")

        # ---- STEP 4: FUSIONS ----
        banner(4, f"Fusion calling ({args.fusion_caller})")
        fusions_tsv = os.path.join(dirs["fusions"], f"{name}.fusions.tsv")
        discarded = os.path.join(dirs["fusions"],
                                 f"{name}.fusions.discarded.tsv")
        if "fusions" in skip:
            print("[SKIP] Fusion calling skipped.")
        else:
            ok = run_arriba(
                aligned["bam"], args.reference, args.gtf, fusions_tsv,
                discarded, resources, log_path=log_path,
                dry_run=args.dry_run,
                extra_args=(shlex.split(args.arriba_extra_args)
                            if args.arriba_extra_args else None))
            if not ok:
                print(f"[FATAL] fusion calling failed for {name}.")
                return 1
            record["fusions_tsv"] = fusions_tsv
            record["discarded_tsv"] = discarded
            manifest["steps_completed"].append(f"fusions:{name}")

        # ---- STEP 5: RNA LIBRARY QC ----
        banner(5, "RNA library QC (is a negative interpretable?)")
        if "rnaqc" in skip:
            print("[SKIP] RNA QC skipped. Nothing else in this run says "
                  "whether an empty fusion table is a negative result or an "
                  "unusable library.")
        else:
            try:
                from rna_qc_report import run_rna_qc
            except ImportError as exc:
                print(f"[WARN] rna_qc_report.py not importable ({exc}); no "
                      f"library QC statement will be made.")
                run_rna_qc = None
            if run_rna_qc is not None:
                qc_result = run_rna_qc(
                    name, aligned.get("star_log"),
                    record.get("sorted_bam"), args.gtf, dirs["rna_qc"],
                    min_reads_millions=args.rna_min_reads_millions,
                    min_unique_mapped_pct=args.rna_min_unique_mapped_pct,
                    panel_label=panel_label, dry_run=args.dry_run)
                record["rna_qc"] = {
                    "html": qc_result.get("html"),
                    "json": qc_result.get("json"),
                    "verdict": (qc_result.get("metrics") or {}).get(
                        "verdict"),
                }
                if qc_result["ok"] and not args.dry_run:
                    verdict = qc_result["metrics"]["verdict"]
                    print(f"[{verdict['state'].upper()}] {verdict['text']}")
                    print(f"[OK] {qc_result['html']}")
                    manifest["steps_completed"].append(f"rnaqc:{name}")

        # ---- STEP 6: FUSION REPORT ----
        banner(6, "Fusion report")
        if "report" in skip:
            print("[SKIP] Fusion report skipped.")
        elif not record.get("fusions_tsv") and not args.dry_run:
            print("[SKIP] No fusion table to report on.")
        else:
            try:
                from fusion_report import run_fusion_report
            except ImportError as exc:
                print(f"[WARN] fusion_report.py not importable ({exc}); the "
                      f"caller's raw TSV is the only output.")
                run_fusion_report = None
            if run_fusion_report is not None:
                report = run_fusion_report(
                    fusions_tsv, name, dirs["report"],
                    min_confidence=args.fusion_min_confidence,
                    panel_label=panel_label, caller=args.fusion_caller,
                    dry_run=args.dry_run)
                record["fusion_report"] = {
                    "html": report.get("html"),
                    "json": report.get("json"),
                    "tsv": report.get("tsv"),
                    # PCGR 2.x takes RNA fusions as a first-class molecular
                    # input, so this is the file that turns a fusion table
                    # into a clinical interpretation -- on its own, or
                    # alongside the DNA library's VCF for one combined
                    # report per specimen. Written every run; whether PCGR
                    # is then invoked is a separate decision.
                    "pcgr_tsv": report.get("pcgr_tsv"),
                    "summary": report.get("summary"),
                }
                if report["ok"] and not args.dry_run:
                    summary = report["summary"]
                    print(f"[OK] {report['html']}")
                    print(f"[OK] {summary['total_called']} call(s), "
                          f"{summary['involving_actionable_gene']} involving "
                          f"an actionable gene.")
                    manifest["steps_completed"].append(f"report:{name}")
                    if report.get("pcgr_tsv"):
                        print(f"[OK] PCGR input written: "
                              f"{report['pcgr_tsv']} "
                              f"({summary.get('pcgr_fusions_written', 0)} "
                              f"fusion(s) at or above "
                              f"'{args.fusion_min_confidence}'). Produce the "
                              f"clinical report with pcgr_report.py "
                              f"--input-rna-fusion, adding --input-vcf from "
                              f"this specimen's DNA run for one combined "
                              f"report.")
                    else:
                        print("[INFO] No fusions cleared the confidence bar, "
                              "so no PCGR input was written. PCGR rejects an "
                              "empty fusion file; no fusions found is a "
                              "result, not an error.")
                elif not report["ok"]:
                    print(f"[WARN] {report['error']}")

        results[name] = record

    manifest["samples"] = results
    finished = datetime.now(timezone.utc)
    manifest["finished_utc"] = finished.isoformat()
    manifest["duration_seconds"] = round(
        (finished - started).total_seconds(), 1)

    path = os.path.join(
        output_dir,
        f"rna_manifest_{started.strftime('%Y%m%dT%H%M%SZ')}.json")
    if not args.dry_run:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, default=str)

    print(f"\n{'=' * 72}")
    print("RNA PIPELINE COMPLETE")
    print(f"{'=' * 72}")
    for name, record in sorted(results.items()):
        summary = (record.get("fusion_report") or {}).get("summary") or {}
        verdict = (record.get("rna_qc") or {}).get("verdict") or {}
        print(f"{name}: {summary.get('total_called', '?')} fusion call(s), "
              f"{summary.get('involving_actionable_gene', '?')} on the "
              f"actionable list; library QC "
              f"{verdict.get('state', 'not assessed')}")
    if not args.dry_run:
        print(f"Manifest: {path}")
    print("\nA negative fusion result is only as good as the library it "
          "came from. Read the RNA QC report before reporting one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
