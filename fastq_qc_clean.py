#!/usr/bin/env python3
"""
fastq_qc_clean.py
=================
QC and clean Illumina NextSeq 1000 paired-end FASTQ files in preparation for a
somatic / cancer bioinformatics pipeline.

What this script DOES
---------------------
  1. Discovers R1/R2 pairs by filename token (default '_R1_' / '_R2_').
  2. Optionally merges lane-split files (L001..L004) per sample.
  3. Runs FastQC on raw reads (optional) and on cleaned reads (optional).
  4. Runs fastp for adapter removal, poly-G trimming, length and quality
     filtering, and writes per-sample HTML + JSON QC reports.
  5. Verifies that outputs actually exist and are non-trivial in size.
  6. Writes a run manifest (JSON) recording tool versions, exact command lines,
     input file sizes/checksums, and the key fastp metrics per sample.

What this script DELIBERATELY DOES NOT DO
-----------------------------------------
  * No sliding-window 3' quality trimming (`--cut_right`) by default. For
    low-VAF somatic calling, aggressive end trimming distorts allele fractions
    at read termini; soft-clipping plus BQSR handles residual low-quality
    bases. Enable with --cut-right only if you have a specific reason.
  * No overlap-based base correction (`--correction`) by default, because it
    interferes with UMI consensus / duplex calling.

Requirements (must be installed and on PATH)
--------------------------------------------
  - Python 3.8+
  - fastp    https://github.com/OpenGene/fastp
  - FastQC   https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  - (optional) multiqc, to aggregate the per-sample reports afterwards

DEBUGGING NOTES (things that surprise people)
---------------------------------------------
  * SAMPLE NAMES depend on --merge-lanes. Without it each lane is its own
    sample and the name KEEPS the lane token ("S1_L001"); with it the lanes
    are concatenated and the name DROPS it ("S1"). If downstream scripts
    report an unexpected sample name, check this flag first.
    See find_paired_fastqs() (no lane stripping) + group_by_lane().

  * THE MANIFEST IS THE CONTRACT. run_manifest_<timestamp>.json is what
    align_reads.py and comprehensive_variant_calling.py read via their
    --manifest option. Per-sample keys they rely on: "sample", "status"
    ("ok" | "failed" | "skipped_existing") and "outputs".{r1,r2}.
    Careful: records with status "skipped_existing" have NO "outputs" key
    -- the readers reconstruct the standard cleaned path instead. Changing
    that convention breaks both downstream scripts.

  * --dry-run runs no tools, but it DOES still create the output directory
    skeleton, and it does not merge lanes (so the fastp command it prints
    for a lane-split sample names only the first lane's file).

  * FLAG ALIASES: --min-read-length/--min-length, --min-base-quality/
    --qualified-quality, --max-unqualified-pct/--unqualified-percent and
    --max-n-bases/--n-base-limit are the same options. The first spelling
    is canonical and matches comprehensive_variant_calling.py; the second
    is kept so older command lines keep working.

  * -q (--min-base-quality) is NOT a mean-quality filter. It only defines
    which bases count as "unqualified" for --max-unqualified-pct. The
    mean-quality read filter is --average-qual.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

# A cleaned FASTQ smaller than this is almost certainly a truncated write.
MIN_PLAUSIBLE_OUTPUT_BYTES = 1024


# ---------------------------------------------------------------------------
# Command execution
# ---------------------------------------------------------------------------
def run_command(cmd, tag, log_path=None, dry_run=False):
    """
    Execute a command, streaming its output line by line to the console and
    (optionally) to a log file.

    Returns the exit code. Does NOT terminate the script: the caller decides
    whether a failure is fatal, so that one bad sample does not abort a batch.
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

        # Popen + line iteration gives real progress on long-running samples,
        # unlike subprocess.run(stdout=PIPE), which buffers until the process
        # exits.
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
        print(f"[{tag}] ERROR: '{cmd[0]}' not found. Is it installed and on PATH?")
        return 127
    except OSError as exc:
        print(f"[{tag}] ERROR: {exc}")
        return 1
    finally:
        if log_fh:
            log_fh.close()


def tool_version(executable):
    """Best-effort capture of a tool's version string, for the run manifest."""
    if shutil.which(executable) is None:
        return None
    for flag in ("--version", "-v", "-V"):
        try:
            res = subprocess.run(
                [executable, flag],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=30,
            )
            out = (res.stdout or "").strip()
            if out:
                return out.splitlines()[0].strip()
        except (OSError, subprocess.SubprocessError):
            continue
    return "unknown"


def file_digest(path, max_bytes=64 * 1024 * 1024):
    """
    SHA-256 of the first `max_bytes` of a file, for provenance.

    Hashing entire FASTQ files is slow for large panels; the leading 64 MB is
    enough to detect an accidentally swapped or re-demultiplexed input. Set
    max_bytes=None to hash the whole file.
    """
    h = hashlib.sha256()
    read = 0
    with open(path, "rb") as fh:
        while True:
            chunk_size = 1024 * 1024
            if max_bytes is not None:
                remaining = max_bytes - read
                if remaining <= 0:
                    break
                chunk_size = min(chunk_size, remaining)
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
            read += len(chunk)
    return {
        "sha256": h.hexdigest(),
        "bytes_hashed": read,
        "partial": max_bytes is not None and read >= max_bytes,
    }


# ---------------------------------------------------------------------------
# FASTQ discovery
# ---------------------------------------------------------------------------
def derive_tokens(r1_pattern, r2_pattern):
    """
    Extract the read-identifying tokens from the two glob patterns.

    Works for '*_R1_*.fastq*' / '*_R2_*.fastq*' (-> '_R1_' / '_R2_') and also
    for '*_1.fq.gz' / '*_2.fq.gz' (-> '_1.' / '_2.') style conventions: we take
    the longest common structure and report the single differing segment.
    """
    r1_bits = r1_pattern.split("*")
    r2_bits = r2_pattern.split("*")
    if len(r1_bits) != len(r2_bits):
        return None, None

    diffs = [(a, b) for a, b in zip(r1_bits, r2_bits) if a != b]
    if len(diffs) != 1:
        return None, None
    t1, t2 = diffs[0]
    if not t1 or not t2:
        return None, None
    return t1, t2


def find_paired_fastqs(input_dir, r1_pattern, r1_token, r2_token, recursive):
    """
    Locate R1/R2 pairs under `input_dir`.

    Returns (pairs, warnings) where pairs is a list of dicts:
        {"sample": str, "r1": path, "r2": path}

    Notes on the two bugs this replaces:
      * The pattern is now anchored at `input_dir` itself for non-recursive
        searches, so FASTQs sitting directly in the input directory are found.
      * glob(recursive=True) is required for '**' to actually descend.
    """
    warnings = []
    if recursive:
        pattern = os.path.join(input_dir, "**", r1_pattern)
    else:
        pattern = os.path.join(input_dir, r1_pattern)
    r1_files = sorted(glob.glob(pattern, recursive=recursive))

    pairs = []
    seen_samples = {}
    for r1 in r1_files:
        directory, base = os.path.split(r1)

        if r1_token not in base:
            warnings.append(f"'{base}' matched the pattern but has no "
                            f"'{r1_token}' token; skipping.")
            continue

        # Split on the LAST occurrence, and only within the basename, so that a
        # parent directory named e.g. 'run_R1_backup' is never rewritten.
        stem = base.rsplit(r1_token, 1)[0]
        r2_base = r2_token.join(base.rsplit(r1_token, 1))
        r2 = os.path.join(directory, r2_base)

        if not os.path.exists(r2):
            warnings.append(f"No matching R2 for {r1} (expected {r2}); skipping.")
            continue

        # Duplicate detection uses the full basename minus the FASTQ
        # extension as identity. This keeps lanes distinct whether the lane
        # token precedes (S1_L001_R1.fastq.gz) or follows
        # (SAME_R1_L001.fastq.gz) the read token, while still catching true
        # duplicates (identical filename twice).
        identity = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base)
        if identity in seen_samples:
            warnings.append(
                f"Duplicate sample name '{identity}': {r1} collides with "
                f"{seen_samples[identity]}. Skipping the duplicate to avoid "
                f"silently overwriting output."
            )
            continue
        seen_samples[identity] = r1

        pairs.append({"sample": stem, "r1": r1, "r2": r2})

    return pairs, warnings


LANE_RE = re.compile(r"_L\d{3}(?=_|$)")


def group_by_lane(pairs):
    """
    Group lane-split pairs (…_L001_…, …_L002_…) under a single sample key.

    Returns (groups, is_lane_split) where groups maps merged_sample_name to a
    list of the constituent pair dicts, in lane order.
    """
    groups = {}
    for pair in pairs:
        merged = LANE_RE.sub("", pair["sample"])
        groups.setdefault(merged, []).append(pair)
    for key in groups:
        groups[key].sort(key=lambda p: p["sample"])
    is_lane_split = any(len(v) > 1 for v in groups.values())
    return groups, is_lane_split


def concatenate(files, destination):
    """
    Concatenate gzip members into one file. Concatenated gzip streams are a
    valid gzip file, so this is safe for .fastq.gz without recompression.
    Plain (uncompressed) FASTQ concatenates trivially too.
    """
    with open(destination, "wb") as out:
        for path in files:
            with open(path, "rb") as fh:
                shutil.copyfileobj(fh, out, length=8 * 1024 * 1024)
    return destination


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description="QC and clean Illumina NextSeq 1000 paired-end FASTQ files "
                    "for a somatic cancer pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io_group = parser.add_argument_group("input / output")
    # Both --input/--output (original) and --input-dir/--output-dir (alias,
    # for uniform chaining with align_reads.py / comprehensive_variant_calling.py)
    # are accepted. The two names write to the same dest.
    io_group.add_argument("-i", "--input", "--input-dir", dest="input_dir",
                          required=True,
                          help="Directory containing the raw FASTQ files.")
    io_group.add_argument("-o", "--output", "--output-dir", dest="output_dir",
                          required=True,
                          help="Directory for cleaned FASTQ files and reports.")
    io_group.add_argument("--r1-pattern", default="*_R1_*.fastq*",
                          help="Glob pattern matching Read 1 files.")
    io_group.add_argument("--r2-pattern", default="*_R2_*.fastq*",
                          help="Glob pattern matching Read 2 files.")
    io_group.add_argument("--recursive", action="store_true",
                          help="Search subdirectories recursively.")
    io_group.add_argument("--merge-lanes", action="store_true",
                          help="Concatenate lane-split files (_L001.._L004) "
                               "into one sample before cleaning.")
    io_group.add_argument("--overwrite", action="store_true",
                          help="Reprocess samples whose cleaned output already "
                               "exists (default: skip them).")
    io_group.add_argument("--dry-run", action="store_true",
                          help="Show what would be run without executing "
                               "fastp or FastQC.")

    filt = parser.add_argument_group("filtering / trimming")
    filt.add_argument("--average-qual", type=int, default=20,
                      help="Drop a read whose MEAN base quality is below this "
                           "(fastp -e/--average_qual). Set to 0 to disable.")
    filt.add_argument("--min-base-quality", "--qualified-quality",
                      dest="min_base_quality", type=int, default=15,
                      help="Phred score at or above which a base counts as "
                           "'qualified' (fastp -q). This is NOT a mean-quality "
                           "filter; it feeds --max-unqualified-pct. Same flag "
                           "name as comprehensive_variant_calling.py "
                           "(--qualified-quality still accepted).")
    filt.add_argument("--max-unqualified-pct", "--unqualified-percent",
                      dest="max_unqualified_pct", type=int, default=40,
                      help="Drop a read if more than this percentage of its "
                           "bases are unqualified (fastp -u). Same flag name "
                           "as comprehensive_variant_calling.py "
                           "(--unqualified-percent still accepted).")
    filt.add_argument("--min-read-length", "--min-length",
                      dest="min_read_length", type=int, default=50,
                      help="Minimum read length after trimming. Same flag "
                           "name as comprehensive_variant_calling.py "
                           "(--min-length still accepted).")
    filt.add_argument("--max-n-bases", "--n-base-limit",
                      dest="max_n_bases", type=int, default=5,
                      help="Drop a read with more than this many N bases. "
                           "Same flag name as comprehensive_variant_calling.py "
                           "(--n-base-limit still accepted).")
    filt.add_argument("--adapter-r1", default=None,
                      help="Explicit Read 1 adapter sequence.")
    filt.add_argument("--adapter-r2", default=None,
                      help="Explicit Read 2 adapter sequence.")
    filt.add_argument("--no-detect-adapter-pe", action="store_true",
                      help="Disable --detect_adapter_for_pe. By default it is "
                           "ON: overlap analysis alone misses read-through on "
                           "non-overlapping (long-insert) fragments.")
    filt.add_argument("--no-poly-g", action="store_true",
                      help="Disable poly-G trimming. It is ON by default "
                           "because NextSeq 1000 uses two-colour chemistry.")
    filt.add_argument("--poly-g-min-len", type=int, default=10,
                      help="Minimum poly-G run length to trim.")
    filt.add_argument("--cut-right", action="store_true",
                      help="Enable 3' sliding-window quality trimming. OFF by "
                           "default: it can distort allele fractions at read "
                           "ends in low-VAF somatic calling.")
    filt.add_argument("--cut-right-mean-quality", type=int, default=20,
                      help="Sliding-window mean quality, if --cut-right is set.")
    filt.add_argument("--correction", action="store_true",
                      help="Enable overlap-based base correction. OFF by "
                           "default: it interferes with UMI consensus calling.")

    umi = parser.add_argument_group("UMI (required for duplex/consensus panels)")
    umi.add_argument("--umi-loc", default=None,
                     choices=["index1", "index2", "read1", "read2",
                              "per_index", "per_read"],
                     help="Where the UMI lives. If your Celemics/FFPE panel "
                          "carries UMIs and you omit this, the UMI bases are "
                          "consumed as insert sequence and downstream "
                          "consensus calling is invalid.")
    umi.add_argument("--umi-len", type=int, default=None,
                     help="UMI length in bases (required with --umi-loc when "
                          "the UMI is in read1/read2).")
    umi.add_argument("--umi-skip", type=int, default=None,
                     help="Bases to skip after the UMI (e.g. a linker).")

    runtime = parser.add_argument_group("runtime")
    runtime.add_argument("--threads", type=int, default=8,
                         help="Threads for fastp (fastp caps usefully near 16).")
    runtime.add_argument("--skip-fastqc", action="store_true",
                         help="Skip FastQC on raw reads.")
    runtime.add_argument("--fastqc-cleaned", action="store_true",
                         help="Also run FastQC on the cleaned reads.")
    runtime.add_argument("--failed-out", action="store_true",
                         help="Write discarded reads to a per-sample file for "
                              "audit.")
    runtime.add_argument("--no-overrepresentation", action="store_true",
                         help="Disable fastp overrepresented-sequence analysis "
                              "(it is ON by default; useful for spotting "
                              "library artefacts, but costs runtime).")
    runtime.add_argument("--reads-to-process", type=int, default=None,
                         help="Only process this many reads from each FASTQ. "
                              "Useful for quick test runs / debugging. Same "
                              "flag as comprehensive_variant_calling.py.")

    # --- Panel / assay profile ---
    # The same profiles the variant caller uses, applied to the subset this
    # stage implements: adapters, read length and the UMI layout. Handing
    # both stages one profile name is the point -- a QIAseq run whose UMI
    # was extracted here and a caller stage configured from somewhere else
    # is how a chain ends up half-configured.
    panel = parser.add_argument_group(
        "panel / assay profile",
        "Configure trimming and UMI handling by naming the assay. Explicit "
        "flags always win over the profile.")
    panel.add_argument("--panel", default=None, metavar="ID",
                       help="Panel profile to take QC settings from (UMI "
                            "location and length, adapters, minimum read "
                            "length). Only the trimming keys are read here; "
                            "the same profile also configures the variant "
                            "caller. --list-panels shows what is available.")
    panel.add_argument("--panel-file", action="append", default=[],
                       metavar="JSON",
                       help="Additional profile file or directory to load, "
                            "repeatable.")
    panel.add_argument("--list-panels", action="store_true",
                       help="List every known panel profile and exit.")
    panel.add_argument("--describe-panel", default=None, metavar="ID",
                       help="Print everything a profile sets, and exit.")
    return parser


def validate_args(args, parser):
    if args.umi_loc in ("read1", "read2", "per_read") and not args.umi_len:
        parser.error("--umi-len is required when --umi-loc is read1/read2/per_read.")
    if args.umi_len and not args.umi_loc:
        parser.error("--umi-len given without --umi-loc.")
    # build_fastp_command() only emits --umi_skip inside the `if args.umi_loc`
    # block, so without this check the flag would be accepted and then
    # silently ignored. Erroring matches comprehensive_variant_calling.py.
    if args.umi_skip is not None and not args.umi_loc:
        parser.error("--umi-skip given without --umi-loc.")
    if args.threads < 1:
        parser.error("--threads must be >= 1.")
    if not 0 <= args.max_unqualified_pct <= 100:
        parser.error("--max-unqualified-pct must be between 0 and 100.")


# ---------------------------------------------------------------------------
# fastp command construction
# ---------------------------------------------------------------------------
def build_fastp_command(args, sample, r1, r2, out_r1, out_r2, reports_dir):
    cmd = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "-o", out_r1,
        "-O", out_r2,
        "--thread", str(args.threads),
        "--length_required", str(args.min_read_length),
        "--qualified_quality_phred", str(args.min_base_quality),
        "--unqualified_percent_limit", str(args.max_unqualified_pct),
        "--n_base_limit", str(args.max_n_bases),
        "--html", os.path.join(reports_dir, f"{sample}.fastp.html"),
        "--json", os.path.join(reports_dir, f"{sample}.fastp.json"),
        "--report_title", f"fastp report: {sample}",
    ]

    # Mean-quality read filter. This is the one that matches "minimum mean base
    # quality"; -q/--qualified_quality_phred above does something different.
    if args.average_qual > 0:
        cmd += ["--average_qual", str(args.average_qual)]

    # Adapters: explicit sequences if supplied, otherwise detection. For PE
    # input fastp does overlap-based trimming and only runs adapter-sequence
    # detection when --detect_adapter_for_pe is given.
    if args.adapter_r1:
        cmd += ["--adapter_sequence", args.adapter_r1]
    if args.adapter_r2:
        cmd += ["--adapter_sequence_r2", args.adapter_r2]
    if not args.no_detect_adapter_pe:
        cmd += ["--detect_adapter_for_pe"]

    # Poly-G. fastp auto-enables this when it recognises a NextSeq/NovaSeq
    # instrument ID in the read header, but header sniffing fails on renamed or
    # reheadered files, so we force it explicitly.
    if not args.no_poly_g:
        cmd += ["--trim_poly_g", "--poly_g_min_len", str(args.poly_g_min_len)]
    else:
        cmd += ["--disable_trim_poly_g"]

    if args.cut_right:
        cmd += ["--cut_right",
                "--cut_right_mean_quality", str(args.cut_right_mean_quality)]
    if args.correction:
        cmd += ["--correction"]

    if args.umi_loc:
        cmd += ["--umi", "--umi_loc", args.umi_loc]
        if args.umi_len:
            cmd += ["--umi_len", str(args.umi_len)]
        if args.umi_skip is not None:
            cmd += ["--umi_skip", str(args.umi_skip)]

    if not args.no_overrepresentation:
        cmd += ["--overrepresentation_analysis"]
    if args.failed_out:
        cmd += ["--failed_out",
                os.path.join(reports_dir, f"{sample}.failed.fastq.gz")]
    if getattr(args, "reads_to_process", None):
        cmd += ["--reads_to_process", str(args.reads_to_process)]

    return cmd


def summarise_fastp_json(json_path):
    """Pull the handful of numbers worth keeping in the run manifest."""
    try:
        with open(json_path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        return {"error": f"could not read fastp JSON: {exc}"}

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    reads_before = before.get("total_reads")
    reads_after = after.get("total_reads")
    retained = None
    if reads_before:
        retained = round(100.0 * (reads_after or 0) / reads_before, 2)

    return {
        "reads_before": reads_before,
        "reads_after": reads_after,
        "percent_reads_retained": retained,
        "q30_rate_before": before.get("q30_rate"),
        "q30_rate_after": after.get("q30_rate"),
        "read1_mean_length_after": after.get("read1_mean_length"),
        "read2_mean_length_after": after.get("read2_mean_length"),
        "gc_content_after": after.get("gc_content"),
        "duplication_rate": data.get("duplication", {}).get("rate"),
        "insert_size_peak": data.get("insert_size", {}).get("peak"),
        "adapter_trimmed_reads":
            data.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
        "low_quality_reads": filtering.get("low_quality_reads"),
        "too_short_reads": filtering.get("too_short_reads"),
        "too_many_N_reads": filtering.get("too_many_N_reads"),
    }


def verify_output(path):
    """Confirm an output file exists and is not a truncated/empty write."""
    if not os.path.exists(path):
        return False, f"missing output: {path}"
    size = os.path.getsize(path)
    if size < MIN_PLAUSIBLE_OUTPUT_BYTES:
        return False, f"suspiciously small output ({size} bytes): {path}"
    return True, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def apply_panel(args, parser, argv):
    """
    Fill in the QC settings this stage implements from --panel.

    Only the trimming and UMI keys are applied: the same profile carries
    caller and report settings, which this script has no arguments for and
    must ignore rather than choke on. panel_profiles.py is imported lazily
    and its absence is tolerated unless --panel was actually asked for --
    a run configured from a profile that could not be loaded would be a
    differently configured run that said nothing about it.
    """
    if not args.panel:
        return
    try:
        import panel_profiles
    except ImportError as exc:
        parser.error(f"--panel needs panel_profiles.py beside this script, "
                     f"and it could not be imported ({exc}).")
        return
    try:
        registry = panel_profiles.available_panels(
            extra_files=args.panel_file,
            warn=lambda msg: print(f"[WARN] {msg}"))
        profile = panel_profiles.resolve_panel(args.panel, registry)
        applied, overridden = panel_profiles.apply_profile(
            args, profile, panel_profiles.explicit_dests(parser, argv),
            only=panel_profiles.QC_SETTINGS)
    except panel_profiles.PanelError as exc:
        parser.error(str(exc))
        return
    print()
    print(panel_profiles.format_application(profile, applied, overridden))


def main():
    # --list-panels / --describe-panel answer without -i/-o, which argparse
    # cannot express, so they are handled before the real parser runs.
    try:
        import panel_profiles
        panel_profiles.handle_panel_queries()
    except ImportError:
        pass

    parser = build_parser()
    args = parser.parse_args()
    # Before validation: the profile can supply the --umi-len that
    # validate_args() is about to insist on.
    apply_panel(args, parser, sys.argv[1:])
    validate_args(args, parser)

    started = datetime.now(timezone.utc)

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    if not os.path.isdir(input_dir):
        print(f"[ERROR] Input directory does not exist: {input_dir}")
        return 1

    # ---- Pre-flight tool checks ----
    if shutil.which("fastp") is None and not args.dry_run:
        print("[ERROR] 'fastp' not found on PATH.")
        return 1
    need_fastqc = (not args.skip_fastqc) or args.fastqc_cleaned
    if need_fastqc and shutil.which("fastqc") is None and not args.dry_run:
        print("[ERROR] 'fastqc' not found on PATH. Install it, or use "
              "--skip-fastqc and drop --fastqc-cleaned.")
        return 1

    # ---- Output layout ----
    cleaned_dir = os.path.join(output_dir, "cleaned_fastq")
    fastqc_raw_dir = os.path.join(output_dir, "fastqc_raw")
    fastqc_clean_dir = os.path.join(output_dir, "fastqc_cleaned")
    reports_dir = os.path.join(output_dir, "fastp_reports")
    logs_dir = os.path.join(output_dir, "logs")
    for directory in (cleaned_dir, reports_dir, logs_dir):
        os.makedirs(directory, exist_ok=True)
    if not args.skip_fastqc:
        os.makedirs(fastqc_raw_dir, exist_ok=True)
    if args.fastqc_cleaned:
        os.makedirs(fastqc_clean_dir, exist_ok=True)

    # ---- Read tokens ----
    r1_token, r2_token = derive_tokens(args.r1_pattern, args.r2_pattern)
    if not r1_token:
        print("[ERROR] Could not derive R1/R2 tokens from the patterns "
              f"'{args.r1_pattern}' and '{args.r2_pattern}'. They must differ "
              "in exactly one literal segment, e.g. '*_R1_*.fastq*' and "
              "'*_R2_*.fastq*'.")
        return 1

    print(f"[INFO] Input      : {input_dir}")
    print(f"[INFO] Output     : {output_dir}")
    print(f"[INFO] Read tokens: '{r1_token}' / '{r2_token}'")

    # ---- Discovery ----
    pairs, warnings = find_paired_fastqs(
        input_dir, args.r1_pattern, r1_token, r2_token, args.recursive)
    for warning in warnings:
        print(f"[WARN] {warning}")
    if not pairs:
        print("[ERROR] No FASTQ pairs found. Check --input, --r1-pattern and "
              "whether you need --recursive.")
        return 1

    # ---- Lane grouping ----
    groups, lane_split = group_by_lane(pairs)
    if lane_split and not args.merge_lanes:
        print("[WARN] Lane-split files detected (_L001.._L004). Without "
              "--merge-lanes each lane is processed as a separate sample and "
              "they are never combined. This is usually not what you want.")
    if args.merge_lanes:
        units = [{"sample": name, "parts": parts}
                 for name, parts in sorted(groups.items())]
    else:
        units = [{"sample": p["sample"], "parts": [p]} for p in pairs]

    print(f"[INFO] {len(units)} sample(s) to process.")

    manifest = {
        "script": os.path.basename(__file__),
        "started_utc": started.isoformat(),
        "host": os.uname().nodename if hasattr(os, "uname") else None,
        "python": sys.version.split()[0],
        "command_line": shlex.join(sys.argv),
        "parameters": vars(args),
        "tool_versions": {
            "fastp": tool_version("fastp"),
            "fastqc": tool_version("fastqc"),
        },
        "input_dir": input_dir,
        "output_dir": output_dir,
        "read_tokens": {"r1": r1_token, "r2": r2_token},
        "discovery_warnings": warnings,
        "samples": [],
    }

    succeeded, failed, skipped = [], [], []
    temp_dir = None

    for index, unit in enumerate(units, start=1):
        sample = unit["sample"]
        parts = unit["parts"]
        log_path = os.path.join(logs_dir, f"{sample}.log")
        record = {
            "sample": sample,
            "inputs": [{"r1": p["r1"], "r2": p["r2"]} for p in parts],
            "lanes_merged": len(parts) > 1,
            "log": log_path,
        }

        print(f"\n{'=' * 72}")
        print(f"[{index}/{len(units)}] Sample: {sample}")
        print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'=' * 72}")

        out_r1 = os.path.join(cleaned_dir, f"{sample}_R1.clean.fastq.gz")
        out_r2 = os.path.join(cleaned_dir, f"{sample}_R2.clean.fastq.gz")

        if os.path.exists(out_r1) and os.path.exists(out_r2) and not args.overwrite:
            print("[SKIP] Cleaned output already exists. Use --overwrite to "
                  "reprocess.")
            record["status"] = "skipped_existing"
            manifest["samples"].append(record)
            skipped.append(sample)
            continue

        # -- Provenance on the raw inputs --
        if not args.dry_run:
            record["input_provenance"] = []
            for part in parts:
                for role in ("r1", "r2"):
                    path = part[role]
                    entry = {"path": path, "size_bytes": os.path.getsize(path)}
                    entry.update(file_digest(path))
                    record["input_provenance"].append(entry)

        # -- Merge lanes if requested --
        r1_in, r2_in = parts[0]["r1"], parts[0]["r2"]
        if len(parts) > 1:
            if temp_dir is None:
                temp_dir = tempfile.mkdtemp(prefix="fastq_qc_merge_",
                                            dir=output_dir)
            print(f"[INFO] Merging {len(parts)} lanes for {sample}.")
            if not args.dry_run:
                r1_in = concatenate([p["r1"] for p in parts],
                                    os.path.join(temp_dir, f"{sample}_R1.fastq.gz"))
                r2_in = concatenate([p["r2"] for p in parts],
                                    os.path.join(temp_dir, f"{sample}_R2.fastq.gz"))

        # -- FastQC on raw reads --
        if not args.skip_fastqc:
            code = run_command(
                ["fastqc", r1_in, r2_in, "-o", fastqc_raw_dir,
                 "--threads", str(min(2, args.threads)), "--quiet"],
                tag=f"FastQC raw [{sample}]",
                log_path=log_path, dry_run=args.dry_run)
            if code != 0:
                # FastQC is diagnostic only; a failure here is not fatal.
                print(f"[WARN] FastQC failed for {sample}; continuing.")
                record["fastqc_raw_exit_code"] = code

        # -- fastp --
        fastp_cmd = build_fastp_command(args, sample, r1_in, r2_in,
                                        out_r1, out_r2, reports_dir)
        record["fastp_command"] = shlex.join(fastp_cmd)
        code = run_command(fastp_cmd, tag=f"fastp [{sample}]",
                           log_path=log_path, dry_run=args.dry_run)

        if code != 0:
            record["status"] = "failed"
            record["error"] = f"fastp exit code {code}"
            manifest["samples"].append(record)
            failed.append(sample)
            print(f"[FAIL] {sample}: fastp exited {code}. Continuing with the "
                  f"remaining samples.")
            continue

        # -- Verify outputs really exist and are plausible --
        if not args.dry_run:
            problems = []
            for path in (out_r1, out_r2):
                ok, message = verify_output(path)
                if not ok:
                    problems.append(message)
            if problems:
                record["status"] = "failed"
                record["error"] = "; ".join(problems)
                manifest["samples"].append(record)
                failed.append(sample)
                print(f"[FAIL] {sample}: {record['error']}")
                continue

            record["outputs"] = {
                "r1": out_r1, "r2": out_r2,
                "r1_bytes": os.path.getsize(out_r1),
                "r2_bytes": os.path.getsize(out_r2),
            }
            record["metrics"] = summarise_fastp_json(
                os.path.join(reports_dir, f"{sample}.fastp.json"))

            metrics = record["metrics"]
            if metrics.get("percent_reads_retained") is not None:
                print(f"[INFO] {sample}: "
                      f"{metrics['reads_before']} -> {metrics['reads_after']} "
                      f"reads ({metrics['percent_reads_retained']}% retained), "
                      f"Q30 {metrics.get('q30_rate_before')} -> "
                      f"{metrics.get('q30_rate_after')}")

        # -- FastQC on cleaned reads --
        if args.fastqc_cleaned:
            code = run_command(
                ["fastqc", out_r1, out_r2, "-o", fastqc_clean_dir,
                 "--threads", str(min(2, args.threads)), "--quiet"],
                tag=f"FastQC cleaned [{sample}]",
                log_path=log_path, dry_run=args.dry_run)
            if code != 0:
                print(f"[WARN] FastQC (cleaned) failed for {sample}; continuing.")
                record["fastqc_cleaned_exit_code"] = code

        record["status"] = "ok"
        manifest["samples"].append(record)
        succeeded.append(sample)

    # ---- Clean up merge temporaries ----
    if temp_dir and os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)

    # ---- Manifest ----
    finished = datetime.now(timezone.utc)
    manifest["finished_utc"] = finished.isoformat()
    manifest["duration_seconds"] = round((finished - started).total_seconds(), 1)
    manifest["result"] = {
        "succeeded": succeeded,
        "failed": failed,
        "skipped": skipped,
    }
    manifest_path = os.path.join(
        output_dir, f"run_manifest_{started.strftime('%Y%m%dT%H%M%SZ')}.json")
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    # ---- Summary ----
    print(f"\n{'=' * 72}")
    print(f"Succeeded: {len(succeeded)}   Failed: {len(failed)}   "
          f"Skipped: {len(skipped)}")
    if failed:
        print(f"Failed samples: {', '.join(failed)}")
    print(f"Cleaned reads : {cleaned_dir}")
    print(f"fastp reports : {reports_dir}")
    print(f"Per-sample logs: {logs_dir}")
    print(f"Run manifest  : {manifest_path}")
    print(f"Aggregate with: multiqc {output_dir} -o {output_dir}")
    print(f"{'=' * 72}")

    # Non-zero exit if anything failed, so a scheduler or Nextflow wrapper
    # notices.
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())