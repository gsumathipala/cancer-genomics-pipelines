#!/usr/bin/env python3
"""
pipeline_orchestrator.py
========================
Orchestrator that chains the three stages of the cancer DNA pipeline
into a single command:

    STAGE 1: fastq_qc_clean.py
             QC + adapter/poly-G trimming -> cleaned FASTQs + QC manifest.

    STAGE 2 (optional): align_reads.py
             Alignment to human reference -> sorted/indexed BAMs.
             Skip this if you only want variant calling output.

    STAGE 3: comprehensive_variant_calling.py
             Duplicate marking -> BQSR -> Mutect2 -> Filter -> COSMIC ->
             SnpEff -> PCGR clinical report (the last of these only when
             --pcgr-refdata-dir is passed through in --stage3-args).
             When a QC manifest is passed via --manifest, the "qc" step is
             auto-skipped and the cleaned reads from STAGE 1 are reused
             (no redundant re-QC).

THE HANDOFF MECHANISM
---------------------
  The manifest written by fastq_qc_clean.py (run_manifest_*.json) records the
  exact path of each cleaned FASTQ. align_reads.py and
  comprehensive_variant_calling.py both accept --manifest and read the
  cleaned FASTQ paths + sample names from it. This gives a clean,
  verifiable handoff between stages -- no fragile filename guessing.

  If STAGE 2 is SKIPPED, the QC manifest goes straight to STAGE 3, which
  aligns the cleaned reads itself (chain: QC -> variant calling).

  If STAGE 2 RUNS, its BAMs would otherwise be wasted: STAGE 3 re-aligns
  from the cleaned FASTQs whenever it is handed a QC manifest, so BWA-MEM2
  would run twice. To prevent that, the orchestrator appends

      --skip-steps align --aligned-dir <stage2-output-dir>/aligned

  to STAGE 3's arguments, pointing it at the BAMs stage 2 just produced.
  STAGE 3 still receives --input-dir (the QC cleaned_fastq directory) and
  --manifest; in manifest mode --input-dir is not actually read, but it is
  passed for consistency with the other invocation modes.

Usage Examples
--------------
  # Full 3-stage pipeline (QC -> align -> variant call) with COSMIC:
  python pipeline_orchestrator.py \\
      --stage1-args "-i raw_fastqs/ -o qc_out/ --threads 16" \\
      --stage2-args "-o aligned_out/ --reference hg38 --threads 16 --two-pass" \\
      --stage3-args "-o results/ --reference hg38 \\
                     --cosmic /data/COSMIC_v98.vcf.gz \\
                     --dbsnp /data/dbsnp.vcf.gz \\
                     --threads 16"

  # QC -> variant calling (skip alignment stage):
  python pipeline_orchestrator.py \\
      --stage1-args "-i raw_fastqs/ -o qc_out/ --threads 16" \\
      --skip-stage2 \\
      --stage3-args "-o results/ --reference hg38 \\
                     --cosmic /data/COSMIC_v98.vcf.gz --threads 16"

  # Reuse a PREVIOUS QC run (don't re-run fastp):
  python pipeline_orchestrator.py \\
      --skip-stage1 \\
      --qc-manifest qc_out/run_manifest_20240101T000000Z.json \\
      --stage3-args "-o results/ --reference hg38 \\
                     --cosmic /data/COSMIC_v98.vcf.gz --threads 16"

NOTES
-----
  * Stage 3 requires the newest comprehensive_variant_calling.py with
    --manifest support (this repo).
  * The exact command for each stage is PRINTED to the console just before
    that stage runs, and each stage writes its own per-sample logs beneath
    its own output directory. --log-dir holds only this orchestrator's
    pipeline_report_*.json summary -- there is no command transcript file.

DEBUGGING NOTES (things that surprise people)
---------------------------------------------
  * THE FULL STAGE-3 COMMAND IS PRINTED before it runs. When a stage
    misbehaves, copy that line and run it directly -- the orchestrator adds
    arguments of its own (--manifest, --input-dir, and when stage 2 ran,
    --skip-steps align --aligned-dir), so the command you typed is NOT the
    command that executed.

  * STAGE 2 -> STAGE 3 HANDOFF: if alignment ran, stage 3 is told to skip
    its own alignment and reuse those BAMs. Without that, BWA-MEM2 would
    run twice (hours of wasted compute) because stage 3 re-aligns from the
    cleaned FASTQs whenever it is handed a QC manifest.

  * add_skip_step() REWRITES the user's --skip-steps rather than appending
    a second one, because argparse's nargs="*" keeps only the LAST
    occurrence. It handles both "--skip-steps=x" and "--skip-steps x y".
    Appending naively silently discards whatever the user asked to skip.

  * PRECEDENCE when several inputs are given:
      --skip-stage1 wins over --stage1-args (stage 1 will not run).
      --qc-manifest wins over --stage1-args (stage 1 is skipped, with a
      [WARN]) -- if stage 1 unexpectedly did not run, look for that line.
      --qc-output-dir resolves to the NEWEST run_manifest_*.json in it,
      relying on the timestamped filenames sorting chronologically.

  * Stage exit codes are fatal: any non-zero stage aborts the chain with
    sys.exit(1), so a later stage never runs on a half-built input.
"""

import argparse
import glob
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))

# Canonical script names -- the orchestrator expects these files to exist
# in the same directory as this file.
SCRIPT_QC = os.path.join(SCRIPTS_DIR, "fastq_qc_clean.py")
SCRIPT_ALIGN = os.path.join(SCRIPTS_DIR, "align_reads.py")
SCRIPT_VARIANT = os.path.join(SCRIPTS_DIR, "comprehensive_variant_calling.py")


def run_stage(script, args_str, tag, dry_run=False):
    """
    Run one pipeline stage as a subprocess.

    The args_str is a space-separated argument string (shlex parsed) so
    that the orchestrator stays simple and the stage's own argparse does
    the real validation.

    RETURNS:
        int: Exit code of the stage (0 = success).
    """
    cmd = [sys.executable, script] + shlex.split(args_str)
    printable = shlex.join(cmd)
    print(f"\n{'#' * 72}")
    print(f"# {tag}")
    print(f"# $ {printable}")
    print(f"{'#' * 72}", flush=True)

    if dry_run:
        print(f"[{tag}] (dry run - not executed)")
        return 0

    # Stream the subprocess output to the console in real time, the same
    # way the pipeline scripts themselves do. Merge stderr so nothing is
    # lost, and ever-dangerous buffering is avoided.
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
    proc.wait()

    if proc.returncode != 0:
        print(f"[{tag}] ERROR: exit code {proc.returncode}")
    return proc.returncode


def find_latest_manifest(qc_output_dir):
    """
    Find the most recent QC run manifest in the QC output directory.

    fastq_qc_clean.py writes manifests named run_manifest_<timestamp>.json.
    This lets the orchestrator pick up the newest QC run so STAGE 2/3 can
    consume it -- even if the user did not pass --qc-manifest.

    RETURNS:
        str: Path to the newest manifest, or None if none found.
    """
    pattern = os.path.join(qc_output_dir, "run_manifest_*.json")
    matches = sorted(glob_manifests(pattern))
    return matches[-1] if matches else None


def glob_manifests(pattern):
    """Minimal glob wrapper; returns a sorted list so that 'newest is last'."""
    return sorted(glob.glob(pattern))


def build_manifest_records(qc_manifest_path, tumour_hint=None, normal_hint=None):
    """
    Summarise what the QC manifest tells us about samples, for the final
    orchestrator report. Returns tuple (tumour_name, normal_name, sample_count).

    `tumour_hint`/`normal_hint` are the --tumour-sample/--normal-sample
    values (if any) pulled from --stage3-args: when Stage 3 was told
    explicitly which sample is which, the report should reflect that
    instead of guessing tumour = first sample / normal = second sample.
    """
    try:
        with open(qc_manifest_path, encoding="utf-8") as fh:
            data = json.load(fh)
        samples = [s for s in data.get("samples", [])
                   if s.get("status") != "failed"]
        names = [s.get("sample") for s in samples]
        tumour = tumour_hint if tumour_hint in names else (
            names[0] if names else None)
        remaining = [n for n in names if n != tumour]
        normal = normal_hint if normal_hint in names else (
            remaining[0] if remaining else None)
        return tumour, normal, len(names)
    except Exception as exc:
        print(f"[WARN] Could not read QC manifest metadata: {exc}")
        return None, None, 0


def main():
    parser = argparse.ArgumentParser(
        description="Chain fastq_qc_clean.py -> align_reads.py -> "
                    "comprehensive_variant_calling.py into one command.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Stage 1 (QC) ---
    s1 = parser.add_argument_group("stage 1 - QC (fastq_qc_clean.py)")
    s1.add_argument("--stage1-args",
                    help="Arguments passed to fastq_qc_clean.py, e.g. "
                         "'-i raw/ -o qc_out/ --threads 16'. Optional if "
                         "--qc-manifest is given (reuse previous QC).")
    s1.add_argument("--skip-stage1", action="store_true",
                    help="Skip running fastq_qc_clean.py. Requires "
                         "--qc-manifest (or --qc-output-dir + a previous run).")

    # --- Stage 2 (Alignment) ---
    s2 = parser.add_argument_group("stage 2 - alignment (align_reads.py)")
    s2.add_argument("--stage2-args",
                    help="Arguments passed to align_reads.py, e.g. "
                         "'-o aligned_out/ --reference hg38 --threads 16'. "
                         "Leave empty to skip stage 2.")
    s2.add_argument("--skip-stage2", action="store_true",
                    help="Skip alignment. The QC manifest is then passed "
                         "directly to variant calling.")

    # --- Stage 3 (Variant calling) ---
    s3 = parser.add_argument_group("stage 3 - variant calling")
    s3.add_argument("--stage3-args", required=True,
                    help="Arguments passed to "
                         "comprehensive_variant_calling.py, e.g. "
                         "'-o results/ --reference hg38 --cosmic ... --threads 16'. "
                         "Do NOT pass --manifest here; the orchestrator adds it.")

    # --- Handoff ---
    h = parser.add_argument_group("handoff / resume")
    h.add_argument("--qc-manifest", default=None,
                   help="Path to an existing fastq_qc_clean.py run manifest. "
                        "If given, stage 1 is skipped and this manifest is "
                        "used for the handoff.")
    h.add_argument("--qc-output-dir", default=None,
                   help="Directory containing a previous QC run. The newest "
                        "run_manifest_*.json in it is used for the handoff.")

    # --- Runtime ---
    r = parser.add_argument_group("runtime")
    r.add_argument("--dry-run", action="store_true",
                   help="Print the commands without running them.")
    r.add_argument("--log-dir", default="pipeline_logs",
                   help="Directory for orchestrator log files.")

    args = parser.parse_args()
    dry_run = args.dry_run

    # ---- Validate combinatorial inputs ----
    # Stage 3 always needs a QC manifest. That manifest comes from either:
    #   1. A previous run (--qc-manifest / --qc-output-dir), or
    #   2. Stage 1 running now (--stage1-args).
    # Stage 2 is optional and consumes the same QC manifest.
    has_qc_handoff = bool(args.qc_manifest or args.qc_output_dir)
    if args.skip_stage1 and not has_qc_handoff:
        parser.error("--skip-stage1 requires --qc-manifest or "
                     "--qc-output-dir (a previous QC run to reuse).")
    if args.skip_stage2 and not has_qc_handoff:
        print("[INFO] Stage 2 skipped; stage 3 will consume the QC manifest "
              "produced by stage 1 (" + (args.stage1_args or "none") + ").")
        if not args.stage1_args:
            parser.error("--skip-stage2 and no --stage1-args means nothing "
                         "will produce input reads for stage 3. Provide "
                         "--stage1-args or --qc-manifest/--qc-output-dir.")

    if not args.stage1_args and not has_qc_handoff:
        parser.error("Provide --stage1-args, OR a prior QC handoff via "
                     "--qc-manifest/--qc-output-dir.")

    started = datetime.now(timezone.utc)
    log_dir = args.log_dir
    os.makedirs(log_dir, exist_ok=True)

    # =====================================================================
    # STAGE 1: QC & CLEANING
    # =====================================================================
    qc_manifest = args.qc_manifest

    if args.skip_stage1:
        # --skip-stage1 always wins: never run stage 1, even if --stage1-args
        # was also given (previously this flag was parsed but never
        # consulted, so it silently did nothing in that combination).
        if args.stage1_args:
            print("[INFO] --skip-stage1 given; ignoring --stage1-args and "
                  "reusing the QC handoff instead.")
        if qc_manifest is None and args.qc_output_dir:
            qc_manifest = find_latest_manifest(args.qc_output_dir)
        if not qc_manifest:
            parser.error("--skip-stage1 requires --qc-manifest or "
                         "--qc-output-dir (a previous QC run to reuse).")
        print("[SKIP] Stage 1 (QC) skipped (--skip-stage1).")
    else:
        if args.stage1_args and qc_manifest:
            print("[WARN] Both --stage1-args and --qc-manifest were given; "
                  "Stage 1 will be skipped and the QC manifest reused as-is. "
                  "Drop --qc-manifest if you want Stage 1 to run.")

        if args.stage1_args and not qc_manifest:
            code = run_stage(SCRIPT_QC, args.stage1_args, "STAGE 1: QC (fastp)",
                             dry_run=dry_run)
            if code != 0:
                print("[FATAL] Stage 1 (QC) failed.")
                sys.exit(1)

            if dry_run:
                # In dry-run mode no manifest exists yet; fall back to the
                # user-provided handoff or print what would be discovered.
                qc_manifest = args.qc_manifest or args.qc_output_dir or \
                    "(produced by stage 1 when run for real)"
            else:
                # Locate the manifest the QC stage just wrote.
                qc_output_dir = extract_arg(args.stage1_args,
                                            ("-o", "--output", "--output-dir"))
                if not qc_output_dir:
                    print("[ERROR] Could not determine the QC output dir from "
                          "--stage1-args. Pass --qc-output-dir explicitly.")
                    sys.exit(1)
                qc_manifest = find_latest_manifest(qc_output_dir)
                if not qc_manifest:
                    print(f"[ERROR] No run_manifest_*.json found in {qc_output_dir}. "
                          f"Did stage 1 run successfully?")
                    sys.exit(1)
        elif qc_manifest is None and args.qc_output_dir:
            qc_manifest = find_latest_manifest(args.qc_output_dir)
            if not qc_manifest:
                print(f"[ERROR] No run_manifest_*.json found in "
                      f"{args.qc_output_dir}.")
                sys.exit(1)

    print(f"\n[INFO] Using QC manifest: {qc_manifest}")

    # =====================================================================
    # STAGE 2: ALIGNMENT (optional)
    # =====================================================================
    if not args.skip_stage2 and args.stage2_args:
        # Pass the --manifest and --input-dir (cleaned reads dir) to align.
        # We must locate the cleaned_fastq directory. It lives next to the
        # manifest as <qc_output>/cleaned_fastq/.
        if dry_run and not os.path.isfile(qc_manifest):
            cleaned_dir = "<QC cleaned_fastq dir>"
        else:
            qc_root = os.path.dirname(os.path.abspath(qc_manifest))
            # Prefer the manifest's own output_dir; fall back to dirname.
            if os.path.isfile(qc_manifest):
                with open(qc_manifest, encoding="utf-8") as fh:
                    _data = json.load(fh)
                if _data.get("output_dir"):
                    qc_root = _data["output_dir"]
            cleaned_dir = os.path.join(qc_root, "cleaned_fastq")
        io = shlex.split(args.stage2_args)
        # If the user did not supply --input-dir, inject it.
        if not has_option(io, ("-i", "--input-dir")):
            io += ["--input-dir", cleaned_dir]
        io += ["--manifest", qc_manifest]
        code = run_stage(SCRIPT_ALIGN, shlex.join(io),
                         "STAGE 2: Alignment (BWA-MEM2)", dry_run=dry_run)
        if code != 0:
            print("[FATAL] Stage 2 (alignment) failed.")
            sys.exit(1)
        align_output = extract_arg(args.stage2_args, ("-o", "--output-dir"))
    else:
        print("[SKIP] Stage 2 (alignment) skipped.")
        align_output = None

    # =====================================================================
    # STAGE 3: VARIANT CALLING
    # =====================================================================
    # Pass the QC manifest to comprehensive_variant_calling.py. Its qc step
    # is auto-skipped when --manifest is present, so cleaned reads from
    # stage 1 are reused directly.
    v3 = shlex.split(args.stage3_args)

    # If Stage 2 (alignment) ran, tell Stage 3 to reuse its BAMs instead of
    # re-aligning from the cleaned FASTQs -- otherwise BWA-MEM2 silently
    # runs twice (once in Stage 2, once again inside Stage 3).
    if align_output:
        aligned_dir = os.path.join(os.path.abspath(align_output), "aligned")
        v3 = add_skip_step(v3, "align")
        if not has_option(v3, ("--aligned-dir",)):
            v3 += ["--aligned-dir", aligned_dir]
        print(f"[INFO] Stage 3 will reuse Stage 2's alignment output: "
              f"{aligned_dir}")

    if not has_option(v3, ("-i", "--input-dir")):
        if dry_run and not os.path.isfile(qc_manifest):
            v3 += ["--input-dir", "<QC cleaned_fastq dir>"]
        else:
            # Point at the cleaned reads dir so the manifest's (possibly
            # relative) cleaned FASTQ paths resolve.
            qc_root = os.path.dirname(os.path.abspath(qc_manifest))
            if os.path.isfile(qc_manifest):
                with open(qc_manifest, encoding="utf-8") as fh:
                    _data = json.load(fh)
                if _data.get("output_dir"):
                    qc_root = _data["output_dir"]
            v3 += ["--input-dir", os.path.join(qc_root, "cleaned_fastq")]
    v3 += ["--manifest", qc_manifest]
    code = run_stage(SCRIPT_VARIANT, shlex.join(v3),
                     "STAGE 3: Variant Calling (Mutect2 + COSMIC + SnpEff)",
                     dry_run=dry_run)
    if code != 0:
        print("[FATAL] Stage 3 (variant calling) failed.")
        sys.exit(1)

    # =====================================================================
    # FINAL REPORT
    # =====================================================================
    finished = datetime.now(timezone.utc)
    tumour_hint = extract_arg(args.stage3_args, ("--tumour-sample",))
    normal_hint = extract_arg(args.stage3_args, ("--normal-sample",))
    if os.path.isfile(qc_manifest):
        tumour, normal, n_samples = build_manifest_records(
            qc_manifest, tumour_hint, normal_hint)
    elif dry_run:
        tumour, normal, n_samples = tumour_hint, normal_hint, "n/a (dry run)"

    report = {
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "duration_seconds": round((finished - started).total_seconds(), 1),
        "stages": {
            "1_qc": "ran" if qc_manifest else "skipped",
            "2_align": "ran" if align_output else "skipped",
            "3_variant_calling": "ran",
        },
        "qc_manifest": qc_manifest,
        "samples": {
            "count": n_samples,
            "tumour": tumour,
            "normal": normal,
        },
        "command_line": shlex.join(sys.argv),
    }

    report_path = os.path.join(
        log_dir,
        f"pipeline_report_{started.strftime('%Y%m%dT%H%M%SZ')}.json")
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print(f"\n{'=' * 72}")
    print("PIPELINE COMPLETE")
    print(f"{'=' * 72}")
    print(f"Duration   : {report['duration_seconds']}s")
    print(f"Samples    : {n_samples} (tumour: {tumour}, "
          f"normal: {normal or '-none-'})")
    print(f"QC manifest: {qc_manifest}")
    # Label this accurately: the directory holds the JSON report only, not
    # a transcript of the stage commands (those are printed, not saved).
    print(f"Log dir    : {log_dir}")
    print(f"Report     : {report_path}")
    print(f"{'=' * 72}")

    return 0


def has_option(tokens, option_names):
    """Return True if any of `option_names` appears in `tokens`, either as
    '-o value' or '-o=value'."""
    for token in tokens:
        for opt in option_names:
            if token == opt or token.startswith(opt + "="):
                return True
    return False


def add_skip_step(tokens, step):
    """
    Ensure `step` is included in a --skip-steps list within `tokens`,
    merging with any --skip-steps the user already passed in --stage3-args.

    WHY THIS IS FIDDLY:
      comprehensive_variant_calling.py declares --skip-steps with argparse's
      nargs="*". Repeating the option does NOT accumulate -- the last
      occurrence wins. So naively appending '--skip-steps align' would throw
      away whatever the user asked to skip. We therefore rewrite the user's
      existing list in place instead of adding a second one.

      There are two spellings to handle, and missing either one silently
      drops the user's request:
        1. "--skip-steps=bqsr"      (one token; argparse binds a single value)
        2. "--skip-steps bqsr qc"   (flag token followed by 0+ value tokens)

    Returns a NEW token list; the input is not mutated.
    """
    out = list(tokens)

    # --- Form 1: "--skip-steps=bqsr" -------------------------------------
    # Rewrite into the space form so we can append our step alongside the
    # user's. ("--skip-steps=a b" is not a way to pass two values, so the
    # equals form always carries exactly one.)
    for i, token in enumerate(out):
        if token.startswith("--skip-steps="):
            existing = token.split("=", 1)[1]
            values = [existing] if existing else []
            if step not in values:
                values.append(step)
            return out[:i] + ["--skip-steps"] + values + out[i + 1:]

    # --- Form 2: "--skip-steps bqsr qc" ----------------------------------
    if "--skip-steps" in out:
        idx = out.index("--skip-steps")
        # Values run until the next option token. Step names never start
        # with "-", so a leading dash reliably marks the end of the list.
        j = idx + 1
        values = []
        while j < len(out) and not out[j].startswith("-"):
            values.append(out[j])
            j += 1
        if step not in values:
            values.append(step)
        return out[:idx] + ["--skip-steps"] + values + out[j:]

    # --- Form 3: the user passed no --skip-steps at all -------------------
    return out + ["--skip-steps", step]


def extract_arg(args_str, option_names):
    """
    Extract the value of an option from an argument string.

    Handles both '-o value' and '-o=value' syntax.
    Returns None if the option is not present.
    """
    tokens = shlex.split(args_str)
    for i, token in enumerate(tokens):
        for opt in option_names:
            if token == opt and i + 1 < len(tokens):
                return tokens[i + 1]
            if token.startswith(opt + "="):
                return token.split("=", 1)[1]
    return None


if __name__ == "__main__":
    sys.exit(main())