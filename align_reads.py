#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
align_reads.py
==============
Align paired-end FASTQ files from Illumina NextSeq 1000 to the human
reference genome (hg38) using BWA-MEM2.

This is a STANDALONE alignment script -- it does NOT perform variant
calling, duplicate marking, or any downstream analysis. It produces
sorted, indexed BAM files ready for further processing (e.g. by
variant_calling.py or GATK).

What this script DOES
---------------------
  1. Validates inputs: FASTQ files exist, reference genome is accessible.
  2. Checks that BWA-MEM2 index, samtools .fai index, and GATK sequence
     dictionary exist (builds them if missing).
  3. Aligns paired-end reads with BWA-MEM2, using the two-pass alignment
     strategy for improved mapping of difficult regions.
  4. Pipes directly into samtools sort (no intermediate unsorted BAM).
  5. Creates a BAM index (.bai) for random access.
  6. Generates alignment statistics (samtools flagstat, idxstats, stats).
  7. Optionally runs FastQC on the resulting BAM for quality assessment.
  8. Writes a JSON manifest recording tool versions, parameters, and
     alignment metrics for reproducibility.

What this script DELIBERATELY DOES NOT DO
-----------------------------------------
  * No duplicate marking (use samtools markdup or GATK MarkDuplicates).
  * No base quality recalibration (use GATK BQSR).
  * No variant calling (use GATK Mutect2 or similar).
  * No adapter trimming (use fastp or fastq_qc_clean.py first).
  * No colour-space conversion or format conversion beyond BAM.

Requirements
------------
  - Python 3.8+
  - BWA-MEM2  (https://github.com/bwa-mem2/bwa-mem2)
  - samtools  (https://github.com/samtools/samtools)
  - (optional) GATK 4+   -- only needed if creating the sequence dictionary
  - (optional) FastQC     -- for alignment QC plots

BWA-MEM2 two-pass alignment
----------------------------
  BWA-MEM2 supports a two-pass mode where:
    Pass 1: Align reads normally, collect statistics on difficult regions.
    Pass 2: Re-align reads in difficult regions with adjusted parameters.

  Two-pass mode improves mapping sensitivity for:
    - Reads spanning splice junctions (RNA-seq, but also genomic).
    - Reads in repetitive or low-complexity regions.
    - Reads with structural variations.

  For cancer DNA, two-pass is recommended because tumour genomes often
  contain rearranged regions that are harder to map.

Usage examples
--------------
  # Basic alignment (single sample):
  python align_reads.py \\
      --input-dir cleaned_fastq/ \\
      --output-dir aligned/ \\
      --reference hg38.fa \\
      --tumour-sample TUMOUR_01 \\
      --tumour-r1 T1_R1.clean.fastq.gz \\
      --tumour-r2 T1_R2.clean.fastq.gz

  # With matched normal:
  python align_reads.py \\
      --input-dir cleaned_fastq/ \\
      --output-dir aligned/ \\
      --reference hg38.fa \\
      --tumour-sample T1 --tumour-r1 T1_R1.gz --tumour-r2 T1_R2.gz \\
      --normal-sample N1 --normal-r1 N1_R1.gz --normal-r2 N1_R2.gz \\
      --threads 16 \\
      --two-pass

  # Using hg38 auto-download:
  python align_reads.py \\
      --input-dir cleaned_fastq/ \\
      --output-dir aligned/ \\
      --reference hg38 \\
      --tumour-sample SAMPLE_01 \\
      --tumour-r1 R1.fastq.gz --tumour-r2 R2.fastq.gz

DEBUGGING NOTES (things that surprise people)
---------------------------------------------
  * OUTPUT LAYOUT IS LOAD-BEARING. BAMs go to <output-dir>/aligned/ named
    <sample>.sorted.bam, because comprehensive_variant_calling.py looks for
    exactly that (see its --aligned-dir option). Renaming the subdirectory
    silently breaks the QC -> align -> variant-calling chain.

  * RESUME BEHAVIOUR: a sample whose .bam AND .bai already exist is skipped
    (status "skipped_existing") unless --overwrite is given. Both files must
    be present -- a BAM without an index is treated as incomplete work.
    This check runs in dry-run mode too, so a preview tells the truth about
    what a real run would skip.

  * --dry-run performs NO I/O of consequence: no alignment, no reference
    download, no hashing of inputs. It still prints every command. If you
    ever see a dry run touch the network or write a big file, that is a bug
    -- the reference download and the BAM writes are both gated on it.

  * BWA STDERR IS DRAINED BY A THREAD (see align_sample). BWA writes
    progress to stderr while its stdout is piped into `samtools sort`. If
    you "simplify" this by reading stderr only after wait(), a long
    alignment fills the ~64 KB pipe buffer and the whole pipeline
    deadlocks. Do not merge stderr into stdout either: stdout is the SAM
    stream and merging would corrupt the BAM.

  * MANIFEST MODE: with --manifest, sample names and FASTQ paths come from
    the QC manifest and --input-dir is not needed (or read). A requested
    --tumour-sample/--normal-sample that is absent from the manifest falls
    back to alphabetical order, with a [WARN] -- watch for that line if the
    wrong sample got treated as the tumour.
"""

import argparse
import gzip
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import threading
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MIN_PLAUSIBLE_OUTPUT_BYTES = 1024  # A valid BAM is at least this size.

# Known reference genomes. Extend if you need GRCh37/hg19.
# Source bucket: gs://gcp-public-data--broad-references (public, no login).
# NOTE: the older https://s3.amazonaws.com/broad-references/... path that
# used to be here is dead (404), and that bucket has no .fasta.gz at all.
KNOWN_REFERENCES = {
    "hg38": {
        "url": "https://storage.googleapis.com/gcp-public-data--broad-references/"
               "hg38/v0/Homo_sapiens_assembly38.fasta",
        # Where install_pipeline.py puts this genome, relative to the
        # reference cache, and the name it gives it.
        "subdir": "hg38",
        "filename": "Homo_sapiens_assembly38.fasta",
    },
}

# Shared cache for downloaded genomes, deliberately outside the output
# directory: the reference and its indices are identical for every run and
# cost ~20 GB and 1-2 hours of bwa-mem2 indexing, so a per-run copy bought
# that wait again for every new -o. Mirrors
# comprehensive_variant_calling.py -- keep the two in step.
DEFAULT_REFERENCE_DIR = os.environ.get(
    "PIPELINE_REFERENCE_DIR", os.path.expanduser("~/data/references"))


# ---------------------------------------------------------------------------
# Command execution
# ---------------------------------------------------------------------------
def run_command(cmd, tag, log_path=None, dry_run=False):
    """
    Run a shell command and stream its output to console and log file.

    KEY DESIGN DECISION:
      We use subprocess.Popen() with line-by-line reading instead of
      subprocess.run(). Why? Because BWA-MEM2 can run for hours on large
      whole-genome samples. subprocess.run() buffers ALL output until the
      process finishes -- you'd see nothing for hours. Popen gives us
      real-time progress so you know it's not stuck.

    Returns:
        int: Exit code (0 = success, non-zero = failure).
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

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Merge stderr into stdout.
            text=True,
            bufsize=1,  # Line-buffered for real-time output.
        )

        # Print each line as it arrives. The "    | " prefix makes it
        # easy to distinguish tool output from our own messages.
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
        # This happens when the tool (bwa-mem2, samtools) is not installed.
        print(f"[{tag}] ERROR: '{cmd[0]}' not found on PATH.")
        print("        Install it or check your PATH environment variable.")
        return 127
    except OSError as exc:
        print(f"[{tag}] ERROR: {exc}")
        return 1
    finally:
        if log_fh:
            log_fh.close()


def tool_version(executable):
    """
    Capture a tool's version string for the run manifest.

    We try --version, -v, and -V because different tools use different
    flags. Returns the first line of output from whichever flag works.
    """
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
    Compute SHA-256 of the first `max_bytes` of a file.

    We hash only the first 64 MB (default) for speed. This is enough to
    detect accidentally swapped or corrupted inputs. Set max_bytes=None
    to hash the entire file (slow for large FASTQs).
    """
    h = hashlib.sha256()
    bytes_read = 0
    with open(path, "rb") as fh:
        while True:
            chunk_size = 1024 * 1024  # 1 MB chunks.
            if max_bytes is not None:
                remaining = max_bytes - bytes_read
                if remaining <= 0:
                    break
                chunk_size = min(chunk_size, remaining)
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
            bytes_read += len(chunk)
    return {
        "sha256": h.hexdigest(),
        "bytes_hashed": bytes_read,
        "partial": max_bytes is not None and bytes_read >= max_bytes,
    }


# ---------------------------------------------------------------------------
# Reference genome handling
# ---------------------------------------------------------------------------
def is_gzipped(path):
    """Return True if `path` begins with the gzip magic bytes (1f 8b)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def remote_content_length(url, timeout=60):
    """
    Size the server reports for `url`, or None when it cannot be learned.

    None always means "unknown", never "mismatch". A machine with no
    network must still be able to use a reference it already holds, so
    every caller treats None as "cannot check" and proceeds.

    ARGS:
        url:     URL to size.
        timeout: Seconds allowed for the request.

    RETURNS:
        int or None: Size in bytes, if the server reported one.
    """
    if shutil.which("curl"):
        cmd = ["curl", "-sIL", "--max-time", str(timeout), url]
    elif shutil.which("wget"):
        cmd = ["wget", "--spider", "-S", "--timeout", str(timeout), url]
    else:
        return None
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True,
                             timeout=timeout + 30)
    except Exception:
        return None
    # Take the LAST Content-Length: with redirects there is one per hop,
    # and only the final hop describes the file we would actually get.
    size = None
    for line in res.stdout.splitlines():
        if line.lower().strip().startswith("content-length:"):
            try:
                size = int(line.split(":", 1)[1].strip())
            except ValueError:
                pass
    return size


def download_verified(url, dest, dry_run=False):
    """
    Download `url` to `dest` so that `dest` exists ONLY if it is complete.

    THE BUG THIS PREVENTS
      Downloading straight to the final path leaves a truncated file
      behind whenever the transfer is interrupted -- a dropped link, a
      full disk, Ctrl-C. The cache check is only "does this path exist",
      so every later run reuses that stump quite happily. A 3 GB
      reference that stopped at 2 GB yields a genome ending partway
      through a chromosome, and everything downstream aligns and calls
      against it without complaint.

      Writing to a temporary name and renaming only after the size is
      confirmed makes the file's existence mean what the cache check has
      always assumed it means. os.replace() is atomic within a
      filesystem, so an interrupted run leaves the .partial behind and
      never a half-written `dest`.

    ARGS:
        url:     Source URL.
        dest:    Final path. Created only on success.
        dry_run: If True, print the command and change nothing.

    RETURNS:
        bool: True on success (or dry run), False otherwise.
    """
    partial = dest + ".partial"
    if shutil.which("wget"):
        cmd = ["wget", "-q", "-O", partial, url]
    elif shutil.which("curl"):
        cmd = ["curl", "-sL", "-o", partial, url]
    else:
        print("[ERROR] Need wget or curl to download the reference.")
        return False

    # NOTE: dry_run MUST be forwarded. Without it a dry run silently
    # pulls ~3 GB from GCS.
    code = run_command(cmd, tag="Download reference", dry_run=dry_run)
    if dry_run:
        return True

    def _discard(reason):
        print(f"[ERROR] {reason}")
        try:
            if os.path.exists(partial):
                os.remove(partial)
        except OSError:
            print(f"[WARN] Could not remove {partial}; delete it by hand.")
        return False

    if code != 0:
        return _discard(f"Download failed from {url}")

    actual = os.path.getsize(partial) if os.path.exists(partial) else 0
    if actual == 0:
        return _discard(f"Download from {url} produced an empty file.")

    expected = remote_content_length(url)
    if expected is not None and actual != expected:
        return _discard(
            f"Download from {url} is {actual} bytes, server said "
            f"{expected}. Treating it as truncated.")

    os.replace(partial, dest)
    return True


def cached_download_is_complete(path, url):
    """
    Decide whether an already-downloaded file can be trusted.

    Existence is not evidence: see download_verified() above for how a
    truncated file gets there. Where the server will tell us the size we
    compare against it; where it will not -- offline, or a server that
    omits Content-Length -- we say so and accept the file, because
    refusing to run offline would be worse than the risk.

    ARGS:
        path: Local file to judge.
        url:  URL it was fetched from.

    RETURNS:
        bool: True if the file should be used as-is.
    """
    try:
        actual = os.path.getsize(path)
    except OSError:
        return False
    if actual == 0:
        print(f"[WARN] Cached reference {path} is empty.")
        return False

    expected = remote_content_length(url)
    if expected is None:
        print(f"[WARN] Could not reach {url} to verify the cached "
              f"reference; using {path} unchecked.")
        return True
    if actual != expected:
        print(f"[WARN] Cached reference {path} is {actual} bytes but the "
              f"server reports {expected} -- it is incomplete.")
        return False
    return True


def reference_index_files(fasta):
    """Every index path ensure_indices() would build for this FASTA."""
    return [fasta + ".fai", os.path.splitext(fasta)[0] + ".dict",
            fasta + ".bwt.2bit.64", fasta + ".amb", fasta + ".ann",
            fasta + ".pac", fasta + ".0123"]


def find_cached_reference(name, reference_dir):
    """
    An already-downloaded copy of a known genome, or None.

    Looks where install_pipeline.py puts it first, then at the flat layout
    this script uses when it downloads the genome itself. Only a plausible
    FASTA counts, so a half-finished download is re-fetched rather than fed
    to bwa-mem2.
    """
    entry = KNOWN_REFERENCES.get(name)
    if not entry or not reference_dir:
        return None

    reference_dir = os.path.expanduser(reference_dir)
    candidates = []
    if entry.get("subdir") and entry.get("filename"):
        candidates.append(os.path.join(reference_dir, entry["subdir"],
                                       entry["filename"]))
    if entry.get("filename"):
        candidates.append(os.path.join(reference_dir, entry["filename"]))
    candidates.append(os.path.join(reference_dir, name, f"{name}.fasta"))
    candidates.append(os.path.join(reference_dir, f"{name}.fasta"))

    for path in candidates:
        # 100 MB rules out an error page or an interrupted transfer without
        # hashing 3 GB on every run.
        if os.path.isfile(path) and os.path.getsize(path) > 100 * 1024 ** 2:
            missing = [os.path.basename(i)
                       for i in reference_index_files(path)
                       if not os.path.exists(i)]
            if missing:
                print(f"[INFO] Shared reference found without all indices "
                      f"({', '.join(missing)}); they will be built once, "
                      f"beside it.")
            return os.path.abspath(path)
    return None


def reference_cache_target(name, reference_dir, output_dir, dry_run=False):
    """
    Where a downloaded genome should be written: (directory, fasta path).

    The shared cache, unless it cannot be created or written; the fallback
    to output_dir/reference is announced, because a silent one would
    quietly reintroduce the per-run download this exists to stop.
    """
    entry = KNOWN_REFERENCES.get(name, {})
    filename = entry.get("filename") or f"{name}.fasta"
    if reference_dir:
        shared = os.path.join(os.path.expanduser(reference_dir),
                              entry.get("subdir") or name)
        if dry_run:
            return shared, os.path.join(shared, filename)
        try:
            os.makedirs(shared, exist_ok=True)
            if os.access(shared, os.W_OK):
                return shared, os.path.join(shared, filename)
            print(f"[WARN] Reference cache {shared} is not writable.")
        except OSError as exc:
            print(f"[WARN] Cannot use reference cache {shared}: {exc}")

    fallback = os.path.join(output_dir, "reference")
    print(f"[WARN] Falling back to {fallback}; this copy is not shared, so "
          f"the next output directory will download and index its own. "
          f"Set --reference-dir to somewhere writable to avoid that.")
    if not dry_run:
        os.makedirs(fallback, exist_ok=True)
    return fallback, os.path.join(fallback, filename)


def ensure_reference(args):
    """
    Verify or download the reference genome FASTA.

    If the user passes a known genome name like 'hg38', we download it
    from the Broad Institute's public S3 bucket. If they pass a file
    path, we verify it exists.

    DRY-RUN CONTRACT (important when debugging):
      Every expensive side effect here is gated on args.dry_run. The hg38
      download is ~3 GB and decompression writes ~3 GB more, so a "preview"
      must never perform them. In dry-run mode we print what WOULD happen
      and return the path the reference would occupy; ensure_indices() is
      dry-run aware too, so a non-existent path is fine downstream.

      This mirrors comprehensive_variant_calling.py's copy of this function
      -- keep the two in step if you change either.

    Returns:
        str: Absolute path to the reference FASTA.
    """
    ref = args.reference
    # Read the dry-run flag off the args namespace rather than adding a
    # parameter, so existing call sites stay unchanged.
    dry_run = getattr(args, "dry_run", False)

    # Check if this is a known genome name (not a file path).
    if ref in KNOWN_REFERENCES and not os.path.isfile(ref):
        entry = KNOWN_REFERENCES[ref]
        url = entry["url"]

        # The shared cache first. An installed genome is already indexed, so
        # finding it here saves the download AND the 1-2 hour bwa-mem2 index
        # that would otherwise follow.
        cached = find_cached_reference(ref, getattr(args, "reference_dir",
                                                    DEFAULT_REFERENCE_DIR))
        if cached:
            print(f"[INFO] Using the shared reference: {cached}")
            print("[INFO] Nothing to download; indices are reused as found.")
            ref = cached
        else:
            ref_dir, local_ref = reference_cache_target(
                ref, getattr(args, "reference_dir", DEFAULT_REFERENCE_DIR),
                args.output_dir, dry_run=dry_run)

            needs_download = True
            if os.path.exists(local_ref):
                if dry_run or cached_download_is_complete(local_ref, url):
                    print(f"[INFO] Using cached reference: {local_ref}")
                    needs_download = False
                else:
                    print("[INFO] Re-downloading the reference.")

            if needs_download:
                print(f"[INFO] Downloading reference genome '{ref}' into "
                      f"{ref_dir} ...")
                print("[INFO] This is a one-off: later runs reuse it from "
                      "there, whatever their output directory.")
                if not download_verified(url, local_ref, dry_run=dry_run):
                    sys.exit(1)

            ref = local_ref

            # Downstream tools need plain-text FASTA. The Broad GCS bucket
            # serves hg38 uncompressed, so this normally no-ops -- it is
            # kept because it still rescues a copy cached from the older
            # gzipped URL, and any future KNOWN_REFERENCES entry that is
            # compressed. is_gzipped() returns False for a missing file, so
            # after a dry-run "download" (which wrote nothing) this branch
            # is simply skipped.
            if is_gzipped(ref):
                if dry_run:
                    print(f"[DRY RUN] Would decompress reference {ref}")
                else:
                    print(f"[INFO] Decompressing reference {ref} ...")
                    tmp_path = ref + ".tmp"
                    try:
                        with gzip.open(ref, "rb") as src, \
                                open(tmp_path, "wb") as dst:
                            shutil.copyfileobj(src, dst,
                                               length=8 * 1024 * 1024)
                        os.replace(tmp_path, ref)
                    except OSError as exc:
                        print(f"[ERROR] Failed to decompress downloaded "
                              f"reference: {exc}")
                        sys.exit(1)

    # Resolve to absolute path and verify existence.
    ref = os.path.abspath(ref)
    if not os.path.isfile(ref):
        # In a dry run the reference legitimately may not exist yet (we just
        # skipped downloading it), so warn instead of aborting -- the same
        # philosophy as the pre-flight tool checks, which are also relaxed
        # for --dry-run.
        if dry_run:
            print(f"[DRY RUN] Reference not present yet: {ref}")
            print("          A real run would download/build it here.")
            return ref
        print(f"[ERROR] Reference FASTA not found: {ref}")
        print("        Provide a valid path or use 'hg38' for auto-download.")
        sys.exit(1)

    print(f"[INFO] Reference: {ref}")
    return ref


def ensure_indices(ref, threads, dry_run=False):
    """
    Build all required index files for the reference genome.

    BWA-MEM2 requires its own index (.amb, .ann, .bwt, .pac, .sa files,
    plus .bwt.2bit.64 for BWA-MEM2 specifically). samtools needs a .fai
    index. GATK needs a .dict sequence dictionary.

    We check for each index file and only build what's missing.
    """
    fai = ref + ".fai"
    dict_file = os.path.splitext(ref)[0] + ".dict"

    # --- BWA-MEM2 index ---
    # BWA-MEM2 uses a .bwt.2bit.64 file as its primary index marker.
    bwt_mem2 = ref + ".bwt.2bit.64"
    if not os.path.exists(bwt_mem2):
        print("[INFO] Building BWA-MEM2 index (this takes 1-2 hours for hg38) ...")
        code = run_command(
            ["bwa-mem2", "index", ref],
            tag="BWA-MEM2 index",
            dry_run=dry_run,
        )
        if code != 0:
            print("[ERROR] BWA-MEM2 indexing failed. Check your reference FASTA.")
            sys.exit(1)
    else:
        print("[INFO] BWA-MEM2 index already exists.")

    # --- samtools faidx ---
    if not os.path.exists(fai):
        print("[INFO] Building samtools .fai index ...")
        code = run_command(
            ["samtools", "faidx", ref],
            tag="samtools faidx",
            dry_run=dry_run,
        )
        if code != 0:
            print("[ERROR] samtools faidx failed.")
            sys.exit(1)
    else:
        print("[INFO] samtools .fai index already exists.")

    # --- GATK CreateSequenceDictionary ---
    if not os.path.exists(dict_file):
        print("[INFO] Building GATK sequence dictionary ...")
        # GATK must be on PATH for this. If it's not installed, skip
        # but warn -- downstream BQSR will need it.
        if shutil.which("gatk") is None:
            print("[WARN] GATK not found on PATH. Skipping .dict creation.")
            print("       You'll need to run 'gatk CreateSequenceDictionary' later.")
        else:
            code = run_command(
                ["gatk", "CreateSequenceDictionary",
                 "-R", ref, "--VERBOSITY", "ERROR"],
                tag="GATK dict",
                dry_run=dry_run,
            )
            if code != 0:
                print("[ERROR] GATK CreateSequenceDictionary failed.")
                sys.exit(1)
    else:
        print("[INFO] GATK .dict already exists.")

    print("[INFO] All reference indices ready.")


# ---------------------------------------------------------------------------
# Alignment (BWA-MEM2)
# ---------------------------------------------------------------------------
def align_sample(sample_name, r1, r2, ref, output_bam, threads,
                 two_pass, log_path, dry_run=False):
    """
    Align paired-end reads to the reference genome and produce a
    sorted, indexed BAM file.

    ALIGNMENT STRATEGY:
      We pipe: bwa-mem2 mem | samtools sort | samtools index

      This avoids writing an unsorted BAM to disk, which would:
        1. Require ~3x the disk space (unsorted + sorted + index).
        2. Take twice as long (write + read back for sorting).

    READ GROUP (@RG) HEADER:
      GATK and other downstream tools require a Read Group header.
      We set it here during alignment so it's baked into the BAM.
      Required fields:
        ID: Unique read group identifier (usually the sample name).
        SM: Sample name (used by GATK for variant calling).
        PL: Platform (ILLUMINA for NextSeq 1000).
        LB: Library prep identifier.
        PU: Platform unit (flowcell + lane + barcode).

    TWO-PASS ALIGNMENT:
      If --two-pass is enabled, BWA-MEM2 runs twice:
        Pass 1: Identify difficult-to-map regions.
        Pass 2: Re-align reads in those regions with tuned parameters.

      This is especially valuable for cancer samples because tumour
      genomes contain structural rearrangements that create novel
      junctions not present in the reference.

    Args:
        sample_name: Name used for the @RG SM tag (e.g. "TUMOUR_01").
        r1: Path to Read 1 FASTQ.
        r2: Path to Read 2 FASTQ.
        ref: Path to reference genome FASTA.
        output_bam: Path where the sorted BAM will be written.
        threads: Number of CPU threads to use.
        two_pass: Enable BWA-MEM2 two-pass alignment.
        log_path: Path to the log file for this sample.
        dry_run: If True, show commands without executing.

    Returns:
        str: Path to the output BAM on success, None on failure.
    """
    # --- Build the Read Group header string ---
    # GATK requires this in the format @RG\\tID:...\\tSM:...\\tPL:...\\tLB:...
    # The double-backslash is because we're building a string that will be
    # passed to the shell -- we need literal \t characters, not tab escapes.
    rg_id = sample_name
    rg_header = (
        f"@RG\\tID:{rg_id}\\tSM:{sample_name}\\tPL:ILLUMINA"
        f"\\tLB:{sample_name}_lib\\tPU:unknown"
    )

    # --- Build BWA-MEM2 command ---
    cmd = [
        "bwa-mem2", "mem",
        "-t", str(threads),
        "-R", rg_header,
    ]

    if two_pass:
        # Two-pass mode: BWA-MEM2 identifies difficult regions in pass 1
        # and re-aligns them in pass 2. This is slower but more sensitive.
        cmd.append("-2")

    cmd += [ref, r1, r2]

    # --- Print what we're about to do ---
    pipeline_desc = (
        f"{shlex.join(cmd)} | samtools sort | samtools index"
    )
    print(f"\n[Align {sample_name}] $ {pipeline_desc}", flush=True)

    if dry_run:
        print(f"[Align {sample_name}] (dry run - not executed)")
        return output_bam

    log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if log_fh:
            log_fh.write(f"\n### Align {sample_name}\n### {pipeline_desc}\n")

        # --- PIPELINE: bwa-mem2 | samtools sort ---
        # Step 1: Start BWA-MEM2. Its stdout (SAM output) is piped to
        #         samtools sort's stdin.
        bwa_proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Drain BWA's stderr from a background thread. BWA streams
        # progressive [M::...] status lines to stderr; reading them only
        # after the process exits lets a long alignment fill the OS pipe
        # buffer (~64 KB), block BWA, and deadlock the pipeline.
        bwa_stderr_lines = []
        stderr_lock = threading.Lock()

        def _drain_bwa_stderr():
            for line in bwa_proc.stderr:
                with stderr_lock:
                    bwa_stderr_lines.append(line.rstrip("\n"))

        stderr_thread = threading.Thread(target=_drain_bwa_stderr, daemon=True)
        stderr_thread.start()

        # Step 2: samtools sort reads from stdin and writes the sorted BAM.
        # We use -(threads - 1) for sort because 1 thread is used by bwa.
        sort_proc = subprocess.Popen(
            ["samtools", "sort",
             "-@", str(max(1, threads - 1)),
             "-o", output_bam,
             "-"],  # "-" means read from stdin.
            stdin=bwa_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        # Close bwa's stdout in the parent process so that if bwa exits,
        # sort will get a broken pipe signal and terminate cleanly.
        bwa_proc.stdout.close()

        # --- Stream sort output (progress messages) ---
        for line in sort_proc.stdout:
            line = line.rstrip("\n")
            print(f"    | {line}", flush=True)
            if log_fh:
                log_fh.write(line + "\n")

        sort_proc.wait()
        bwa_proc.wait()

        # --- Collect BWA-MEM2 stderr (alignment statistics) ---
        # BWA writes useful stats to stderr: read counts, mapping rates, etc.
        stderr_thread.join()
        if bwa_stderr_lines:
            print("    | [bwa-mem2 stderr]")
            for stderr_line in bwa_stderr_lines:
                print(f"    |   {stderr_line}")
            if log_fh:
                log_fh.write("[bwa-mem2 stderr]\n"
                             + "\n".join(bwa_stderr_lines) + "\n")

        # --- Check exit codes ---
        if bwa_proc.returncode != 0:
            print(f"[Align {sample_name}] ERROR: bwa-mem2 failed "
                  f"(exit code {bwa_proc.returncode}).")
            print("  Possible causes:")
            print("    - Reference genome not indexed (run bwa-mem2 index).")
            print("    - FASTQ files corrupted or empty.")
            print("    - Insufficient memory (BWA-MEM2 needs ~32GB for hg38).")
            return None

        if sort_proc.returncode != 0:
            print(f"[Align {sample_name}] ERROR: samtools sort failed "
                  f"(exit code {sort_proc.returncode}).")
            return None

        # --- Index the sorted BAM ---
        # A .bai index is required for random access (e.g. viewing regions
        # in IGV, or for GATK tools that process specific intervals).
        print(f"\n[Index {sample_name}] Creating BAM index ...")
        index_code = run_command(
            ["samtools", "index", output_bam],
            tag=f"Index {sample_name}",
            log_path=log_path,
        )
        if index_code != 0:
            print(f"[Align {sample_name}] ERROR: samtools index failed.")
            return None

        return output_bam

    except FileNotFoundError as exc:
        print(f"[Align {sample_name}] ERROR: {exc}")
        return None
    finally:
        if log_fh:
            log_fh.close()


# ---------------------------------------------------------------------------
# Alignment statistics
# ---------------------------------------------------------------------------
def collect_stats(bam_path, sample_name, stats_dir, log_path, dry_run=False):
    """
    Generate alignment statistics for the BAM file.

    We run three samtools commands:
      1. flagstat: Count reads by flag (mapped, unmapped, duplicate, etc.)
      2. idxstats: Per-chromosome alignment counts.
      3. stats: Detailed alignment statistics (identity, insert size, etc.)

    These metrics are essential for QC:
      - Mapping rate should be >95% for tumour DNA.
      - Duplicate rate should be <30% (higher suggests PCR over-amplification).
      - Insert size peak should match library preparation (typically 200-400bp).

    Returns:
        dict: Alignment statistics.
    """
    if dry_run:
        return {"dry_run": True}

    stats = {}
    os.makedirs(stats_dir, exist_ok=True)

    # --- samtools flagstat ---
    # Produces a human-readable summary like:
    #   12345678 + 0 in total (QC-passed reads + QC-failed reads)
    #   12345678 + 0 mapped (99.50% : 0.00%)
    flagstat_path = os.path.join(stats_dir, f"{sample_name}.flagstat.txt")
    cmd = ["samtools", "flagstat", bam_path]
    try:
        res = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(flagstat_path, "w") as fh:
                fh.write(res.stdout)
            # Parse key metrics from flagstat output.
            for line in res.stdout.splitlines():
                if "in total" in line:
                    stats["total_reads"] = line.split("+")[0].strip()
                elif "mapped (" in line and "primary" not in line:
                    stats["mapped_reads"] = line.split("+")[0].strip()
                    # Extract percentage: "mapped (99.50% : 0.00%)"
                    pct = line.split("(")[1].split("%")[0]
                    stats["mapping_rate"] = pct + "%"
                elif "duplicates" in line and "primary" not in line:
                    stats["duplicates"] = line.split("+")[0].strip()
    except Exception as exc:
        print(f"[WARN] flagstat failed: {exc}")

    # --- samtools idxstats ---
    # Per-chromosome read counts. Useful for checking coverage uniformity.
    idxstats_path = os.path.join(stats_dir, f"{sample_name}.idxstats.txt")
    cmd = ["samtools", "idxstats", bam_path]
    try:
        res = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(idxstats_path, "w") as fh:
                fh.write(res.stdout)
    except Exception as exc:
        print(f"[WARN] idxstats failed: {exc}")

    # --- samtools stats ---
    # Detailed statistics: insert size distribution, alignment identity, etc.
    stats_path = os.path.join(stats_dir, f"{sample_name}.stats.txt")
    cmd = ["samtools", "stats", bam_path]
    try:
        res = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(stats_path, "w") as fh:
                fh.write(res.stdout)
            # Extract insert size peak from IS line.
            for line in res.stdout.splitlines():
                if line.startswith("IS") and "insert size" not in line:
                    parts = line.split("\t")
                    if len(parts) >= 7:
                        # IS lines: IS\tinsert_size\tcount_fq\tcount_rc\t...
                        try:
                            size = int(parts[1])
                            count = int(parts[2])
                            if count > stats.get("_max_is_count", 0):
                                stats["insert_size_peak"] = size
                                stats["_max_is_count"] = count
                        except ValueError:
                            pass
    except Exception as exc:
        print(f"[WARN] stats failed: {exc}")

    # Clean up internal tracking key.
    stats.pop("_max_is_count", None)

    return stats


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def build_parser():
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Align paired-end FASTQ files to the human reference "
                    "genome using BWA-MEM2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Input / Output ---
    io = parser.add_argument_group("input / output")
    # NOT required=True: in --manifest mode the manifest carries the full
    # cleaned-FASTQ paths, so an input dir is never consulted. Requiring it
    # here would reject a perfectly self-sufficient manifest-driven command
    # (comprehensive_variant_calling.py treats it the same way). The
    # requirement for the direct-input mode is enforced in validate_args().
    io.add_argument("-i", "--input-dir",
                    help="Directory containing cleaned FASTQ files. Required "
                         "unless --manifest is given.")
    io.add_argument("-o", "--output-dir", required=True,
                    help="Directory for aligned BAM files and statistics.")
    io.add_argument("-r", "--reference", required=True,
                    help="Path to reference genome FASTA, or a known name "
                         "like 'hg38' to download automatically.")
    io.add_argument("--reference-dir", default=DEFAULT_REFERENCE_DIR,
                    help="Shared directory holding downloaded genomes and "
                         "their indices, used when --reference names a "
                         "genome rather than a FASTA. Kept out of the "
                         "output directory so every run reuses one copy "
                         "instead of downloading and indexing its own "
                         f"(default: {DEFAULT_REFERENCE_DIR}, or "
                         "$PIPELINE_REFERENCE_DIR).")
    io.add_argument("--manifest", default=None,
                    help="Path to a fastq_qc_clean.py run manifest JSON. "
                         "When provided, tumour/normal sample names and "
                         "cleaned FASTQ paths are read from it (any explicit "
                         "--tumour-r1/r2 etc. are overridden). This script "
                         "then produces a BAM layout compatible with the "
                         "variant-calling scripts.")

    # --- Sample specification ---
    # NOTE: these are NOT required=True because, when --manifest is given,
    # sample names and FASTQ paths come from the manifest instead. They are
    # required only in the direct-input mode (enforced in validate_args).
    sample = parser.add_argument_group("sample specification")
    sample.add_argument("--tumour-sample", default=None,
                        help="Sample name for the tumour (used in BAM headers).")
    sample.add_argument("--tumour-r1", default=None,
                        help="Read 1 FASTQ for the tumour (filename or path).")
    sample.add_argument("--tumour-r2", default=None,
                        help="Read 2 FASTQ for the tumour.")
    sample.add_argument("--normal-sample", default=None,
                        help="Sample name for the matched normal.")
    sample.add_argument("--normal-r1", default=None,
                        help="Read 1 FASTQ for the matched normal.")
    sample.add_argument("--normal-r2", default=None,
                        help="Read 2 FASTQ for the matched normal.")

    # --- Alignment options ---
    align = parser.add_argument_group("alignment options")
    align.add_argument("--two-pass", action="store_true",
                       help="Enable BWA-MEM2 two-pass alignment. Slower but "
                            "more sensitive for difficult regions (recommended "
                            "for cancer samples with structural rearrangements).")
    align.add_argument("--threads", type=int, default=8,
                       help="Number of CPU threads. BWA-MEM2 uses all threads "
                            "for alignment; samtools sort uses (threads - 1).")
    align.add_argument("--skip-stats", action="store_true",
                       help="Skip alignment statistics generation.")
    align.add_argument("--skip-fastqc", action="store_true",
                       help="Skip FastQC on the aligned BAM.")
    align.add_argument("--overwrite", action="store_true",
                       help="Reprocess a sample whose sorted+indexed BAM "
                            "already exists (default: skip it). Matches "
                            "fastq_qc_clean.py's --overwrite behaviour.")
    align.add_argument("--dry-run", action="store_true",
                       help="Show commands without executing them.")
    return parser


def validate_args(args, parser):
    """Validate argument combinations that argparse can't catch."""
    if args.threads < 1:
        parser.error("--threads must be >= 1.")

    # If a manifest is given, sample info is read from it, so explicit
    # sample args are optional.
    if args.manifest:
        if not os.path.isfile(args.manifest):
            parser.error(f"Manifest not found: {args.manifest}")
        return  # nothing else to validate

    # Direct-input mode: --input-dir is needed both as the fallback location
    # for bare FASTQ filenames (resolved at the bottom of this function) and
    # as documentation of where the reads came from.
    if not args.input_dir:
        parser.error("--input-dir is required (unless --manifest is given).")

    # Without a manifest, the tumour args are REQUIRED.
    if not (args.tumour_sample and args.tumour_r1 and args.tumour_r2):
        parser.error("Please provide --tumour-sample, --tumour-r1 and "
                     "--tumour-r2 (or use --manifest to reuse a "
                     "fastq_qc_clean.py run).")

    # If normal is specified, all normal args must be present.
    if bool(args.normal_sample) != bool(args.normal_r1 and args.normal_r2):
        parser.error("If --normal-sample is provided, you must also provide "
                     "--normal-r1 and --normal-r2 (and vice versa).")

    # Resolve FASTQ paths: try the input directory if the path doesn't exist.
    for role in ("tumour_r1", "tumour_r2", "normal_r1", "normal_r2"):
        val = getattr(args, role)
        if val is not None and not os.path.isfile(val):
            candidate = os.path.join(args.input_dir, os.path.basename(val))
            if os.path.isfile(candidate):
                setattr(args, role, candidate)
            else:
                parser.error(f"FASTQ not found: {val} (also tried {candidate})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def load_qc_manifest(manifest_path):
    """
    Load a fastq_qc_clean.py run manifest and extract cleaned FASTQ paths.

    This mirrors the function in comprehensive_variant_calling.py so the
    two scripts behave identically when chained:
        fastq_qc_clean.py  ->  align_reads.py  /  comprehensive_variant_calling.py

    NOTE: fastq_qc_clean.py omits the "outputs" key for records it skipped
    ("skipped_existing") even though the cleaned files DO exist, at the
    standard location <qc_output>/cleaned_fastq/{sample}_R1.clean.fastq.gz.
    We reconstruct that path from the manifest's output_dir in that case.

    RETURNS:
        dict: {sample_name: {"r1": path, "r2": path, "lanes_merged": bool}}
    """
    with open(manifest_path, encoding="utf-8") as fh:
        data = json.load(fh)

    qc_output_dir = data.get("output_dir") or os.path.dirname(manifest_path)
    samples = {}
    for record in data.get("samples", []):
        if record.get("status") == "failed":
            print(f"[WARN] Sample '{record.get('sample')}' failed in QC; "
                  f"excluding from downstream.")
            continue
        name = record.get("sample")
        outputs = record.get("outputs", {})
        if outputs.get("r1") and outputs.get("r2"):
            r1, r2 = outputs["r1"], outputs["r2"]
        else:
            r1 = os.path.join(qc_output_dir, "cleaned_fastq",
                              f"{name}_R1.clean.fastq.gz")
            r2 = os.path.join(qc_output_dir, "cleaned_fastq",
                              f"{name}_R2.clean.fastq.gz")
            if not (os.path.isfile(r1) and os.path.isfile(r2)):
                # The QC output may have been moved since it ran; the
                # cleaned reads normally live next to the manifest itself.
                alt_dir = os.path.join(os.path.dirname(manifest_path),
                                       "cleaned_fastq")
                alt_r1 = os.path.join(alt_dir, f"{name}_R1.clean.fastq.gz")
                alt_r2 = os.path.join(alt_dir, f"{name}_R2.clean.fastq.gz")
                if os.path.isfile(alt_r1) and os.path.isfile(alt_r2):
                    r1, r2 = alt_r1, alt_r2
                else:
                    print(f"[WARN] Sample '{name}' has no cleaned outputs in "
                          f"manifest and none found at {qc_output_dir}/cleaned_fastq; "
                          f"excluding.")
                    continue
        samples[name] = {
            "r1": r1,
            "r2": r2,
            "lanes_merged": record.get("lanes_merged", False),
        }
    return samples


def main():
    """Main entry point: parse arguments, run alignment pipeline."""
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    started = datetime.now(timezone.utc)

    # ---- Pre-flight checks ----
    # Verify required tools are installed. Not enforced for --dry-run: a
    # dry run should preview commands even on machines without the tools.
    required_tools = {"bwa-mem2": "bwa-mem2", "samtools": "samtools"}
    if not args.dry_run:
        missing = [label for label, exe in required_tools.items()
                   if shutil.which(exe) is None]
        if missing:
            print(f"[ERROR] Required tools not found on PATH: "
                  f"{', '.join(missing)}")
            sys.exit(1)

    # ---- Create output directories ----
    output_dir = os.path.abspath(args.output_dir)
    # NOTE: the BAM subdirectory is named "aligned" to stay compatible with
    # the directory layout expected by comprehensive_variant_calling.py
    # (which uses {output}/aligned/*.sorted.bam).
    bam_dir = os.path.join(output_dir, "aligned")
    stats_dir = os.path.join(output_dir, "stats")
    logs_dir = os.path.join(output_dir, "logs")
    qc_dir = os.path.join(output_dir, "fastqc")

    for d in (bam_dir, stats_dir, logs_dir):
        os.makedirs(d, exist_ok=True)
    if not args.skip_fastqc:
        os.makedirs(qc_dir, exist_ok=True)

    # ---- Reference genome ----
    ref = ensure_reference(args)
    ensure_indices(ref, args.threads, dry_run=args.dry_run)

    # ---- Sample setup ----
    # If a QC manifest is provided, sample names and cleaned FASTQ paths
    # are read from it. This is the automated handoff from fastq_qc_clean.py.
    if args.manifest:
        manifest_samples = load_qc_manifest(args.manifest)
        if not manifest_samples:
            print("[ERROR] No usable samples found in manifest.")
            sys.exit(1)

        # Pick the tumour entry: user-specified name or the first sample.
        if args.tumour_sample and args.tumour_sample in manifest_samples:
            tumour = {"name": args.tumour_sample,
                      "r1": manifest_samples[args.tumour_sample]["r1"],
                      "r2": manifest_samples[args.tumour_sample]["r2"]}
        else:
            first_key = sorted(manifest_samples.keys())[0]
            if args.tumour_sample:
                print(f"[WARN] Requested tumour sample "
                      f"'{args.tumour_sample}' not in manifest; using first "
                      f"entry '{first_key}' as the tumour.")
            tumour = {"name": first_key,
                      "r1": manifest_samples[first_key]["r1"],
                      "r2": manifest_samples[first_key]["r2"]}

        normal = None
        if args.normal_sample and args.normal_sample in manifest_samples:
            normal = {"name": args.normal_sample,
                      "r1": manifest_samples[args.normal_sample]["r1"],
                      "r2": manifest_samples[args.normal_sample]["r2"]}
        else:
            if args.normal_sample:
                print(f"[WARN] Requested normal sample "
                      f"'{args.normal_sample}' not in manifest "
                      f"(--normal is auto-picked from the remaining samples).")
            if len(manifest_samples) >= 2:
                candidates = [k for k in sorted(manifest_samples.keys())
                              if k != tumour["name"]]
                if candidates:
                    nkey = candidates[0]
                    normal = {"name": nkey,
                              "r1": manifest_samples[nkey]["r1"],
                              "r2": manifest_samples[nkey]["r2"]}

        print(f"[INFO] Loaded cleaned sample(s) from manifest: "
              f"{args.manifest}")
    else:
        tumour = {
            "name": args.tumour_sample,
            "r1": os.path.abspath(args.tumour_r1),
            "r2": os.path.abspath(args.tumour_r2),
        }
        normal = None
        if args.normal_sample:
            normal = {
                "name": args.normal_sample,
                "r1": os.path.abspath(args.normal_r1),
                "r2": os.path.abspath(args.normal_r2),
            }

    print(f"\n{'=' * 72}")
    print("ALIGNMENT PIPELINE")
    print(f"{'=' * 72}")
    print(f"Tumour sample : {tumour['name']}")
    if normal:
        print(f"Normal sample : {normal['name']}")
    else:
        print("Mode          : tumour-only (no matched normal)")
    print(f"Reference     : {ref}")
    print(f"Two-pass      : {'yes' if args.two_pass else 'no'}")
    print(f"Threads       : {args.threads}")
    print(f"Output        : {output_dir}")
    print(f"{'=' * 72}")

    # ---- Build run manifest ----
    manifest = {
        "script": os.path.basename(__file__),
        "started_utc": started.isoformat(),
        "host": os.uname().nodename if hasattr(os, "uname") else None,
        "python": sys.version.split()[0],
        "command_line": shlex.join(sys.argv),
        "parameters": vars(args),
        "tool_versions": {name: tool_version(exe)
                          for name, exe in required_tools.items()},
        "reference": ref,
        "tumour": tumour["name"],
        "normal": normal["name"] if normal else None,
        "samples": [],
    }

    # ---- Alignment loop ----
    # We process each sample (tumour + optional normal) through the
    # same alignment pipeline.
    samples_to_align = [("tumour", tumour)]
    if normal:
        samples_to_align.append(("normal", normal))

    succeeded = []
    failed = []
    skipped = []

    for role, sample in samples_to_align:
        sample_name = sample["name"]
        log_path = os.path.join(logs_dir, f"{sample_name}.log")

        print(f"\n{'#' * 72}")
        print(f"# Aligning {role}: {sample_name}")
        print(f"{'#' * 72}")

        bam_path = os.path.join(bam_dir, f"{sample_name}.sorted.bam")
        record = {
            "role": role,
            "sample": sample_name,
            "r1": sample["r1"],
            "r2": sample["r2"],
            "output_bam": bam_path,
            "log": log_path,
        }

        if (not args.overwrite and os.path.exists(bam_path)
                and os.path.exists(bam_path + ".bai")):
            print(f"[SKIP] {sample_name}: aligned BAM already exists. Use "
                  f"--overwrite to realign.")
            record["status"] = "skipped_existing"
            manifest["samples"].append(record)
            skipped.append(sample_name)
            continue

        # --- Verify inputs exist and record provenance (skipped on a dry
        # run: there is nothing to hash, and this should preview commands
        # even when the input files don't exist yet) ---
        if not args.dry_run:
            missing = [f for f in (sample["r1"], sample["r2"])
                      if not os.path.exists(f)]
            if missing:
                print(f"[ERROR] FASTQ not found: {missing[0]}")
                record["status"] = "failed"
                record["error"] = f"FASTQ not found: {missing[0]}"
                manifest["samples"].append(record)
                failed.append(sample_name)
                continue

            record["input_provenance"] = []
            for read_file in (sample["r1"], sample["r2"]):
                entry = {
                    "path": read_file,
                    "size_bytes": os.path.getsize(read_file),
                }
                entry.update(file_digest(read_file))
                record["input_provenance"].append(entry)

        # --- Run alignment (align_sample() itself no-ops and just prints
        # the commands when dry_run is set, so this always runs) ---
        result = align_sample(
            sample_name,
            sample["r1"],
            sample["r2"],
            ref,
            bam_path,
            args.threads,
            args.two_pass,
            log_path,
            dry_run=args.dry_run,
        )

        if result is None:
            record["status"] = "failed"
            record["error"] = "Alignment failed. Check logs."
            manifest["samples"].append(record)
            failed.append(sample_name)
            continue

        # --- Verify output BAM exists and is plausible ---
        if not args.dry_run:
            if os.path.exists(bam_path):
                bam_size = os.path.getsize(bam_path)
                if bam_size < MIN_PLAUSIBLE_OUTPUT_BYTES:
                    print(f"[ERROR] Output BAM suspiciously small: "
                          f"{bam_size} bytes")
                    record["status"] = "failed"
                    record["error"] = f"BAM too small: {bam_size} bytes"
                    manifest["samples"].append(record)
                    failed.append(sample_name)
                    continue
                record["bam_size_bytes"] = bam_size
            else:
                record["status"] = "failed"
                record["error"] = "Output BAM not found after alignment."
                manifest["samples"].append(record)
                failed.append(sample_name)
                continue

        # --- Collect alignment statistics ---
        if not args.skip_stats:
            print(f"\n[Stats {sample_name}] Collecting alignment stats ...")
            alignment_stats = collect_stats(
                bam_path, sample_name, stats_dir, log_path,
                dry_run=args.dry_run,
            )
            record["stats"] = alignment_stats
            print(f"[Stats {sample_name}] Mapping rate: "
                  f"{alignment_stats.get('mapping_rate', 'N/A')}, "
                  f"Insert size peak: "
                  f"{alignment_stats.get('insert_size_peak', 'N/A')}")

        # --- FastQC on BAM (optional) ---
        if not args.skip_fastqc and shutil.which("fastqc"):
            run_command(
                ["fastqc", bam_path, "-o", qc_dir,
                 "--threads", "2", "--quiet"],
                tag=f"FastQC {sample_name}",
                log_path=log_path,
                dry_run=args.dry_run,
            )

        record["status"] = "ok"
        manifest["samples"].append(record)
        succeeded.append(sample_name)

    # ---- Write manifest ----
    finished = datetime.now(timezone.utc)
    manifest["finished_utc"] = finished.isoformat()
    manifest["duration_seconds"] = round(
        (finished - started).total_seconds(), 1
    )
    manifest["result"] = {
        "succeeded": succeeded,
        "failed": failed,
        "skipped": skipped,
    }

    manifest_path = os.path.join(
        output_dir,
        f"alignment_manifest_{started.strftime('%Y%m%dT%H%M%SZ')}.json",
    )
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    # ---- Print summary ----
    print(f"\n{'=' * 72}")
    print("ALIGNMENT COMPLETE")
    print(f"{'=' * 72}")
    print(f"Duration   : {manifest['duration_seconds']}s")
    print(f"Succeeded  : {', '.join(succeeded) or 'none'}")
    print(f"Skipped    : {', '.join(skipped) or 'none'}")
    print(f"Failed     : {', '.join(failed) or 'none'}")
    print(f"BAM files  : {bam_dir}")
    print(f"Statistics : {stats_dir}")
    print(f"Logs       : {logs_dir}")
    print(f"Manifest   : {manifest_path}")
    if not args.skip_fastqc:
        print(f"FastQC     : {qc_dir}")
    print(f"{'=' * 72}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
