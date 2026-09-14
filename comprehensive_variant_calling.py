#!/usr/bin/env python3
"""
comprehensive_variant_calling.py
================================
COMPREHENSIVE somatic variant calling pipeline for cancer DNA sequencing.

This is a COMPLETE, SELF-CONTAINED pipeline that takes raw Illumina
NextSeq 1000 paired-end FASTQ files and produces annotated somatic
variant calls -- no other scripts needed.

PIPELINE STEPS
--------------
  STEP 1: QC & Cleaning (fastp)
      - Adapter trimming, poly-G removal (NextSeq 1000 specific),
        quality filtering, UMI extraction (optional).

  STEP 2: Alignment (BWA-MEM2)
      - Two-pass alignment to hg38 reference genome.
      - Pipes directly into samtools sort (no intermediate files).

  STEP 3: Duplicate Marking (GATK MarkDuplicatesSpark)
      - Identifies PCR/optical duplicates from library amplification.

  STEP 4: Base Quality Score Recalibration (GATK BQSR)
      - Corrects systematic errors in base quality scores.

  STEP 5: Somatic Variant Calling (GATK Mutect2)
      - Calls SNVs and small indels in tumour vs. normal (or tumour-only).
      - Produces F1R2 tarballs for strand-bias modelling.

  STEP 6: Strand Bias Modelling (GATK LearnReadOrientationModel)
      - Critical for FFPE samples where oxidation artefacts are common.

  STEP 7: Contamination Estimation (GATK GetPileupSummaries +
          CalculateContamination)
      - Estimates the fraction of foreign DNA in the library from the
        allele balance at common germline SNPs. Runs only when
        --contamination-resource is given.
      - Without it FilterMutectCalls runs with --contamination-estimate
        0.0, so the "contamination" filter is inert: it still tags a
        handful of calls, but it is testing against an ASSUMED zero
        rather than a measurement, and cannot detect real contamination.

  STEP 8: Microsatellite Instability (MSIsensor2)
      - Scores MSI from the tumour BAM alone. Runs only when
        --msi-models is given.
      - PCGR will not do this: it restricts MSI to WGS/WES tumour-control
        runs and silently omits the section otherwise.

  STEP 9: Variant Filtering (GATK FilterMutectCalls)
      - Applies the filters recommended by the Mutect2 team, using the
        read-orientation model from STEP 6 and, when available, the
        contamination tables from STEP 7.
      - Then applies the --min-depth floor, which GATK has no equivalent
        of: FilterMutectCalls will PASS a call standing on two reads.

  STEP 10: COSMIC Annotation
      - Overlaps variants with COSMIC (Catalogue Of Somatic Mutations
        In Cancer) to identify known cancer driver mutations.

  STEP 11: SnpEff Annotation
      - Adds gene names, consequence types (missense, nonsense, etc.),
        and predicted impact (HIGH/MODERATE/LOW/MODIFIER).

  STEP 12: Target Coverage Check
      - Optional. Measures every region of --coverage-bed against the depth
        a call must reach and names the stretches that fall short, per
        sample, as HTML/JSON/TSV under <output-dir>/coverage.
      - It answers what the VCF cannot: a region nobody sequenced and a
        region that is wild type both yield no variant, so without this a
        capture dropout reads as a negative result all the way into the
        clinical report.
      - Runs only when --coverage-bed is given; skipped otherwise, and the
        run says so. Implemented in the sibling module coverage_report.py.

  STEP 13: Clinical Interpretation Report (PCGR)
      - Optional. Therapeutic actionability tiered by AMP/ASCO/CAP, as a
        self-contained HTML report.
      - TMB, MSI and COSMIC mutational-signature fitting (SBS3/HRD, MMR,
        APOBEC) are each OFF unless asked for: --pcgr-estimate-tmb,
        --pcgr-estimate-msi, --pcgr-estimate-signatures. PCGR does not
        warn when they are off, it just omits those sections.
      - Runs only when --pcgr-refdata-dir is given; skipped otherwise.
        Implemented in the sibling module pcgr_report.py.

  STEP 14: Summary Statistics
      - VCF stats, variant counts by type, COSMIC overlap counts.

COSMIC DATABASE
---------------
  COSMIC (https://cancer.sanger.ac.uk/cosmic) is the world's largest
  database of somatic mutations in human cancer. It contains:
    - Over 20 million coding mutations across millions of tumour samples.
    - Information on mutation frequency, sample count, and mutation type.
    - Gene-level summaries (e.g. TP53 is mutated in ~30% of cancers).

  This pipeline uses COSMIC in two ways:
    1. During Mutect2 calling: known COSMIC sites help improve sensitivity
       for low-VAF somatic variants.
    2. Post-call annotation: each called variant is annotated with COSMIC
       mutation IDs (COSM12345), sample counts, and mutation frequency.

  To obtain COSMIC data:
    1. Register at https://cancer.sanger.ac.uk/cosmic (free for academics).
    2. Download the "COSMIC Coding Mutation VCF" for your genome build.
    3. Pass the path via --cosmic /path/to/COSMIC_v98.vcf.gz

Requirements
------------
  - Python 3.8+
  - fastp       (QC and adapter trimming)
  - BWA-MEM2    (alignment)
  - samtools    (BAM manipulation)
  - GATK 4+     (variant calling, BQSR, Mutect2)
  - (optional) bcftools    (COSMIC annotation)
  - (optional) SnpEff      (variant annotation)
  - (optional) FastQC      (QC reports)

Usage Examples
--------------
  # Full pipeline with COSMIC:
  python comprehensive_variant_calling.py \\
      --input-dir raw_fastqs/ \\
      --output-dir results/ \\
      --reference hg38 \\
      --tumour-sample TUMOUR_01 \\
      --tumour-r1 T1_R1_001.fastq.gz \\
      --tumour-r2 T1_R2_001.fastq.gz \\
      --normal-sample NORMAL_01 \\
      --normal-r1 N1_R1_001.fastq.gz \\
      --normal-r2 N1_R2_001.fastq.gz \\
      --cosmic /data/cosmic/COSMIC_v98.vcf.gz \\
      --dbsnp /data/dbsnp/dbsnp_155.hg38.vcf.gz \\
      --germline-resource /data/gnomad/gnomad.vcf.gz \\
      --threads 16

  # Tumour-only mode (no matched normal):
  python comprehensive_variant_calling.py \\
      --input-dir raw_fastqs/ \\
      --output-dir results/ \\
      --reference hg38 \\
      --tumour-sample TUMOUR_01 \\
      --tumour-r1 T1_R1_001.fastq.gz \\
      --tumour-r2 T1_R2_001.fastq.gz \\
      --cosmic /data/cosmic/COSMIC_v98.vcf.gz \\
      --threads 16

  # Dry run (show commands without executing):
  python comprehensive_variant_calling.py ... --dry-run

  # Skip QC (if FASTQs are already cleaned):
  python comprehensive_variant_calling.py ... --skip-steps qc

  # Skip specific steps:
  python comprehensive_variant_calling.py ... --skip-steps qc bqsr annotate

  # Resume after a crash, reusing whatever finished (same arguments as the
  # interrupted run -- --resume refuses if they changed):
  python comprehensive_variant_calling.py ... --resume

Output Structure
----------------
  results/
  ├── aligned/                  # Sorted, indexed BAMs
  │   ├── TUMOUR_01.sorted.bam
  │   └── NORMAL_01.sorted.bam
  ├── dedup/                    # Duplicate-marked BAMs
  │   ├── TUMOUR_01.dedup.bam
  │   └── NORMAL_01.dedup.bam
  ├── bqsr/                     # BQSR-corrected BAMs
  │   ├── TUMOUR_01.bqsr.bam
  │   └── NORMAL_01.bqsr.bam
  ├── cleaned_fastq/            # QC'd FASTQs from fastp
  ├── mutect2/                  # Mutect2 raw calls
  │   ├── TUMOUR_01.mutect2.vcf.gz
  │   ├── TUMOUR_01.mutect2.filtered.vcf.gz
  │   └── TUMOUR_01.mutect2.cosmic.vcf.gz
  ├── annotated/                # SnpEff annotated VCF
  │   └── TUMOUR_01.annotated.vcf
  ├── pcgr/                     # PCGR clinical report (HTML + TSV)
  ├── coverage/                 # Per-sample target coverage report
  │                             # (--coverage-bed; HTML/JSON/TSV)
  ├── metrics/                  # Duplicate rates, BQSR tables
  ├── stats/                    # Alignment statistics
  ├── logs/                     # Per-sample logs
  ├── fastp_reports/            # fastp HTML/JSON reports
  ├── reference/                # Only if the shared reference cache
  │                            # was unwritable -- see --reference-dir
  └── run_manifest_*.json       # Full pipeline record

FFPE MATERIAL
-------------
  Formalin fixation does two things this pipeline has to answer for.

  It DEAMINATES CYTOSINE, producing C>T (and G>A on the other strand)
  changes that are damage, not biology. STEP 6's read-orientation model is
  the defence: deamination artefacts appear on one read orientation, and
  --ob-priors lets FilterMutectCalls use that. On the FFPE panel this was
  developed against it removed 3,222 calls, 1,824 of them C>T or G>A.
  Never turn that step off on FFPE.

  It also FRAGMENTS DNA, which shows up as a short insert size (108 bp
  here, shorter than the reads themselves) and, more awkwardly, as
  foldback artefacts -- a damaged fragment end folds back on itself and
  is sequenced as its own reverse complement. See the foldback filter.

  Residual damage survives both. Among PASS SNVs on that sample, 59% of
  C>T/G>A calls sat below 10% VAF against 37% for every other
  substitution: the deamination tail is thinner after filtering but still
  there. On FFPE, treat a low-VAF C>T as damage until something else
  argues otherwise, and prefer --min-allele-fraction over trusting the
  raw PASS set.

TARGETED PANELS NEED --intervals
-------------------------------
  Nothing restricts calling to the capture regions unless you say so.
  Given no --intervals, Mutect2 traverses the WHOLE GENOME, and every
  capture protocol leaves stray off-target reads scattered across it, so
  it calls on them: measured on a 4.5 Mb panel, the unrestricted run
  produced 635,414 raw candidates and 127,081 PASS calls, the great
  majority of them off target and many standing on one or two reads.

  Two flags address the two halves of that problem:
    --intervals    the panel BED. Restricts BQSR and Mutect2. Use the
                   manufacturer's file, in the reference's contig naming.
    --min-depth    tags thin calls 'low_depth' after filtering. GATK's
                   thresholds are likelihood-based and have no depth
                   floor, so this is not redundant with them.

  Both are off by default, because guessing a target region or a depth
  cutoff on someone else's assay would be worse than doing nothing.

DEBUGGING NOTES (things that surprise people)
---------------------------------------------
  * --auto-discover ASSUMES ONE PATIENT PER DIRECTORY. The first pair
    found becomes the tumour; the SECOND becomes that tumour's matched
    normal; anything after is ignored. Nothing warns about any of it.
    Two patients in one directory therefore produce a somatic call set
    that is the difference between two people -- germline variants where
    they differ reported as somatic, shared real mutations subtracted --
    and the output looks entirely ordinary. Give each patient its own
    directory, or name the samples explicitly and use --normal-sample
    only for a genuine matched normal from the same person. The webapp
    refuses an ambiguous directory and offers a batch mode; this script,
    used directly, does not.

  * --skip-steps SEMANTICS DIFFER PER STEP. Two kinds exist:
      - "reuse the previous artefact": qc, align, dedup, bqsr. Skipping
        these makes the pipeline fall back to the prior stage's output, so
        the variable that would have been produced still points at a real
        file. Skipping "align" additionally REQUIRES the BAM (and .bai) to
        already exist -- see --aligned-dir.
      - "just don't do it": mutect2, contamination, msi, filter, cosmic,
        annotate. Skipping "filter" makes filtered_vcf fall back to the
        raw Mutect2 VCF; skipping "contamination" leaves both tables None,
        and FilterMutectCalls then simply omits the two arguments.
    If you add a step, decide which kind it is and wire the fallback in the
    matching `else:` branch, or downstream code will reference a file that
    was never created.

  * --resume vs --skip-steps: BOTH end up in skip_steps, but they mean
    opposite things about which file feeds the next stage. A step named
    in --skip-steps falls back to the PREVIOUS stage's artefact (skipping
    bqsr hands the dedup BAM to Mutect2). A step skipped by --resume
    already produced its OWN artefact, so that is what moves forward
    (resuming bqsr hands the recalibrated BAM to Mutect2). That is what
    the `resumed_steps` set is for; when adding a step, mirror the
    existing `elif "<step>" in resumed_steps:` branches or the run will
    silently use the wrong input.

  * WHICH FAILURES ARE FATAL: qc, align, dedup, bqsr, mutect2 and filter
    all sys.exit(1) -- they produce artefacts the rest of the run needs.
    contamination, msi, cosmic, annotate and pcgr are enrichment only, so
    they
    warn and continue (they also return None, meaning "tool not installed
    / skipped", which the callers deliberately distinguish from False,
    meaning "failed").

  * PCGR (step 10) is off unless --pcgr-refdata-dir is given, and it does
    NOT consume the SnpEff output: PCGR runs its own VEP annotation and
    validates the INFO fields of its input, so it is handed the same call
    set SnpEff was (COSMIC-annotated if COSMIC ran, else the filtered
    VCF). It lives in the sibling module pcgr_report.py, imported lazily
    so this script still runs if that file is absent.

  * PCGR's TMB is only trustworthy when it has depth/AF INFO tags to
    filter on, and Mutect2 supplies those as FORMAT fields instead. Pass
    --pcgr-lift-tags to copy them into INFO (writing a separate VCF into
    the pcgr/ directory; the called VCF is untouched) and have the tag
    names forwarded automatically. Without it, PCGR still reports but its
    TMB is computed over unfiltered calls and reads high -- run_pcgr()
    warns loudly when that is the case. See pcgr_report.py.

  * --dry-run performs NO expensive I/O: no reference download, no
    reference decompression, no lane concatenation, no tool execution. It
    prints the commands instead. Every one of those was a real bug at some
    point; if you add a side effect, gate it on args.dry_run.

  * A reference named rather than given as a path (--reference hg38) is
    looked up in the SHARED cache first -- --reference-dir, default
    ~/data/references, which is where install_pipeline.py puts hg38 and its
    indices. Only a genome missing from there is downloaded, and it is
    downloaded into that cache, never into output_dir. Writing it per-run
    was a real cost: a 3 GB download and a 1-2 hour bwa-mem2 index repeated
    for every new output directory, producing bytes identical to the ones
    the last run had already built.

  * A VCF CANNOT SAY "NOT SEQUENCED". A region with no coverage yields no
    variant, exactly like a region that is wild type, and PCGR reports
    both as silence -- so a capture dropout reads as a negative result.
    --coverage-bed measures the panel against the depth a call must reach
    and names the stretches that fall short, per sample, in
    <output-dir>/coverage. It is optional and the run proceeds without it;
    what it costs is the ability to tell those two silences apart.

  * COSMIC CONTIG NAMES ARE ENSEMBL-STYLE (1, 2, MT), while hg38 is
    UCSC-style (chr1, chr2, chrM). bcftools annotate matches on contig
    NAME, so an un-renamed COSMIC file annotates nothing and does it
    silently -- exit 0, no diagnostic, overlap count zero.
    warn_on_contig_mismatch() checks for this before step 8 runs; rename
    the COSMIC VCF with `bcftools annotate --rename-chrs` once, up front.

  * COSMIC IDs LIVE IN THE VCF *ID* COLUMN, not an INFO field. That is why
    annotate_cosmic() uses `bcftools annotate -c ID,INFO` and
    count_cosmic_overlaps() filters on `ID~"COS[VMN]"` (COSV is the
    current scheme; COSM/COSN are legacy). Using `-c INFO` alone
    transfers no COSMIC identifiers at all and the overlap count silently
    reads zero.

  * SAMPLE NAMES MUST MATCH THE @RG SM TAG. Mutect2 is given -tumor/-normal
    by name; if the BAM's read group says something else, GATK errors out.
    The names come from the manifest / auto-discovery, and align_sample()
    writes them into @RG, so they agree as long as the same source is used
    for both.

  * dirs["tmp"] doubles as the GATK --tmp-dir and as the scratch space for
    lane-merged FASTQs. The merged FASTQs are deleted on the success path
    only -- after a crash they are left behind on purpose, for inspection.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse      # Command-line argument parsing
import difflib       # Longest common substring for the foldback check
import glob          # File pattern matching (find FASTQs)
import gzip          # Gunzip of the auto-downloaded reference
import hashlib       # SHA-256 hashing for file provenance
import json          # Writing the run manifest
import os            # File system operations (paths, dirs)
import re            # Regular expressions (pattern matching)
import shlex         # Safe command string formatting
import shutil        # Tool version detection (shutil.which)
import subprocess    # Running external tools (BWA, samtools, GATK)
import sys           # System exit, version info
import threading     # Background drain of BWA stderr (deadlock prevention)
from datetime import datetime, timezone  # Timestamps for manifest


# =============================================================================
# CONSTANTS
# =============================================================================
# A valid BAM/FASTQ file is at least this size. Anything smaller is likely
# a truncated or failed write. 1 KB is conservative but catches most issues.
MIN_PLAUSIBLE_OUTPUT_BYTES = 1024

# Reference genomes we can auto-download. Currently only hg38 is supported.
# To add hg19/GRCh37, add an entry with the Broad Institute's public URL.
# Source bucket: gs://gcp-public-data--broad-references (public, no login).
# NOTE: the older https://s3.amazonaws.com/broad-references/... path that
# used to be here is dead (404), and that bucket has no .fasta.gz at all.
KNOWN_REFERENCES = {
    "hg38": {
        "url": "https://storage.googleapis.com/gcp-public-data--broad-references/"
               "hg38/v0/Homo_sapiens_assembly38.fasta",
        # Where install_pipeline.py puts this genome, relative to the
        # reference cache, and the name it gives it. Matching the installer
        # is the whole point: a run that names 'hg38' then finds the
        # already-indexed copy instead of spending two hours rebuilding one.
        "subdir": "hg38",
        "filename": "Homo_sapiens_assembly38.fasta",
    },
}

# Where downloaded genomes live when the run does not name a FASTA of its
# own. This is deliberately NOT inside output_dir: the reference and its
# indices are identical for every run, they cost ~20 GB and 1-2 hours of
# bwa-mem2 indexing to produce, and putting them under the output directory
# meant every new output directory paid that price again. $PIPELINE_REFERENCE_DIR
# or --reference-dir moves the cache; the default matches install_pipeline.py.
DEFAULT_REFERENCE_DIR = os.environ.get(
    "PIPELINE_REFERENCE_DIR", os.path.expanduser("~/data/references"))


# =============================================================================
# SECTION 1: COMMAND EXECUTION
# =============================================================================
# These functions run external tools (fastp, BWA, GATK, etc.) and handle
# their output, errors, and logging. This is the foundation of the pipeline.

def run_command(cmd, tag, log_path=None, dry_run=False):
    """
    Execute a shell command with real-time output streaming.

    WHY THIS DESIGN?
      - We use subprocess.Popen() instead of subprocess.run() because
        bioinformatics tools can run for hours. Popen streams output
        line-by-line so you see progress in real-time.
      - We merge stderr into stdout (stderr=subprocess.STDOUT) because
        most tools write progress/status to stderr.
      - We return the exit code instead of raising exceptions, so one
        failed sample doesn't abort the entire pipeline.

    ARGS:
        cmd:      List of command arguments (e.g. ["bwa-mem2", "mem", ...]).
        tag:      Short label for log messages (e.g. "BWA TUMOUR_01").
        log_path: Path to append log output to. None = no logging.
        dry_run:  If True, print the command but don't execute it.

    RETURNS:
        int: Exit code (0 = success, non-zero = failure).
    """
    # Format the command as a human-readable string for logging.
    printable = shlex.join(cmd)
    print(f"\n[{tag}] $ {printable}", flush=True)

    if dry_run:
        print(f"[{tag}] (dry run - not executed)")
        return 0

    # Open log file for appending (create if needed).
    log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if log_fh:
            log_fh.write(f"\n### {tag}\n### {printable}\n")

        # Popen with line-buffered text mode. The key insight is bufsize=1:
        # this enables line buffering so we get output as soon as a newline
        # arrives, not when the buffer fills up.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Merge stderr into stdout.
            text=True,
            bufsize=1,
        )

        # Print each line with a "    | " prefix for visual distinction.
        # This makes it easy to scan terminal output and find tool messages.
        for line in proc.stdout:
            line = line.rstrip("\n")
            print(f"    | {line}", flush=True)
            if log_fh:
                log_fh.write(line + "\n")

        # Wait for the process to finish and get its exit code.
        proc.wait()

        if proc.returncode != 0:
            print(f"[{tag}] ERROR: exit code {proc.returncode}")
        return proc.returncode

    except FileNotFoundError:
        # This error means the tool binary wasn't found on PATH.
        # Common causes: tool not installed, wrong conda env, typo in name.
        print(f"[{tag}] ERROR: '{cmd[0]}' not found on PATH.")
        print(f"        Is it installed? Check: which {cmd[0]}")
        return 127
    except OSError as exc:
        print(f"[{tag}] ERROR: {exc}")
        return 1
    finally:
        # Always close the log file, even if an error occurred.
        if log_fh:
            log_fh.close()


def tool_version(executable):
    """
    Capture a tool's version string for the run manifest.

    We try three common flags (--version, -v, -V) because different tools
    use different conventions. Returns the first line of whichever works.
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
    Compute SHA-256 hash of the first `max_bytes` of a file.

    PURPOSE: Provenance tracking. If a FASTQ gets corrupted or swapped,
    we can detect it by comparing hashes. Hashing only the first 64 MB
    keeps it fast while still being unique enough for detection.
    """
    h = hashlib.sha256()
    bytes_read = 0
    with open(path, "rb") as fh:
        while True:
            chunk_size = 1024 * 1024  # Read 1 MB at a time.
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


def concatenate(files, destination):
    """
    Concatenate gzip members into one file. Concatenated gzip streams are a
    valid gzip file, so this is safe for .fastq.gz without recompression.
    Plain (uncompressed) FASTQ concatenates trivially too.

    Mirrors fastq_qc_clean.py so lane-split reads can be merged identically.
    """
    with open(destination, "wb") as out:
        for path in files:
            with open(path, "rb") as fh:
                shutil.copyfileobj(fh, out, length=8 * 1024 * 1024)
    return destination


# =============================================================================
# SECTION 2: REFERENCE GENOME HANDLING
# =============================================================================

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
    FASTA counts: an empty or truncated file is ignored so a half-finished
    download is re-fetched rather than fed to bwa-mem2.
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
        # the cost of hashing 3 GB on every run.
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

    The shared cache, unless it cannot be created or written -- a read-only
    or unwritable data directory falls back to output_dir/reference, which
    is where every run used to put it. The fallback is announced, because a
    silent one would quietly reintroduce the per-run download this exists
    to stop.
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

    The user can provide:
      - A file path: /data/references/hg38.fa
      - A known name: 'hg38' (triggers auto-download from Broad Institute)

    DRY-RUN CONTRACT (important when debugging):
      Every expensive side effect in here is gated on args.dry_run. The
      hg38 download is ~3 GB and the decompression writes ~3 GB more, so a
      "preview" that performed them would be a nasty surprise. In dry-run
      mode we print what WOULD happen and hand back the path the reference
      would occupy; the caller (ensure_indices) is also dry-run aware, so a
      non-existent path is fine downstream.

    RETURNS:
        str: Absolute path to the reference FASTA.
    """
    ref = args.reference
    # ensure_reference() is called with the whole args namespace, so read
    # the dry-run flag from there rather than adding a parameter (keeps the
    # call sites in main() unchanged).
    dry_run = getattr(args, "dry_run", False)

    # Check if this is a known genome name (not a file path).
    if ref in KNOWN_REFERENCES and not os.path.isfile(ref):
        entry = KNOWN_REFERENCES[ref]
        url = entry["url"]

        # Look in the shared cache before deciding to download anything.
        # An installed genome is already indexed, so finding it here saves
        # the download AND the 1-2 hour bwa-mem2 index that would follow.
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
                        print(f"[ERROR] Failed to decompress reference: "
                              f"{exc}")
                        sys.exit(1)

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
        sys.exit(1)

    print(f"[INFO] Reference: {ref}")
    return ref


def ensure_indices(ref, threads, dry_run=False):
    """
    Create all required index files for the reference genome.

    THREE INDICES ARE NEEDED:
      1. BWA-MEM2 index (.bwt.2bit.64, .amb, .ann, .pac, .sa)
         Used by: bwa-mem2 mem
         Purpose: Fast k-mer lookup during alignment.

      2. samtools .fai index
         Used by: samtools faidx, many GATK tools
         Purpose: Fast chromosome coordinate lookup.

      3. GATK .dict sequence dictionary
         Used by: GATK BaseRecalibrator, Mutect2, etc.
         Purpose: Reference metadata for interval parsing.

    We check for each and only build what's missing. This saves hours
    on subsequent runs since indexing takes 1-2 hours for hg38.
    """
    fai = ref + ".fai"
    dict_file = os.path.splitext(ref)[0] + ".dict"

    # --- BWA-MEM2 index ---
    bwt_mem2 = ref + ".bwt.2bit.64"
    if not os.path.exists(bwt_mem2):
        print("[INFO] Building BWA-MEM2 index (1-2 hours for hg38) ...")
        code = run_command(["bwa-mem2", "index", ref],
                           tag="BWA-MEM2 index", dry_run=dry_run)
        if code != 0:
            print("[ERROR] BWA-MEM2 indexing failed.")
            sys.exit(1)

    # --- samtools .fai index ---
    if not os.path.exists(fai):
        print("[INFO] Building samtools .fai index ...")
        code = run_command(["samtools", "faidx", ref],
                           tag="samtools faidx", dry_run=dry_run)
        if code != 0:
            print("[ERROR] samtools faidx failed.")
            sys.exit(1)

    # --- GATK sequence dictionary ---
    if not os.path.exists(dict_file):
        if shutil.which("gatk") is None:
            print("[WARN] GATK not found. Skipping .dict creation.")
        else:
            print("[INFO] Building GATK sequence dictionary ...")
            code = run_command(
                ["gatk", "CreateSequenceDictionary",
                 "-R", ref, "--VERBOSITY", "ERROR"],
                tag="GATK dict", dry_run=dry_run,
            )
            if code != 0:
                print("[ERROR] GATK dict failed.")
                sys.exit(1)

    print("[INFO] All reference indices ready.")


# =============================================================================
# SECTION 3: FASTQ QC AND CLEANING (fastp)
# =============================================================================

def run_fastp(sample_name, r1, r2, out_r1, out_r2, reports_dir,
              args, log_path, dry_run=False):
    """
    Run fastp for adapter removal, quality filtering, and poly-G trimming.

    WHY FASTP?
      - It's the fastest QC tool available (10x faster than Trimmomatic).
      - It handles poly-G trimming natively, which is critical for
        NextSeq 1000 (2-colour chemistry generates poly-G tails).
      - It produces HTML and JSON reports for easy QC review.

    WHAT IT DOES:
      1. Removes adapter sequences (auto-detected or user-specified).
      2. Trims poly-G runs from read ends (NextSeq-specific artefact).
      3. Filters reads by mean quality and minimum length.
      4. Reports before/after statistics.

    ARGS:
        sample_name: Used for output file naming and report titles.
        r1, r2:      Input FASTQ paths.
        out_r1, out_r2: Output cleaned FASTQ paths.
        reports_dir: Directory for HTML/JSON reports.
        args:        Parsed command-line arguments.
        log_path:    Path to append log output.
        dry_run:     If True, show command without executing.

    RETURNS:
        int: fastp exit code (0 = success).
    """
    cmd = [
        "fastp",
        "-i", r1,              # Input Read 1
        "-I", r2,              # Input Read 2
        "-o", out_r1,          # Output cleaned Read 1
        "-O", out_r2,          # Output cleaned Read 2
        "--thread", str(args.threads),
        "--length_required", str(args.min_read_length),
        "--qualified_quality_phred", str(args.min_base_quality),
        "--unqualified_percent_limit", str(args.max_unqualified_pct),
        "--n_base_limit", str(args.max_n_bases),
        "--html", os.path.join(reports_dir, f"{sample_name}.fastp.html"),
        "--json", os.path.join(reports_dir, f"{sample_name}.fastp.json"),
        "--report_title", f"fastp: {sample_name}",
    ]

    # Mean-quality filter: drop reads whose average Phred score is too low.
    if args.average_qual > 0:
        cmd += ["--average_qual", str(args.average_qual)]

    # Adapter handling. For PE data, fastp can auto-detect adapters by
    # examining read overlap. --detect_adapter_for_pe enables this.
    if args.adapter_r1:
        cmd += ["--adapter_sequence", args.adapter_r1]
    if args.adapter_r2:
        cmd += ["--adapter_sequence_r2", args.adapter_r2]
    if not args.no_detect_adapter_pe:
        cmd += ["--detect_adapter_for_pe"]

    # Poly-G trimming: NextSeq 1000 uses 2-colour chemistry where G bases
    # are signaled by "no signal" -- so when the insert is shorter than
    # the read, you get a long poly-G tail. This must be trimmed.
    if not args.no_poly_g:
        cmd += ["--trim_poly_g", "--poly_g_min_len", str(args.poly_g_min_len)]
    else:
        cmd += ["--disable_trim_poly_g"]

    # 3' sliding-window quality trimming. OFF by default for low-VAF
    # somatic calling: aggressive end trimming distorts allele fractions.
    if args.cut_right:
        cmd += ["--cut_right",
                "--cut_right_mean_quality",
                str(getattr(args, "cut_right_mean_quality", 20))]

    # Overlap-based base correction. OFF by default: it interferes with
    # UMI consensus calling.
    if args.correction:
        cmd += ["--correction"]

    # Overrepresented-sequence analysis (ON by default; detect library
    # artefacts like primers or adapter dimers). fastp has no dedicated
    # "disable" flag -- the analysis is simply omitted when disabled, the
    # same behaviour as fastq_qc_clean.py.
    if not args.no_overrepresentation:
        cmd += ["--overrepresentation_analysis"]

    # UMI handling for duplex/consensus panels.
    if args.umi_loc:
        cmd += ["--umi", "--umi_loc", args.umi_loc]
        if args.umi_len:
            cmd += ["--umi_len", str(args.umi_len)]
        if args.umi_skip is not None:
            cmd += ["--umi_skip", str(args.umi_skip)]

    # Write discarded reads to a per-sample file for audit.
    if args.failed_out:
        cmd += ["--failed_out",
                os.path.join(reports_dir, f"{sample_name}.failed.fastq.gz")]

    # Add optional --reads_to_process if present (useful for debugging).
    if hasattr(args, "reads_to_process") and args.reads_to_process:
        cmd += ["--reads_to_process", str(args.reads_to_process)]

    return run_command(cmd, tag=f"fastp [{sample_name}]",
                       log_path=log_path, dry_run=dry_run)


def summarise_fastp_json(json_path):
    """
    Extract key metrics from the fastp JSON report.

    Returns a dict with read counts, quality scores, and filtering results.
    These go into the run manifest for reproducibility.
    """
    try:
        with open(json_path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        return {"error": str(exc)}

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


# =============================================================================
# SECTION 4: ALIGNMENT (BWA-MEM2)
# =============================================================================

def align_sample(sample_name, r1, r2, ref, output_bam, threads,
                 two_pass, log_path, dry_run=False):
    """
    Align paired-end reads to the reference genome using BWA-MEM2.

    ALIGNMENT PIPELINE:
      bwa-mem2 mem ... | samtools sort -o output.bam - | samtools index

    WHY PIPE DIRECTLY?
      - Avoids writing an unsorted BAM to disk (~3x disk savings).
      - Sort happens in memory, which is faster than disk I/O.

    TWO-PASS MODE (--two-pass / -2 flag):
      - Pass 1: Standard alignment, identifies difficult regions.
      - Pass 2: Re-aligns reads in difficult regions with tuned params.
      - Recommended for cancer: tumours have rearrangements that create
        novel junctions not in the reference.

    READ GROUP HEADER (@RG):
      GATK REQUIRES these fields for variant calling:
        ID: Sample identifier (used for tracking).
        SM: Sample name (appears in VCF output).
        PL: Platform (ILLUMINA for NextSeq).
        LB: Library prep ID.
        PU: Platform unit (flowcell + lane + barcode).

    ARGS:
        sample_name: Used for @RG SM tag and log messages.
        r1, r2:      Input FASTQ paths.
        ref:         Reference genome FASTA path.
        output_bam:  Output sorted BAM path.
        threads:     CPU threads to use.
        two_pass:    Enable BWA-MEM2 two-pass alignment.
        log_path:    Log file path.
        dry_run:     If True, show commands without executing.

    RETURNS:
        str: Output BAM path on success, None on failure.
    """
    # Build the Read Group header. The \\t produces literal \t characters
    # that BWA-MEM2 interprets as tab delimiters in the SAM header.
    rg_header = (
        f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
        f"\\tLB:{sample_name}_lib\\tPU:unknown"
    )

    # Build the BWA-MEM2 command.
    cmd = ["bwa-mem2", "mem", "-t", str(threads), "-R", rg_header]
    if two_pass:
        cmd.append("-2")  # Enable two-pass mode.
    cmd += [ref, r1, r2]

    print(f"\n[Align {sample_name}] Pipeline: bwa-mem2 | samtools sort | index",
          flush=True)

    if dry_run:
        print(f"[Align {sample_name}] (dry run)")
        return output_bam

    log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if log_fh:
            log_fh.write(f"\n### Align {sample_name}\n")

        # --- PIPELINE START ---
        # Step 1: Start BWA-MEM2. Its stdout (SAM format) will be piped
        # to samtools sort's stdin.
        bwa_proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Drain BWA's stderr from a background thread. BWA streams
        # progressive status lines to stderr; reading them only after the
        # process exits lets a long alignment fill the OS pipe buffer
        # (~64 KB), block BWA, and deadlock the pipeline.
        bwa_stderr_lines = []
        stderr_lock = threading.Lock()

        def _drain_bwa_stderr():
            for line in bwa_proc.stderr:
                with stderr_lock:
                    bwa_stderr_lines.append(line.rstrip("\n"))

        stderr_thread = threading.Thread(target=_drain_bwa_stderr, daemon=True)
        stderr_thread.start()

        # Step 2: samtools sort reads from stdin and writes sorted BAM.
        # We use (threads - 1) for sort because BWA uses 1 thread already.
        sort_proc = subprocess.Popen(
            ["samtools", "sort",
             "-@", str(max(1, threads - 1)),
             "-o", output_bam,
             "-"],  # "-" = read from stdin
            stdin=bwa_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        # Close bwa's stdout in parent so SIGPIPE propagates if sort exits.
        bwa_proc.stdout.close()

        # Stream sort output (progress messages).
        for line in sort_proc.stdout:
            line = line.rstrip("\n")
            print(f"    | {line}", flush=True)
            if log_fh:
                log_fh.write(line + "\n")

        sort_proc.wait()
        bwa_proc.wait()

        # BWA writes useful stats to stderr: mapping rates, read counts.
        stderr_thread.join()
        if bwa_stderr_lines:
            print("    | [bwa-mem2 stats]")
            for sline in bwa_stderr_lines:
                print(f"    |   {sline}")
            if log_fh:
                log_fh.write("[bwa-mem2 stderr]\n"
                             + "\n".join(bwa_stderr_lines) + "\n")

        # --- CHECK EXIT CODES ---
        # A non-zero exit from BWA usually means: wrong reference, corrupt
        # FASTQ, or out-of-memory (BWA needs ~32GB RAM for hg38).
        if bwa_proc.returncode != 0:
            print(f"[Align {sample_name}] ERROR: bwa-mem2 failed "
                  f"(exit {bwa_proc.returncode}).")
            return None
        if sort_proc.returncode != 0:
            print(f"[Align {sample_name}] ERROR: samtools sort failed.")
            return None

        # Step 3: Index the sorted BAM for random access.
        code = run_command(["samtools", "index", output_bam],
                           tag=f"Index {sample_name}", log_path=log_path)
        if code != 0:
            print(f"[Align {sample_name}] ERROR: samtools index failed.")
            return None

        return output_bam

    except FileNotFoundError as exc:
        print(f"[Align {sample_name}] ERROR: {exc}")
        return None
    finally:
        if log_fh:
            log_fh.close()


def collect_stats(bam_path, sample_name, stats_dir, log_path, dry_run=False):
    """
    Generate alignment statistics for the BAM file.

    Mirrors align_reads.py's collect_stats() so both stages write the same
    flagstat/idxstats/stats files under <output>/stats/ and record the same
    manifest fields (mapping_rate, duplicates, insert_size_peak), instead of
    the earlier stub that ran the commands and threw the output away.

    RETURNS:
        dict: Alignment statistics.
    """
    if dry_run:
        return {"dry_run": True}

    stats = {}
    os.makedirs(stats_dir, exist_ok=True)

    # --- samtools flagstat ---
    flagstat_path = os.path.join(stats_dir, f"{sample_name}.flagstat.txt")
    try:
        res = subprocess.run(
            ["samtools", "flagstat", bam_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(flagstat_path, "w") as fh:
                fh.write(res.stdout)
            for line in res.stdout.splitlines():
                if "in total" in line:
                    stats["total_reads"] = line.split("+")[0].strip()
                elif "mapped (" in line and "primary" not in line:
                    stats["mapped_reads"] = line.split("+")[0].strip()
                    pct = line.split("(")[1].split("%")[0]
                    stats["mapping_rate"] = pct + "%"
                elif "duplicates" in line and "primary" not in line:
                    stats["duplicates"] = line.split("+")[0].strip()
    except Exception as exc:
        print(f"[WARN] flagstat failed: {exc}")

    # --- samtools idxstats ---
    idxstats_path = os.path.join(stats_dir, f"{sample_name}.idxstats.txt")
    try:
        res = subprocess.run(
            ["samtools", "idxstats", bam_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(idxstats_path, "w") as fh:
                fh.write(res.stdout)
    except Exception as exc:
        print(f"[WARN] idxstats failed: {exc}")

    # --- samtools stats ---
    stats_path = os.path.join(stats_dir, f"{sample_name}.stats.txt")
    try:
        res = subprocess.run(
            ["samtools", "stats", bam_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=3600,
        )
        if res.returncode == 0:
            with open(stats_path, "w") as fh:
                fh.write(res.stdout)
            for line in res.stdout.splitlines():
                if line.startswith("IS") and "insert size" not in line:
                    parts = line.split("\t")
                    if len(parts) >= 7:
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

    stats.pop("_max_is_count", None)
    return stats


# =============================================================================
# SECTION 5: DUPLICATE MARKING
# =============================================================================

def mark_duplicates(input_bam, output_bam, metrics_path, tmp_dir,
                    use_spark, threads, log_path, dry_run=False):
    """
    Mark PCR and optical duplicates in the BAM file.

    WHAT ARE DUPLICATES?
      During library preparation, DNA molecules are amplified by PCR.
      If the same molecule is amplified multiple times, it creates
      "duplicates" -- reads with identical start/end positions and
      sequences. These inflate coverage and can cause false variant calls.

    OPTICAL DUPLICATES:
      On patterned flowcells (like NextSeq 1000), clusters that are
      very close together can be confused by the image analysis.
      These are "optical" duplicates and are marked separately.

    WHY MARK, NOT REMOVE?
      GATK recommends marking (not removing) because:
        - Duplicate reads still carry useful base quality information.
        - GATK Mutect2 handles duplicates internally.
        - Removing reads reduces effective coverage and hurts sensitivity.

    SPARK vs PICARD:
      - GATK MarkDuplicatesSpark: Uses Apache Spark for parallelism.
        Much faster on multi-core machines (~2-5x speedup).
      - Picard MarkDuplicates: Traditional single-threaded implementation.
        Used as fallback if Spark is not available.

    ARGS:
        input_bam:   Sorted BAM from alignment step.
        output_bam:  Output path for duplicate-marked BAM.
        metrics_path: Path for duplicate metrics text file.
        tmp_dir:     Temporary directory for intermediate files.
        use_spark:   Use Spark-accelerated version (faster).
        threads:     CPU threads for Spark.
        log_path:    Log file path.
        dry_run:     If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    os.makedirs(tmp_dir, exist_ok=True)

    if use_spark:
        cmd = [
            "gatk", "MarkDuplicatesSpark",
            "-I", input_bam,
            "-O", output_bam,
            "-M", metrics_path,
            "--tmp-dir", tmp_dir,
            # 'spark.local.cores' is not a Spark property. It was passed
            # via --conf, where Spark silently ignores keys it does not
            # know, so it never did anything and --threads never reached
            # this step. The knob that works is the Spark master URL:
            # local[N] runs N worker threads in-process.
            "--spark-master", f"local[{threads}]" if threads and threads > 0
                              else "local[*]",
            "--verbosity", "ERROR",
        ]
    else:
        # Picard fallback: slower but no Spark dependency.
        cmd = [
            "gatk", "MarkDuplicates",
            "-I", input_bam,
            "-O", output_bam,
            "-M", metrics_path,
            "--CREATE_INDEX", "TRUE",
            "--VALIDATION_STRINGENCY", "LENIENT",
            "--TMP_DIR", tmp_dir,
        ]

    code = run_command(cmd, tag="MarkDuplicates",
                       log_path=log_path, dry_run=dry_run)
    if code != 0:
        print("[ERROR] MarkDuplicates failed.")
        return False

    # Ensure BAM index exists (GATK MarkDuplicates doesn't always create it).
    index_path = output_bam + ".bai"
    if not os.path.exists(index_path):
        run_command(["samtools", "index", output_bam],
                    tag="Index dedup", log_path=log_path,
                    dry_run=dry_run)
    return True


# =============================================================================
# SECTION 5b: TARGET INTERVALS
# =============================================================================

def _interval_args(intervals, padding=0):
    """
    Build the GATK -L / --interval-padding arguments, or nothing at all.

    Centralised so every tool that takes intervals spells them the same
    way, and so "no intervals" stays a single empty-list case rather than
    an `if` at each call site.

    GATK accepts .bed, .interval_list, .list and bare "chr:start-end"
    strings here; the BED convention (0-based, half-open) is handled by
    GATK itself, which reads the extension rather than guessing.

    Padding matters on capture panels: probes cover the target but reads
    extend past its edges, so a variant at the first or last base is
    supported by reads that a zero-padded interval would exclude.

    ARGS:
        intervals: Path to an interval file, or None.
        padding:   Bases to extend each interval by, on both sides.

    RETURNS:
        list: Arguments to append to a GATK command line (possibly empty).
    """
    if not intervals:
        return []
    args = ["-L", intervals]
    if padding:
        args += ["--interval-padding", str(padding)]
    return args


# =============================================================================
# SECTION 6: BASE QUALITY SCORE RECALIBRATION (BQSR)
# =============================================================================

def run_bqsr(bam, ref, known_sites, output_bam, report_path,
             tmp_dir, log_path, dry_run=False,
             intervals=None, interval_padding=0):
    """
    GATK Base Quality Score Recalibration (BQSR).

    WHAT IS BQSR?
      Illumina sequencers report base quality scores (Phred scores)
      that estimate the probability of each base being incorrect.
      However, these scores have systematic biases:
        - Quality degrades over the read length.
        - Certain sequence contexts have predictable errors.
        - Machine calibration drifts over time.

      BQSR uses known variant sites (dbSNP, Mills indels) to model
      these biases and recalibrate the quality scores. This improves
      variant calling accuracy, especially for low-coverage data.

    WHY USE KNOWN SITES?
      - dbSNP: Known germline variants. If a base differs from the
        reference at a known variant site, it's probably real (not a
        sequencing error), so the quality score should be trusted.
      - Mills indels: Known indel sites from the 1000 Genomes Project.

    TWO-STEP PROCESS:
      1. BaseRecalibrator: Builds a recalibration model (table file).
      2. ApplyBQSR: Applies the model to produce a new BAM with
         corrected quality scores.

    ARGS:
        bam:          Input duplicate-marked BAM.
        ref:          Reference genome FASTA.
        known_sites:  List of VCF files with known variant sites.
        output_bam:   Output path for recalibrated BAM.
        report_path:  Path for the recalibration report table.
        tmp_dir:      Temporary directory.
        log_path:     Log file path.
        dry_run:      If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    recal_table = os.path.join(tmp_dir, os.path.basename(output_bam) + ".recal")

    # --- Step 1: BaseRecalibrator ---
    base_cmd = [
        "gatk", "BaseRecalibrator",
        "-R", ref,
        "-I", bam,
        "-O", recal_table,
        "--tmp-dir", tmp_dir,
        "--verbosity", "ERROR",
    ]
    for vcf in known_sites:
        base_cmd += ["--known-sites", vcf]

    # Restricting BQSR to the capture targets is not just a speed-up on
    # panel data: off-target reads are sparse and error-prone, and folding
    # them into the error model biases the recalibration that every
    # downstream call depends on. ApplyBQSR is deliberately NOT restricted
    # -- it must rewrite the whole BAM, or reads outside the targets would
    # silently vanish from the output.
    base_cmd += _interval_args(intervals, interval_padding)

    code = run_command(base_cmd, tag="BaseRecalibrator",
                       log_path=log_path, dry_run=dry_run)
    if code != 0:
        print("[ERROR] BaseRecalibrator failed.")
        return False

    # --- Step 2: ApplyBQSR ---
    apply_cmd = [
        "gatk", "ApplyBQSR",
        "-R", ref,
        "-I", bam,
        "-O", output_bam,
        "--bqsr-recal-file", recal_table,
        "--tmp-dir", tmp_dir,
        "--verbosity", "ERROR",
    ]
    code = run_command(apply_cmd, tag="ApplyBQSR",
                       log_path=log_path, dry_run=dry_run)
    if code != 0:
        print("[ERROR] ApplyBQSR failed.")
        return False

    # Index the recalibrated BAM.
    run_command(["samtools", "index", output_bam],
                tag="Index BQSR", log_path=log_path, dry_run=dry_run)
    return True


# =============================================================================
# SECTION 7: SOMATIC VARIANT CALLING (GATK Mutect2)
# =============================================================================

def run_mutect2(tumour_bam, tumour_name, ref, output_vcf,
                normal_bam=None, normal_name=None,
                germline_resource=None, panel_of_normals=None,
                f1r2_dir=None, tmp_dir=None,
                log_path=None, dry_run=False,
                intervals=None, interval_padding=0, extra_args=None):
    """
    GATK Mutect2: The gold standard for somatic SNV and indel calling.

    WHAT IS MUTECT2?
      Mutect2 is a Bayesian somatic variant caller that:
        1. Builds a model of the read data at each position.
        2. Compares tumour (and optionally normal) against the reference.
        3. Calculates the probability that any variant is somatic vs germline.

    TUMOUR-NORMAL vs TUMOUR-ONLY MODE:
      - Tumour-normal (preferred): Uses a matched normal BAM to subtract
        germline variants automatically. Much cleaner calls.
      - Tumour-only: Relies on --germline-resource (gnomAD) to identify
        and filter germline variants. More false positives.

    F1R2 TARBALLS:
      Mutect2 produces F1R2 files that record the first-in-pair (F1) and
      second-in-pair (R2) strand orientation of each variant call. This
      data is used by LearnReadOrientationModel to correct strand bias
      artefacts, which are especially common in FFPE samples.

    COSMIC INTEGRATION:
      COSMIC is applied AFTER calling by annotate_cosmic() (SnpSift) and in
      the dedicated STEP 10. GATK Mutect2 has no COSMIC-aware flag, so the
      COSMIC VCF is deliberately NOT passed here -- doing so would silently
      no-op. See annotate_cosmic().

    ARGS:
        tumour_bam:        Tumour BAM (BQSR-corrected).
        tumour_name:       Tumour sample name (must match @RG SM tag).
        ref:               Reference genome FASTA.
        output_vcf:        Output VCF path (.vcf.gz recommended).
        normal_bam:        Normal BAM (optional, for paired mode).
        normal_name:       Normal sample name (must match @RG SM tag).
        germline_resource: gnomAD VCF for germline filtering.
        panel_of_normals:  Panel of Normals VCF for artefact subtraction.
        extra_args:        Already-split extra GATK arguments, appended
                           verbatim after everything this function builds.
        f1r2_dir:          Directory for F1R2 tarballs.
        tmp_dir:           Temporary directory.
        log_path:          Log file path.
        dry_run:           If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)

    cmd = [
        "gatk", "Mutect2",
        "-R", ref,
        "-I", tumour_bam,
        "-tumor", tumour_name,
        "-O", output_vcf,
        "--tmp-dir", tmp_dir or "/tmp",
        "--verbosity", "ERROR",
    ]

    # F1R2 output for strand-bias modelling.
    if f1r2_dir:
        os.makedirs(f1r2_dir, exist_ok=True)
        f1r2_path = os.path.join(f1r2_dir, f"{tumour_name}.f1r2.tar.gz")
        cmd += ["--f1r2-tar-gz", f1r2_path]

    # Paired normal (strongly recommended).
    if normal_bam and normal_name:
        cmd += ["-I", normal_bam, "-normal", normal_name]

    # Germline resource (gnomAD) for tumour-only mode.
    if germline_resource:
        cmd += ["--germline-resource", germline_resource]

    # Panel of Normals for recurrent artefact subtraction.
    if panel_of_normals:
        cmd += ["-pon", panel_of_normals]

    # Target restriction. Without it Mutect2 walks the entire genome and
    # calls on the off-target reads that any capture protocol leaves
    # scattered everywhere -- on the panel this was validated against,
    # that difference is roughly 635k raw candidates versus a few
    # thousand, and the surviving PASS calls include singletons sitting on
    # one or two stray reads.
    cmd += _interval_args(intervals, interval_padding)

    # Appended LAST on purpose. GATK takes the last occurrence of a
    # repeated argument, so a setting passed through here overrides the
    # equivalent one built above rather than being silently ignored --
    # which is what an escape hatch has to do to be worth having.
    if extra_args:
        cmd += list(extra_args)

    code = run_command(cmd, tag="Mutect2", log_path=log_path, dry_run=dry_run)
    return code == 0


# =============================================================================
# SECTION 8: STRAND BIAS MODELLING
# =============================================================================

def run_learn_read_orientation(f1r2_tar, output_model, log_path=None,
                               dry_run=False):
    """
    LearnReadOrientationModel: Build a strand-bias correction model.

    WHY STRAND BIAS MATTERS:
      In FFPE (formalin-fixed, paraffin-embedded) samples, DNA damage
      during fixation causes C>T and G>A artefacts. These artefacts
      appear predominantly on one strand, creating strand bias.

      This tool uses the F1R2 data from Mutect2 to learn the strand
      orientation patterns and build a model that can distinguish
      real variants from FFPE artefacts.

    ARGS:
        f1r2_tar:     F1R2 tarball from Mutect2.
        output_model: Output path for the read orientation model.
        log_path:     Log file path.
        dry_run:      If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    cmd = [
        "gatk", "LearnReadOrientationModel",
        "-I", f1r2_tar,
        "-O", output_model,
        "--verbosity", "ERROR",
    ]
    code = run_command(cmd, tag="LearnReadOrientation",
                       log_path=log_path, dry_run=dry_run)
    return code == 0


# =============================================================================
# SECTION 8b: CONTAMINATION ESTIMATION
# =============================================================================

# Below this many usable common-SNP sites, a contamination estimate is
# reported but flagged as weak. This is a pragmatic floor, not a GATK
# threshold -- GATK publishes none. Whole-genome and exome runs clear it
# comfortably; targeted panels frequently do not, because they only cover
# the common SNPs that happen to fall inside their capture regions.
MIN_CONTAMINATION_SITES = 500

# Minimum read depth for a common-SNP site to count toward that total.
MIN_CONTAMINATION_SITE_DEPTH = 10


def _count_informative_pileup_sites(pileup_table, dry_run=False):
    """
    Count common-SNP sites in a GetPileupSummaries table with real depth.

    The table has one row per site: contig, position, ref_count, alt_count,
    other_alt_count, allele_frequency. Depth is the sum of the three counts.
    Rows beginning with '#' are metadata and 'contig' is the header.

    Returns 0 when the file cannot be read; the count is advisory and must
    never be the thing that stops a run.
    """
    if dry_run:
        return 0
    count = 0
    try:
        with open(pileup_table, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("contig"):
                    continue
                fields = line.split("\t")
                if len(fields) < 5:
                    continue
                try:
                    depth = int(fields[2]) + int(fields[3]) + int(fields[4])
                except ValueError:
                    continue
                if depth >= MIN_CONTAMINATION_SITE_DEPTH:
                    count += 1
    except OSError as exc:
        print(f"[WARN] Could not read pileup table: {exc}")
        return 0
    return count


def estimate_contamination(bam, common_variants, pileup_table,
                           contamination_table, segmentation_table,
                           tmp_dir, log_path=None, dry_run=False):
    """
    Estimate cross-sample contamination so FilterMutectCalls can use it.

    WHY THIS MATTERS
      Without a contamination estimate FilterMutectCalls does not guess:
      it runs with --contamination-estimate 0.0. The "contamination"
      FILTER is then not disabled -- measured on a real panel run it
      still tagged 10 calls out of 635,414, all of which had already
      failed other filters -- but it is testing every variant against an
      ASSUMED zero, so it cannot detect contamination that is actually
      there. Supplying the tables replaces the assumption with a
      measurement.

      This matters because a few percent of another person's DNA in the
      library produces a tail of low-VAF calls that look exactly like
      real subclonal somatic variants. In TUMOUR-ONLY mode there is no
      matched normal to contradict them, so this is one of the few
      handles available on that failure mode -- and the one that catches
      a sample swap or index hopping on a patterned flowcell.

      This is NOT the same thing as the read-orientation model in STEP 6.
      That one models FFPE/oxidation artefacts (the "orientation" filter);
      this models foreign DNA. Both are needed; neither substitutes.

    HOW IT WORKS
      GetPileupSummaries counts ref/alt reads at a set of COMMON biallelic
      germline SNPs whose population allele frequencies are known.
      CalculateContamination then looks at sites where the sample should be
      homozygous: any consistent minor allele there is foreign DNA, and its
      fraction is the contamination estimate.

      The same VCF is passed as both -V (the sites and their AFs) and -L
      (the intervals to walk), which is what GATK's own somatic workflow
      does -- walking the whole genome to visit ~8k sites would be waste.

    ADVISORY, NOT FATAL
      Returns (None, None) on any failure, and the caller then runs
      FilterMutectCalls without the tables -- i.e. exactly the behaviour
      this pipeline had before the step existed. A contamination estimate
      is worth having but is never worth losing the variant calls over.

    ARGS:
        bam:                  Recalibrated tumour BAM.
        common_variants:      Common biallelic germline SNP VCF, bgzipped
                              and tabix-indexed, with AF in INFO
                              (small_exac_common_3.hg38.vcf.gz).
        pileup_table:         Output path for the pileup summary.
        contamination_table:  Output path for the contamination estimate.
        segmentation_table:   Output path for the tumour segmentation.
        tmp_dir:              Scratch directory for GATK.
        log_path:             Log file path.
        dry_run:              If True, show commands without executing.

    RETURNS:
        tuple: (contamination_table, segmentation_table) on success, or
               (None, None) if the step could not be completed.
    """
    os.makedirs(tmp_dir, exist_ok=True)

    code = run_command(
        ["gatk", "GetPileupSummaries",
         "-I", bam,
         "-V", common_variants,
         "-L", common_variants,
         "-O", pileup_table,
         "--tmp-dir", tmp_dir,
         "--verbosity", "ERROR"],
        tag="GetPileupSummaries", log_path=log_path, dry_run=dry_run,
    )
    if code != 0:
        print("[WARN] GetPileupSummaries failed; continuing without a "
              "contamination estimate.")
        return None, None

    # HOW MANY SITES ACTUALLY INFORMED THE ESTIMATE?
    # A contamination figure carries no information about its own support,
    # and 0.0% from a handful of sites prints identically to 0.0% from
    # thousands. Targeted panels are the common case for a thin result:
    # they cover only the fraction of the common-SNP set that happens to
    # fall inside the capture regions. Count the usable sites so the
    # caller can judge the number instead of trusting it.
    informative = _count_informative_pileup_sites(pileup_table, dry_run)

    code = run_command(
        ["gatk", "CalculateContamination",
         "-I", pileup_table,
         "-O", contamination_table,
         "--tumor-segmentation", segmentation_table,
         "--tmp-dir", tmp_dir,
         "--verbosity", "ERROR"],
        tag="CalculateContamination", log_path=log_path, dry_run=dry_run,
    )
    if code != 0:
        print("[WARN] CalculateContamination failed; continuing without a "
              "contamination estimate.")
        return None, None

    if dry_run:
        return contamination_table, segmentation_table

    # Surface the number rather than burying it in a file nobody opens.
    # A malformed or empty table is not worth failing over, so any parse
    # problem is reported and the tables are still handed back -- GATK
    # wrote them, and FilterMutectCalls is the authority on reading them.
    try:
        with open(contamination_table, encoding="utf-8") as fh:
            rows = [ln.split("\t") for ln in fh.read().splitlines() if ln]
        for row in rows[1:]:
            if len(row) >= 3:
                fraction, error = float(row[1]), float(row[2])
                print(f"[INFO] Contamination: {fraction:.4%} "
                      f"(+/- {error:.4%}) from {informative} usable sites")
                if fraction >= 0.05:
                    print("[WARN] Contamination >= 5%. Low-VAF calls in "
                          "this sample are unreliable; check for a sample "
                          "swap or index hopping before reporting them.")
                elif informative < MIN_CONTAMINATION_SITES:
                    print(f"[WARN] Only {informative} common-SNP sites had "
                          f"usable depth (< {MIN_CONTAMINATION_SITES}).")
                    print("       Treat this estimate as WEAK EVIDENCE, not "
                          "as a clean result: a low figure here means the "
                          "test had little to work with, not necessarily "
                          "that the sample is uncontaminated.")
                    print("       Expected on targeted panels, which cover "
                          "only the common SNPs that fall inside their "
                          "capture regions.")
    except (OSError, ValueError, IndexError) as exc:
        print(f"[WARN] Could not read contamination table: {exc}")

    return contamination_table, segmentation_table


# =============================================================================
# SECTION 9: MUTECT2 CALL FILTERING
# =============================================================================

def filter_mutect2(mutect_vcf, ref, output_vcf,
                   orientation_model=None, contamination_table=None,
                   segmentation_table=None, min_allele_fraction=0.0,
                   log_path=None, dry_run=False):
    """
    FilterMutectCalls: Apply hard filters to Mutect2 output.

    WHY FILTER?
      Mutect2 emits every candidate it can support and defers judgement;
      this is the step that decides. FilterMutectCalls learns a threshold
      from the data and writes a reason into the FILTER column. PASS means
      the variant survived every filter. The tags actually produced are
      GATK's own -- among the ones that matter here:
        - weak_evidence:    log-odds below the learned threshold.
        - germline:         likely germline rather than somatic. Does the
                            heavy lifting in tumour-only mode.
        - orientation:      FFPE/oxidation artefact. Only reachable when
                            --ob-priors is supplied (STEP 6).
        - contamination:    explained by foreign DNA. Meaningful only when
                            --contamination-table is supplied (STEP 7).
                            Without it GATK does not disable the filter,
                            it just assumes a fraction of zero -- so a
                            few calls may still carry the tag while real
                            contamination goes undetected.
        - panel_of_normals: seen in the PoN, so recurrent artefact.
        - strand_bias, base_qual, map_qual, position, slippage,
          clustered_events, multiallelic, haplotype: read-level and
          locus-level artefact signatures.

      Note the pattern in the two marked above: those filters are only
      as good as the evidence file behind them. A VCF with few or no
      `orientation` or `contamination` calls does not mean the sample is
      clean -- it may mean the test was never really run.

    ARGS:
        mutect_vcf:          Raw Mutect2 VCF.
        ref:                 Reference genome FASTA.
        output_vcf:          Output filtered VCF path.
        orientation_model:   Read orientation model from
                             LearnReadOrientationModel (STEP 6).
        contamination_table: Contamination estimate from
                             CalculateContamination (STEP 7).
        min_allele_fraction: Minimum tumour allele fraction; calls below it
                             are tagged low_allele_frac. GATK's default is
                             0.0, i.e. no VAF floor at all. 0 disables.
        segmentation_table:  Tumour segmentation from the same step. GATK
                             uses it to apply the estimate per segment
                             rather than as one genome-wide number, which
                             matters where copy number varies.
        log_path:            Log file path.
        dry_run:             If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    cmd = [
        "gatk", "FilterMutectCalls",
        "-V", mutect_vcf,
        "-O", output_vcf,
        "-R", ref,
        "--verbosity", "ERROR",
    ]
    if orientation_model:
        cmd += ["--ob-priors", orientation_model]
    if min_allele_fraction and min_allele_fraction > 0:
        cmd += ["--min-allele-fraction", str(min_allele_fraction)]
    if contamination_table:
        cmd += ["--contamination-table", contamination_table]
    # Segmentation is only meaningful alongside the estimate it segments.
    if contamination_table and segmentation_table:
        cmd += ["--tumor-segmentation", segmentation_table]

    code = run_command(cmd, tag="FilterMutectCalls",
                       log_path=log_path, dry_run=dry_run)
    return code == 0


# =============================================================================
# SECTION 8c: MICROSATELLITE INSTABILITY
# =============================================================================

# Upstream calls MSI-H at >= 20% for tumour-only data. The paired caller uses
# 3.5%, and quoting that threshold against a tumour-only score is a common way
# to turn a stable sample into a false MSI-H call.
MSI_HIGH_THRESHOLD = 20.0

# Below this many covered sites the score is reported but flagged. Panels
# carry only the microsatellites that fall inside their capture regions, and
# a handful of unstable loci out of very few swings the percentage wildly.
MSI_MIN_SITES = 200


def run_msisensor2(bam, models_dir, out_prefix, intervals=None, threads=1,
                   log_path=None, dry_run=False):
    """
    Score microsatellite instability from the tumour BAM alone.

    WHY A SEPARATE TOOL
      PCGR will not do this. It restricts MSI prediction to WGS/WES
      tumour-control runs and, on a targeted or tumour-only query, logs a
      warning and omits the analysis -- so asking PCGR for MSI on a panel
      returns a report with no MSI section rather than an MSI answer.
      MSIsensor2 is built for the tumour-only case and reads the BAM.

    WHAT THE NUMBER MEANS
      The score is the percentage of covered microsatellite loci whose
      length distribution is unstable. Upstream calls MSI-H at >= 20% for
      TUMOUR-ONLY data. The older paired MSIsensor used 3.5%, and applying
      that figure here would call almost anything unstable.

    THE SITE COUNT MATTERS AS MUCH AS THE SCORE
      A panel only carries the microsatellites inside its capture regions.
      With a hundred-odd covered loci, a couple of noisy ones move the
      score by whole percentage points, so the count is returned alongside
      it and a thin one is flagged rather than quietly averaged away.

    FFPE
      Formalin fixation fragments DNA and deaminates cytosine, and damaged
      reads across a homopolymer can widen its apparent length
      distribution. Treat a borderline score on FFPE material with more
      suspicion than the same score on fresh-frozen.

    ARGS:
        bam:        Analysis-ready tumour BAM.
        models_dir: MSIsensor2 model directory for this genome build.
        out_prefix: Output path prefix; the tool adds _dis and _somatic.
        intervals:  Optional BED restricting the scan to the capture target.
        threads:    Parallel workers.
        log_path:   Log file path.
        dry_run:    If True, show the command without running it.

    RETURNS:
        dict with score, sites and unstable counts, or None when the step
        could not run. Never raises: MSI is enrichment, not a dependency.
    """
    if not models_dir:
        return None
    if shutil.which("msisensor2") is None:
        print("[WARN] msisensor2 not found on PATH; skipping MSI.")
        return None
    if not os.path.isdir(models_dir):
        print(f"[WARN] MSI models directory not found: {models_dir}")
        return None

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    cmd = ["msisensor2", "msi",
           "-M", models_dir,
           "-t", bam,
           "-o", out_prefix,
           "-b", str(max(1, threads))]
    if intervals:
        cmd += ["-e", intervals]

    code = run_command(cmd, tag="MSIsensor2", log_path=log_path,
                       dry_run=dry_run)
    if dry_run:
        return None
    if code != 0:
        print("[WARN] msisensor2 failed; continuing without an MSI score.")
        return None

    # Output is a two-line TSV: header, then total/unstable/percentage.
    try:
        with open(out_prefix, encoding="utf-8") as fh:
            rows = [ln.split() for ln in fh.read().splitlines() if ln.strip()]
    except OSError as exc:
        print(f"[WARN] Could not read the MSI result: {exc}")
        return None
    if len(rows) < 2 or len(rows[1]) < 3:
        print("[WARN] MSI result was not in the expected format.")
        return None

    try:
        sites, unstable, score = (int(rows[1][0]), int(rows[1][1]),
                                  float(rows[1][2]))
    except ValueError:
        print("[WARN] Could not parse the MSI result.")
        return None

    status = "MSI-H" if score >= MSI_HIGH_THRESHOLD else "MSS / MSI-low"
    print(f"[INFO] MSI: {score:.2f}% -- {unstable} unstable of {sites} "
          f"covered sites -> {status} "
          f"(tumour-only threshold {MSI_HIGH_THRESHOLD:.0f}%)")
    if sites < MSI_MIN_SITES:
        print(f"[WARN] Only {sites} microsatellite sites had usable "
              f"coverage (< {MSI_MIN_SITES}). A panel carries only the loci "
              f"inside its targets, and with this few a couple of noisy "
              f"ones move the score by whole points. Treat it as indicative.")
    return {"score_percent": score, "sites_covered": sites,
            "sites_unstable": unstable, "status": status,
            "threshold_percent": MSI_HIGH_THRESHOLD,
            "low_site_count": sites < MSI_MIN_SITES}


# =============================================================================
# SECTION 9b: DEPTH FLOOR
# =============================================================================

def apply_depth_floor(vcf_path, tumour_name, min_depth,
                      log_path=None, dry_run=False):
    """
    Flag calls whose tumour depth is below `min_depth` as low_depth.

    WHY THIS IS NEEDED SEPARATELY
      FilterMutectCalls has no minimum-depth filter. Its thresholds are
      likelihood-based, and a variant supported by one or two reads can
      clear them: on the panel this was validated against, PASS calls
      appeared at DP=1 with AF=0.667, i.e. two reads total. Those are not
      subclonal variants, they are noise that survived because nothing
      asked how much evidence stood behind them.

    SOFT FILTER, NOT DELETION
      The record is tagged, never removed. A variant caller's output is
      evidence, and silently dropping rows makes it impossible to tell
      "not called" from "called and discarded". Downstream steps that want
      only confident calls select on FILTER="PASS", which the tag already
      excludes them from.

    WHOSE DEPTH
      FORMAT/DP is per sample, so the tumour column has to be located by
      NAME. In tumour-normal mode Mutect2 writes both samples and their
      order is not guaranteed, so indexing blindly could threshold on the
      normal's depth instead. If the name is not in the header this
      returns without filtering rather than guessing.

      Records with a missing DP are left untouched: bcftools evaluates a
      missing value as not-matching, so absent evidence is not treated as
      evidence of low depth.

    ARGS:
        vcf_path:    Filtered VCF, bgzipped and indexed. Rewritten in place.
        tumour_name: Sample name whose DP is tested.
        min_depth:   Minimum depth; 0 or less disables the step entirely.
        log_path:    Log file path.
        dry_run:     If True, show commands without executing.

    RETURNS:
        bool: True on success or when the step is a no-op, False on failure.
    """
    if not min_depth or min_depth <= 0:
        return True
    if shutil.which("bcftools") is None:
        print("[WARN] bcftools not found; skipping the depth floor.")
        return True
    if not dry_run and not os.path.exists(vcf_path):
        print(f"[WARN] {vcf_path} not found; skipping the depth floor.")
        return True

    # Locate the tumour column by name -- see WHOSE DEPTH above.
    index = 0
    if not dry_run:
        try:
            res = subprocess.run(
                ["bcftools", "query", "-l", vcf_path],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                text=True, timeout=120,
            )
            samples = res.stdout.split()
        except Exception as exc:
            print(f"[WARN] Could not read sample names ({exc}); "
                  f"skipping the depth floor.")
            return True
        if tumour_name not in samples:
            print(f"[WARN] Sample '{tumour_name}' is not in {vcf_path} "
                  f"(found: {', '.join(samples) or 'none'}).")
            print("       Skipping the depth floor rather than guessing "
                  "which column holds the tumour.")
            return True
        index = samples.index(tumour_name)

    tmp_out = vcf_path + ".depthfloor.tmp.vcf.gz"
    expression = f"FORMAT/DP[{index}] < {min_depth}"
    code = run_command(
        ["bcftools", "filter",
         "-s", "low_depth",       # tag name written into FILTER
         "-m", "+",               # add to existing filters, do not replace
         "-e", expression,        # exclude-expression = what gets tagged
         "-O", "z", "-o", tmp_out, vcf_path],
        tag=f"Depth floor (DP < {min_depth})",
        log_path=log_path, dry_run=dry_run,
    )
    if code != 0:
        print("[WARN] Depth floor failed; leaving the VCF unfiltered.")
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
        return False
    if dry_run:
        return True

    os.replace(tmp_out, vcf_path)
    run_command(["bcftools", "index", "-f", "-t", vcf_path],
                tag="Index depth-floored VCF",
                log_path=log_path, dry_run=dry_run)

    # Report the cost, so the number is visible rather than inferred.
    try:
        res = subprocess.run(
            ["bcftools", "view", "-H", "-i", 'FILTER~"low_depth"', vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=300,
        )
        tagged = len(res.stdout.splitlines())
        print(f"[INFO] Depth floor: {tagged} call(s) tagged low_depth "
              f"(tumour DP < {min_depth}).")
    except Exception:
        pass
    return True


# =============================================================================
# SECTION 9c: FOLDBACK (PALINDROMIC) INSERTION ARTEFACTS
# =============================================================================

# Default threshold: how much of an inserted sequence must reverse-complement
# match the flanking reference before it is called an artefact. Measured on a
# real panel, genuine and artefactual insertions separate cleanly and nothing
# sits near the boundary -- artefacts had a MEDIAN match of 100% while ordinary
# insertions matched 24% in the forward direction and essentially nothing in
# reverse. 0.7 sits in the empty space between the two populations.
FOLDBACK_MIN_MATCH = 0.7

# Insertions shorter than this are not tested: a few bases will match
# something by chance in any 300 bp window, so the test says nothing.
FOLDBACK_MIN_INSERTION = 10

# Reference pulled either side of the breakpoint. The foldback point is not
# always exactly at the called position, so the window has to be wider than
# the insertion itself.
FOLDBACK_FLANK = 150


def _reverse_complement(seq):
    """Reverse complement, leaving anything that is not ACGTN alone."""
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _longest_common_substring_len(a, b):
    """
    Length of the longest exact substring shared by a and b.

    difflib is used rather than a hand-rolled DP table because the common
    case here is a near-exact match, which it finds quickly, and because a
    quadratic table over a 300 bp window costs more than it needs to.
    """
    if not a or not b:
        return 0
    if a in b:
        return len(a)
    return difflib.SequenceMatcher(None, a, b, autojunk=False) \
        .find_longest_match(0, len(a), 0, len(b)).size


def tag_foldback_insertions(vcf_path, ref, min_match=FOLDBACK_MIN_MATCH,
                            log_path=None, dry_run=False):
    """
    Tag insertions that are reverse-complement copies of their own flanks.

    WHAT THIS CATCHES
      A fragment whose single-stranded end folds back on itself and is
      extended during library preparation carries a sequence followed by
      its own reverse complement. Aligned, that renders as an INSERTION of
      the reverse complement of the adjacent reference. The reads are real
      and the alignment is reasonable; the variant is not.

      This is not a rare curiosity. On the panel this was written for,
      95.6% of insertions longer than 10 bp were reverse-complement copies
      of the flanking reference (median match: 100%), against 0.2% in the
      forward direction -- and they accounted for 206 of 270 PASS
      frameshift and stop_gained calls, including frameshifts in BRCA1,
      BRCA2, PALB2, PTEN and ATM that a clinical reporter had tiered as
      oncogenic. Nothing else in the pipeline removes them: they are not
      strand-biased, not low-depth, not in tandem repeats, and they reach
      allele fractions above 0.3, so no VAF floor separates them either.

    WHY IT IS A TAG AND NOT A DELETION
      Genuine foldback inversions occur in cancer genomes, particularly
      where BRCA-mediated repair has failed, so the pattern is evidence of
      artefact rather than proof of one. Every tested insertion also gets
      an INFO/FBMATCH value recording how much of it matched, so a
      reviewer can see the evidence instead of trusting the threshold.

    ARGS:
        vcf_path:  Filtered VCF, bgzipped and indexed. Rewritten in place.
        ref:       Reference FASTA, indexed (.fai).
        min_match: Matched fraction at or above which the FILTER tag is
                   added. 0 or less disables the step entirely.
        log_path:  Log file path.
        dry_run:   If True, report intent and change nothing.

    RETURNS:
        bool: True on success or when the step is a no-op.
    """
    if not min_match or min_match <= 0:
        return True
    for tool in ("bcftools", "samtools", "bgzip", "tabix"):
        if shutil.which(tool) is None:
            print(f"[WARN] {tool} not found; skipping the foldback check.")
            return True
    if dry_run:
        print(f"[Foldback] (dry run) would test insertions >= "
              f"{FOLDBACK_MIN_INSERTION} bp in {vcf_path} against "
              f"+/-{FOLDBACK_FLANK} bp of reference")
        return True
    if not os.path.exists(vcf_path):
        print(f"[WARN] {vcf_path} not found; skipping the foldback check.")
        return True

    # --- collect insertions worth testing -------------------------------
    try:
        res = subprocess.run(
            ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
             vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=600,
        )
        if res.returncode != 0:
            print("[WARN] Could not read the VCF; skipping the foldback check.")
            return True
    except Exception as exc:
        print(f"[WARN] Could not read the VCF ({exc}); skipping foldback.")
        return True

    candidates = []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom, pos, vref, valt = parts[0], parts[1], parts[2], parts[3]
        # Only the first ALT is tested. A multiallelic record tagged on its
        # first allele is still the right call: the site is suspect.
        alt = valt.split(",")[0]
        if len(alt) - len(vref) < FOLDBACK_MIN_INSERTION:
            continue
        candidates.append((chrom, int(pos), vref, alt))

    if not candidates:
        print("[INFO] Foldback check: no insertions long enough to test.")
        return True

    # --- one batched faidx call, not one per variant --------------------
    tmp_dir = os.path.dirname(os.path.abspath(vcf_path))
    regions_path = os.path.join(tmp_dir, ".foldback_regions.txt")
    with open(regions_path, "w", encoding="utf-8") as fh:
        for chrom, pos, _, _ in candidates:
            start = max(1, pos - FOLDBACK_FLANK)
            fh.write(f"{chrom}:{start}-{pos + FOLDBACK_FLANK}\n")
    try:
        res = subprocess.run(
            ["samtools", "faidx", ref, "-r", regions_path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=1800,
        )
        if res.returncode != 0:
            print("[WARN] samtools faidx failed; skipping the foldback check.")
            return True
    finally:
        try:
            os.remove(regions_path)
        except OSError:
            pass

    windows, current = [], []
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            if current:
                windows.append("".join(current))
            current = []
        else:
            current.append(line.strip().upper())
    if current:
        windows.append("".join(current))
    if len(windows) != len(candidates):
        print(f"[WARN] Reference context count mismatch "
              f"({len(windows)} vs {len(candidates)}); skipping foldback.")
        return True

    # --- score, then annotate -------------------------------------------
    ann_path = os.path.join(tmp_dir, ".foldback_ann.tab")
    tagged = 0
    with open(ann_path, "w", encoding="utf-8") as fh:
        for (chrom, pos, vref, alt), window in zip(candidates, windows):
            inserted = alt[len(vref):].upper()
            if not window or not inserted:
                continue
            hit = _longest_common_substring_len(
                _reverse_complement(inserted), window)
            fraction = hit / len(inserted)
            if fraction >= min_match:
                tagged += 1
            fh.write(f"{chrom}\t{pos}\t{vref}\t{alt}\t{fraction:.3f}\n")

    hdr_path = os.path.join(tmp_dir, ".foldback_hdr.txt")
    with open(hdr_path, "w", encoding="utf-8") as fh:
        fh.write('##INFO=<ID=FBMATCH,Number=1,Type=Float,Description='
                 '"Fraction of the inserted sequence matching the reverse '
                 'complement of the flanking reference. Near 1.0 indicates '
                 'a library-prep foldback artefact.">\n')

    def _cleanup():
        for path in (ann_path, ann_path + ".gz", ann_path + ".gz.tbi",
                     hdr_path):
            try:
                if os.path.exists(path):
                    os.remove(path)
            except OSError:
                pass

    for cmd, tag in (
        (["bgzip", "-f", ann_path], "Foldback bgzip"),
        (["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", ann_path + ".gz"],
         "Foldback tabix"),
    ):
        if run_command(cmd, tag=tag, log_path=log_path) != 0:
            print("[WARN] Could not prepare the foldback annotations; "
                  "leaving the VCF unchanged.")
            _cleanup()
            return True

    annotated = vcf_path + ".foldback.tmp.vcf.gz"
    code = run_command(
        ["bcftools", "annotate", "-a", ann_path + ".gz", "-h", hdr_path,
         "-c", "CHROM,POS,REF,ALT,INFO/FBMATCH",
         "-O", "z", "-o", annotated, vcf_path],
        tag="Foldback annotate", log_path=log_path,
    )
    if code != 0:
        print("[WARN] Foldback annotation failed; leaving the VCF unchanged.")
        _cleanup()
        if os.path.exists(annotated):
            os.remove(annotated)
        return True

    filtered = vcf_path + ".foldback2.tmp.vcf.gz"
    code = run_command(
        ["bcftools", "filter", "-s", "foldback", "-m", "+",
         "-e", f"INFO/FBMATCH>={min_match}",
         "-O", "z", "-o", filtered, annotated],
        tag=f"Foldback filter (FBMATCH >= {min_match})", log_path=log_path,
    )
    os.remove(annotated)
    if code != 0:
        print("[WARN] Foldback filtering failed; leaving the VCF unchanged.")
        _cleanup()
        if os.path.exists(filtered):
            os.remove(filtered)
        return True

    os.replace(filtered, vcf_path)
    run_command(["bcftools", "index", "-f", "-t", vcf_path],
                tag="Index foldback-tagged VCF", log_path=log_path)
    _cleanup()

    print(f"[INFO] Foldback check: {tagged} of {len(candidates)} insertion(s) "
          f"tagged 'foldback' (reverse-complement match >= {min_match:.2f}).")
    if tagged:
        print("       These are library-prep artefacts, not variants. They "
              "are tagged, not removed; INFO/FBMATCH records the evidence "
              "for every insertion tested.")
    return True


# =============================================================================
# SECTION 10: COSMIC DATABASE ANNOTATION
# =============================================================================

def _vcf_header_contigs(path):
    """
    Return the set of contig names declared in a VCF header.

    Used only for the compatibility check below, so a failure to read the
    header is reported as "unknown" (empty set) rather than an error: the
    check is advisory and must never be the thing that stops a run.
    """
    if shutil.which("bcftools") is None:
        return set()
    try:
        res = subprocess.run(
            ["bcftools", "view", "-h", path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=120,
        )
        if res.returncode != 0:
            return set()
    except Exception:
        return set()

    contigs = set()
    for line in res.stdout.splitlines():
        if not line.startswith("##contig="):
            continue
        match = re.search(r"[<,]ID=([^,>]+)", line)
        if match:
            contigs.add(match.group(1))
    return contigs


def warn_on_contig_mismatch(vcf_in, cosmic_vcf):
    """
    Warn when the called VCF and the COSMIC VCF use different contig
    naming conventions.

    WHY THIS EXISTS
      bcftools annotate matches records by contig NAME. COSMIC distributes
      its GRCh38 VCFs with Ensembl-style names (1, 2, ... X, Y, MT) while
      this pipeline aligns to hg38, whose contigs are UCSC-style (chr1,
      chr2, ... chrX, chrY, chrM). Hand bcftools the two together and it
      transfers nothing at all -- and it does so SILENTLY: exit code 0, no
      diagnostic, an output VCF that looks fine and carries not one COSMIC
      identifier. The overlap count then reads zero, which is easy to
      misread as "this tumour has no known drivers".

      That is a bad failure to discover late, so spend two header reads on
      it up front.

    NOT FATAL. COSMIC annotation is enrichment, like SnpEff and PCGR, and a
    mismatch might even be deliberate. This warns and lets the run proceed,
    matching how the rest of the enrichment steps behave.

    RETURNS:
        bool: True if a mismatch was detected and warned about.
    """
    called = _vcf_header_contigs(vcf_in)
    cosmic = _vcf_header_contigs(cosmic_vcf)

    # Either header was unreadable or declared no contigs -- nothing can be
    # concluded, so say nothing.
    if not called or not cosmic:
        return False

    if called & cosmic:
        return False

    print("[WARN] CONTIG NAMING MISMATCH -- COSMIC annotation will do "
          "nothing.")
    print(f"       Called VCF uses : {', '.join(sorted(called)[:4])} ...")
    print(f"       COSMIC VCF uses : {', '.join(sorted(cosmic)[:4])} ...")
    print("       The two share no contig names, so bcftools annotate will")
    print("       match no records and transfer no COSMIC IDs (silently,")
    print("       with exit code 0). COSMIC ships Ensembl-style names while")
    print("       hg38 is UCSC-style. Rename the COSMIC file once:")
    print("         { for c in $(seq 1 22) X Y; do echo \"$c chr$c\"; done; \\")
    print("           echo \"MT chrM\"; } > chr_map.txt")
    print("         bcftools annotate --rename-chrs chr_map.txt -Oz \\")
    print("             -o cosmic.chr.vcf.gz cosmic.vcf.gz")
    print("         tabix -p vcf cosmic.chr.vcf.gz")
    print("       then pass the renamed file to --cosmic. Continuing anyway.")
    return True


def annotate_cosmic(vcf_in, vcf_out, cosmic_vcf, log_path=None,
                    dry_run=False):
    """
    Annotate variants with COSMIC database using bcftools.

    WHAT IS COSMIC?
      COSMIC (Catalogue Of Somatic Mutations In Cancer) is the world's
      largest database of somatic mutations found in human cancer,
      maintained by the Wellcome Sanger Institute.

      It contains:
        - Over 20 million coding mutations.
        - Data from >700,000 tumour samples.
        - 30+ cancer types.
        - Mutation frequency, sample count, and mutation type.

    WHY ANNOTATE WITH COSMIC?
      1. IDENTIFY KNOWN DRIVERS: Variants already in COSMIC are likely
         real somatic mutations (not artefacts).
      2. PRIORITISE VARIANTS: Variants in COSMIC with high sample counts
         are more likely to be clinically relevant.
      3. FILTER ARTIFACTS: Variants absent from COSMIC (but at COSMIC
         positions) may be artefacts to investigate further.

    WHAT THIS FUNCTION ADDS:
      - COSMIC mutation ID (e.g. COSM12345).
      - Number of samples with this mutation.
      - Total samples analysed at this position.
      - Mutation type (SNV, indel, etc.).

    REFERENCE BUILD WARNING:
      COSMIC VCFs are genome-build specific. If your reference is hg38
      but your COSMIC VCF is hg19, bcftools will fail with coordinate
      mismatches. Make sure they match!

    ARGS:
        vcf_in:    Input VCF (filtered Mutect2 calls).
        vcf_out:   Output VCF with COSMIC annotations.
        cosmic_vcf: Path to COSMIC VCF file.
        log_path:  Log file path.
        dry_run:   If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure.
    """
    if shutil.which("bcftools") is None:
        print("[WARN] bcftools not found; skipping COSMIC annotation.")
        return None

    if not os.path.exists(cosmic_vcf):
        print(f"[ERROR] COSMIC VCF not found: {cosmic_vcf}")
        return None

    # Catch the silent no-op before spending time on it: if the two files
    # disagree about contig naming, annotate matches nothing. Advisory only.
    if not dry_run:
        warn_on_contig_mismatch(vcf_in, cosmic_vcf)

    # bcftools annotate transfers fields from an annotation VCF.
    # -a: annotation VCF file.
    # -c ID,INFO: copy the ID column -- COSMIC's COSMxxxxx mutation IDs live
    #             there, NOT in an INFO field -- plus all INFO fields.
    # -O z: output compressed VCF.
    cmd = [
        "bcftools", "annotate",
        "-a", cosmic_vcf,
        "-c", "ID,INFO",
        "-o", vcf_out,
        "-O", "z",
        vcf_in,
    ]

    print(f"\n[COSMIC] $ {shlex.join(cmd)}", flush=True)
    if dry_run:
        print("[COSMIC] (dry run)")
        return True

    log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if log_fh:
            log_fh.write(f"\n### COSMIC annotation\n### {shlex.join(cmd)}\n")

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
            print(f"[COSMIC] ERROR: exit code {proc.returncode}")
            return False

        # Index the annotated VCF for downstream tools.
        run_command(["bcftools", "index", vcf_out],
                    tag="Index COSMIC VCF", log_path=log_path,
                    dry_run=dry_run)
        return True

    except Exception as exc:
        print(f"[COSMIC] ERROR: {exc}")
        return False
    finally:
        if log_fh:
            log_fh.close()


def count_cosmic_overlaps(vcf_path, cosmic_vcf, log_path=None,
                          dry_run=False):
    """
    Count how many called variants overlap COSMIC entries.

    This produces a simple summary:
      - Total variants called.
      - Variants found in COSMIC.
      - Variants NOT in COSMIC (potentially novel).
      - Overlap percentage.

    Useful for QC: a high overlap (e.g. >20% for whole-exome) suggests
    the calling pipeline is working well. A low overlap may indicate
    quality issues or that the sample is a cancer type with few COSMIC
    entries.

    ARGS:
        vcf_path:   Filtered VCF with COSMIC annotations.
        cosmic_vcf: COSMIC VCF for overlap counting.
        log_path:   Log file path.
        dry_run:    If True, show commands without executing.

    RETURNS:
        dict: Overlap statistics.
    """
    if not os.path.exists(vcf_path) or shutil.which("bcftools") is None:
        return {}

    if dry_run:
        return {"dry_run": True}

    stats = {}

    # Count total variants in the VCF.
    try:
        res = subprocess.run(
            ["bcftools", "view", "-H", "-c", "1", vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=300,
        )
        if res.returncode != 0:
            stats["total_variants"] = "error"
        else:
            stats["total_variants"] = len(res.stdout.splitlines())
    except Exception:
        stats["total_variants"] = "error"

    # Count variants with a COSMIC ID attached. annotate_cosmic() copies the
    # COSMIC VCF's ID column (via "-c ID,INFO"), so the mutation ID ends up
    # in the VCF ID field, not an INFO tag.
    #
    # THREE ID PREFIXES, NOT ONE. COSMIC changed its identifier scheme at
    # v90: current releases use genomic "COSV" IDs, while "COSM" (coding)
    # and "COSN" (non-coding) are the legacy forms, retained in older files
    # and in INFO/LEGACY_ID of newer ones. Matching only "COSM" -- as this
    # did -- reports 0 overlaps against any COSMIC v90+ file even when the
    # annotation itself worked perfectly. Keep all three in the class.
    # No shell involved, so quotes in the filter expression survive intact.
    try:
        res = subprocess.run(
            ["bcftools", "view", "-H", "-i", 'ID~"COS[VMN]"', vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=300,
        )
        if res.returncode != 0:
            stats["cosmic_overlapping"] = "error"
        else:
            cosmic_count = len(res.stdout.splitlines())
            stats["cosmic_overlapping"] = cosmic_count
            if isinstance(stats.get("total_variants"), int) and stats["total_variants"] > 0:
                pct = round(100.0 * cosmic_count / stats["total_variants"], 2)
                stats["cosmic_overlap_pct"] = pct
    except Exception:
        stats["cosmic_overlapping"] = "error"

    return stats


# =============================================================================
# SECTION 11: SnpEff ANNOTATION
# =============================================================================

def run_snpeff(vcf_in, vcf_out, genome="hg38", memory="8g", log_path=None,
               dry_run=False):
    """
    Annotate variants with SnpEff for gene and consequence information.

    WHAT IS SnpEff?
      SnpEff is a variant annotation tool that adds:
        - GENE NAME: Which gene does this variant affect?
        - CONSEQUENCE TYPE: Missense, nonsense, synonymous, intronic, etc.
        - IMPACT PREDICTION: HIGH / MODERATE / LOW / MODIFIER.

    IMPACT CATEGORIES:
      - HIGH:    Disruptive (nonsense, frameshift, splice donor/acceptor).
      - MODERATE: Potentially damaging (missense, in-frame indel).
      - LOW:     Probably benign (synonymous, intronic).
      - MODIFIER: Regulatory or non-coding (upstream, downstream, intergenic).

    REQUIREMENT:
      The SnpEff database for your genome must be downloaded first:
        snpEff download hg38

    HEAP SIZE IS NOT OPTIONAL
      bioconda ships snpEff behind a Python launcher whose default is
      "-Xms512m -Xmx1g". One gigabyte cannot hold the hg38 database, so
      the run dies partway through with

          java.lang.OutOfMemoryError: Java heap space

      after having already written tens of megabytes of valid-looking VCF.
      The launcher forwards any argument starting with "-Xm" to the JVM,
      so an explicit -Xmx is passed here. It must come BEFORE the genome
      name: everything the launcher does not recognise as a JVM option is
      handed to snpEff itself as a positional argument.

    A FAILED RUN LEAVES NO OUTPUT
      Because snpEff streams to stdout, a crash leaves a TRUNCATED VCF on
      disk -- and a truncated VCF is indistinguishable from a complete one
      to anything that only checks "does the file exist and is it
      non-empty", which is exactly what the --resume detector does. The
      partial file is therefore deleted on any failure, so a later resume
      redoes the step instead of silently annotating from a VCF that stops
      in the middle of chr2.

    ARGS:
        vcf_in:    Input VCF (COSMIC-annotated if available).
        vcf_out:   Output annotated VCF.
        genome:    Genome name for SnpEff database lookup.
        memory:    JVM max heap for snpEff (e.g. "8g"). See above.
        log_path:  Log file path.
        dry_run:   If True, show commands without executing.

    RETURNS:
        bool: True on success, False on failure, None if skipped.
    """
    if shutil.which("snpEff") is None:
        print("[WARN] SnpEff not found; skipping annotation.")
        return None

    def _discard_partial():
        """Remove a truncated output so --resume cannot mistake it for done."""
        try:
            if os.path.exists(vcf_out):
                os.remove(vcf_out)
                print(f"[SnpEff] Removed incomplete output {vcf_out}")
        except OSError as exc:
            print(f"[SnpEff] WARNING: could not remove {vcf_out}: {exc}")
            print("[SnpEff] Delete it by hand before using --resume: it is "
                  "truncated, and resume would treat it as finished.")

    # SnpEff writes to stdout by default. We redirect to a file.
    cmd = ["snpEff", f"-Xmx{memory}", genome, "-noStats", "-noLog", vcf_in]
    print(f"\n[SnpEff] $ {shlex.join(cmd)} > {vcf_out}", flush=True)
    if dry_run:
        print("[SnpEff] (dry run)")
        return True

    try:
        with open(vcf_out, "w", encoding="utf-8") as out_fh:
            proc = subprocess.Popen(
                cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True,
            )
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                print(f"[SnpEff] ERROR: exit code {proc.returncode}")
                print(f"[SnpEff] stderr: {stderr.strip()}")
                if "OutOfMemoryError" in (stderr or ""):
                    print(f"[SnpEff] The {memory} heap was not enough. "
                          f"Raise it with --snpeff-memory (e.g. "
                          f"--snpeff-memory 16g).")
                _discard_partial()
                return False
        return True
    except Exception as exc:
        print(f"[SnpEff] ERROR: {exc}")
        _discard_partial()
        return False


# =============================================================================
# SECTION 12: VCF STATISTICS
# =============================================================================

def vcf_stats(vcf_path, dry_run=False):
    """
    Run bcftools stats on a VCF and extract summary metrics.

    Returns a dict with key statistics like:
      - Number of records (variants).
      - Number of SNVs, indels, multiallelic sites.
      - Ti/Tv ratio (transition/transversion; expected ~2.1 for exomes).
    """
    if not os.path.exists(vcf_path) or shutil.which("bcftools") is None:
        return {"error": "VCF or bcftools not found"}
    if dry_run:
        return {"dry_run": True}

    try:
        res = subprocess.run(
            ["bcftools", "stats", vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=300,
        )
        stats = {}
        for line in res.stdout.splitlines():
            if line.startswith("SN"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    key = parts[2].strip().rstrip(":")
                    stats[key] = parts[-1].strip()
        return stats
    except Exception as exc:
        return {"error": str(exc)}


# =============================================================================
# SECTION 13: FASTQ DISCOVERY
# =============================================================================

LANE_RE = re.compile(r"_L\d{3}(?=_|$)")


def derive_tokens(r1_pattern, r2_pattern):
    """
    Robustly derive R1/R2 tokens from the two glob patterns.

    This mirrors the logic in fastq_qc_clean.py so that both scripts derive
    IDENTICAL sample names from the same input files -- a requirement for
    easy chaining (run QC first, feed the output to variant calling).

    Works for '*_R1_*.fastq*' / '*_R2_*.fastq*' (-> '_R1_' / '_R2_') and
    also for '*_1.fq.gz' / '*_2.fq.gz' (-> '_1.' / '_2.') style conventions:
    we split both patterns on '*' and take the SINGLE differing literal
    segment as the token.

    Returns:
        tuple: (r1_token, r2_token) or (None, None) if patterns are
               incompatible.
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


def find_paired_fastqs(input_dir, r1_pattern, r2_pattern, recursive=False):
    """
    Discover R1/R2 FASTQ pairs in the input directory.

    By default, looks for:
      - Read 1: *_R1_*.fastq* (e.g. Sample_R1_001.fastq.gz)
      - Read 2: *_R2_*.fastq* (e.g. Sample_R2_001.fastq.gz)

    IMPORTANT: This now uses the same robust token-derivation logic as
    fastq_qc_clean.py, and derives sample names identically. For lane-split
    files (_L001, _L002, ...), lane identifiers are stripped so that all
    lanes of a sample group under one name.

    RETURNS:
        tuple: (pairs, warnings)
          - pairs: list of dicts with keys: sample, r1, r2, is_lane_split
          - warnings: list of warning strings.
    """
    warnings = []

    # Build the glob pattern.
    if recursive:
        pattern = os.path.join(input_dir, "**", r1_pattern)
    else:
        pattern = os.path.join(input_dir, r1_pattern)
    r1_files = sorted(glob.glob(pattern, recursive=recursive))

    # Derive tokens using the SAME robust logic as fastq_qc_clean.py.
    r1_token, r2_token = derive_tokens(r1_pattern, r2_pattern)
    if r1_token is None or r2_token is None:
        warnings.append(
            f"Could not derive R1/R2 tokens from '{r1_pattern}' and "
            f"'{r2_pattern}'. They must differ in exactly one literal "
            f"segment, e.g. '*_R1_*.fastq*' and '*_R2_*.fastq*'.")
        return [], warnings

    pairs = []
    seen = {}

    for r1 in r1_files:
        directory, base = os.path.split(r1)

        # Split on the token and rebuild the R2 filename.
        if r1_token not in base:
            warnings.append(
                f"'{base}' matched the pattern but has no '{r1_token}' "
                f"token; skipping.")
            continue
        # Rebuild R2 by replacing only the read token, based on the last
        # occurrence, exactly as fastq_qc_clean.py does.
        r2_base = r2_token.join(base.rsplit(r1_token, 1))
        r2 = os.path.join(directory, r2_base)

        if not os.path.exists(r2):
            warnings.append(f"No R2 match for {r1} (expected {r2}); skipping.")
            continue

        # Lane-strip the sample name for grouping (S1_L001 and S1_L002 both
        # become "S1"), matching the fastq_qc_clean.py convention.
        stem = LANE_RE.sub("", base.rsplit(r1_token, 1)[0])

        # Detect if this is a lane-split occurrence BEFORE stripping.
        is_lane = bool(LANE_RE.search(base.rsplit(r1_token, 1)[0]))

        # Duplicate detection uses the full basename minus the FASTQ
        # extension as identity. This keeps lanes distinct whether the
        # lane token precedes (S2_L001_R1.fastq.gz) or follows
        # (SAME_R1_L001.fastq.gz) the read token, while still catching true
        # duplicates (identical filename twice).
        identity = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base)
        if identity in seen:
            warnings.append(
                f"Duplicate sample name '{identity}': {r1} collides with "
                f"{seen[identity]}. Skipping the duplicate to avoid silently "
                f"overwriting output.")
            continue
        seen[identity] = r1

        pairs.append({
            "sample": stem,
            "r1": r1,
            "r2": r2,
            "is_lane_split": is_lane,
        })

    # Group lane-split files under a single merged sample entry. Each lane's
    # r1/r2 is retained under the "lanes" list so the caller can concatenate
    # them (see main(), auto-discover mode). A single-lane sample keeps its
    # normal r1/r2 keys.
    merged = {}
    for pair in pairs:
        merged.setdefault(pair["sample"], []).append(pair)
    final_pairs = []
    for sample, group in merged.items():
        if len(group) > 1:
            final_pairs.append({
                "sample": sample,
                "r1": group[0]["r1"],
                "r2": group[0]["r2"],
                "is_lane_split": True,
                "lanes": [{"r1": p["r1"], "r2": p["r2"]} for p in group],
            })
        else:
            final_pairs.append(group[0])

    return final_pairs, warnings


# =============================================================================
# SECTION 13b: RESUME SUPPORT
# =============================================================================
# Re-running this pipeline after a crash otherwise repeats duplicate
# marking, BQSR and calling from scratch -- hours of work whose outputs are
# already sitting on disk. --resume detects finished steps and folds them
# into the existing --skip-steps machinery, so it reuses the same
# well-tested fallback paths rather than introducing a parallel one.
#
# Resume is OPT-IN, unlike the skip-if-output-exists behaviour in
# fastq_qc_clean.py and align_reads.py. Those stages are per-sample and
# their output depends on a handful of QC settings; a stage-3 artefact
# depends on the reference, known-sites resources, filter thresholds and
# more. Silently reusing a BQSR BAM built with a different --dbsnp would
# produce a quietly wrong call set, so the user has to ask for it.

# Parameters that do NOT change what lands on disk. Everything else is
# compared against the previous run before any output is reused.
_RESUME_IGNORED_PARAMS = frozenset({
    "threads", "dry_run", "resume", "resume_force", "skip_steps",
    "input_dir", "output_dir", "manifest",
    # PCGR only reads the finished call set and writes its own directory,
    # so its settings cannot invalidate an earlier step's output.
    "pcgr_refdata_dir", "vep_dir", "pcgr_assay", "pcgr_tumour_site",
    "pcgr_target_size_mb", "pcgr_estimate_tmb", "pcgr_estimate_msi",
    "pcgr_estimate_signatures",
    "pcgr_lift_tags", "pcgr_tumor_dp_tag", "pcgr_tumor_af_tag",
    "pcgr_legacy_v1", "pcgr_extra_args",
    # The coverage check only reads the finished BAMs and writes its own
    # directory, so changing the BED or its thresholds invalidates nothing
    # earlier. --reference-dir likewise: it says where the reference is
    # kept, not which reference is used.
    "coverage_bed", "coverage_min_depth", "coverage_min_mapq",
    "coverage_min_baseq", "reference_dir",
    # Panel profiles are applied BEFORE this comparison runs, so what gets
    # compared is the settings the profile produced -- interval padding,
    # the VAF floor, the skipped steps -- rather than its name. Comparing
    # the name as well would refuse a resume for a profile that was renamed
    # while setting exactly the same values, and would MISS a profile that
    # kept its name while its contents changed. The values are the honest
    # thing to compare; these are the flags that only say where they came
    # from.
    "panel", "panel_file", "panel_bed", "list_panels", "describe_panel",
    "new_panel_template", "save_panel_as", "panel_id",
})


def _has_index(path):
    """
    True if an index for `path` sits beside it.

    An index is the cheapest reliable "this file was written completely"
    signal we have: a crashed GATK or bcftools leaves the data file behind
    but never gets as far as indexing it.
    """
    for suffix in (".tbi", ".csi", ".bai", ".idx"):
        if os.path.exists(path + suffix):
            return True
    # samtools also writes foo.bai alongside foo.bam (extension replaced).
    root, ext = os.path.splitext(path)
    return ext == ".bam" and os.path.exists(root + ".bai")


def _complete(path, need_index=False):
    """True if `path` exists, is non-empty, and (optionally) is indexed."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    return _has_index(path) if need_index else True


def detect_completed_steps(dirs, tumour, normal, args):
    """
    Work out which pipeline steps already have complete output on disk.

    A step counts as done only when EVERY artefact it owns is present --
    for a tumour/normal pair that means both samples. A half-finished step
    is reported as not done, so it is simply re-run.

    RETURNS:
        (done, detail): a set of step names, and a {step: explanation} map
        for printing.
    """
    done, detail = set(), {}
    samples = [tumour] + ([normal] if normal else [])

    def all_samples(builder, need_index=False):
        """True when the artefact exists for every sample in the run."""
        return all(_complete(builder(s["name"]), need_index) for s in samples)

    # --- qc: cleaned FASTQs for every sample -------------------------
    if all(_complete(os.path.join(dirs["cleaned"], f"{s['name']}_R{r}.clean.fastq.gz"))
           for s in samples for r in (1, 2)):
        done.add("qc")
        detail["qc"] = "cleaned FASTQs present"

    # --- align / dedup / bqsr: indexed BAMs --------------------------
    for step, folder, suffix in (("align", "align", ".sorted.bam"),
                                 ("dedup", "dedup", ".dedup.bam"),
                                 ("bqsr", "bqsr", ".bqsr.bam")):
        if all_samples(lambda n, f=folder, x=suffix:
                       os.path.join(dirs[f], f"{n}{x}"), need_index=True):
            done.add(step)
            detail[step] = f"indexed {suffix} present"

    # --- the tumour-only, single-artefact steps ----------------------
    t = tumour["name"]
    checks = [
        ("mutect2", os.path.join(dirs["mutect2"], f"{t}.mutect2.vcf.gz"), True),
        ("filter", os.path.join(dirs["mutect2"],
                                f"{t}.mutect2.filtered.vcf.gz"), True),
        ("cosmic", os.path.join(dirs["mutect2"],
                                f"{t}.mutect2.cosmic.vcf.gz"), True),
        ("annotate", os.path.join(dirs["annotated"], f"{t}.annotated.vcf"),
         False),
    ]
    for step, path, need_index in checks:
        if _complete(path, need_index):
            done.add(step)
            detail[step] = f"{os.path.basename(path)} present"

    # --- contamination: BOTH tables, neither indexed -----------------
    # CalculateContamination writes the estimate and the segmentation
    # separately, and FilterMutectCalls is given them as a pair, so one
    # without the other is a half-finished step and must be redone.
    contam = os.path.join(dirs["metrics"], f"{t}.contamination.table")
    segments = os.path.join(dirs["metrics"], f"{t}.segments.table")
    if _complete(contam) and _complete(segments):
        done.add("contamination")
        detail["contamination"] = "contamination + segments tables present"

    # --- pcgr: any report file in the pcgr directory -----------------
    if glob.glob(os.path.join(dirs["pcgr"], "*.html")):
        done.add("pcgr")
        detail["pcgr"] = "report present"

    # COSMIC, contamination and PCGR are only meaningful when they were
    # asked for; do not claim they are "done" when they were never going
    # to run.
    if not args.cosmic:
        done.discard("cosmic")
        detail.pop("cosmic", None)
    if not getattr(args, "contamination_resource", None):
        done.discard("contamination")
        detail.pop("contamination", None)
    if not args.pcgr_refdata_dir:
        done.discard("pcgr")
        detail.pop("pcgr", None)

    return done, detail


def previous_run_parameters(output_dir):
    """
    Load the parameters recorded by the most recent run in `output_dir`.

    Manifest filenames carry a zero-padded UTC timestamp, so lexical order
    is chronological order.

    RETURNS:
        (params, manifest_path) or (None, None) when there is nothing to
        compare against.
    """
    manifests = sorted(glob.glob(os.path.join(output_dir,
                                              "pipeline_manifest_*.json")))
    if not manifests:
        return None, None
    try:
        with open(manifests[-1], encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None, manifests[-1]
    return data.get("parameters"), manifests[-1]


def compare_run_parameters(previous, current):
    """
    Diff two parameter sets, ignoring settings that cannot change output.

    RETURNS:
        list of (name, previous_value, current_value), sorted by name.
    """
    diffs = []
    for key in sorted(set(previous) | set(current)):
        if key in _RESUME_IGNORED_PARAMS:
            continue
        before, after = previous.get(key, "<unset>"), current.get(key, "<unset>")
        # JSON round-trips lists as lists but tuples as lists too; compare
        # on the JSON-ish shape so --known-indels does not false-positive.
        if isinstance(before, (list, tuple)):
            before = list(before)
        if isinstance(after, (list, tuple)):
            after = list(after)
        if before != after:
            diffs.append((key, before, after))
    return diffs


def plan_resume(dirs, tumour, normal, args, output_dir):
    """
    Decide which steps --resume may skip, refusing on parameter drift.

    Reusing an artefact built with different settings is the one genuinely
    dangerous failure mode here: the run would finish successfully and
    report a call set assembled from mismatched stages. So when the
    previous run's parameters differ, this aborts and makes the user
    choose, rather than guessing.

    RETURNS:
        set: step names to add to skip_steps.
    """
    previous, manifest_path = previous_run_parameters(output_dir)
    if previous:
        diffs = compare_run_parameters(previous, vars(args))
        if diffs:
            print(f"\n[ERROR] --resume: parameters differ from the previous "
                  f"run recorded in {os.path.basename(manifest_path)}:")
            for name, before, after in diffs:
                print(f"          {name}: {before!r} -> {after!r}")
            if not args.resume_force:
                print("        Reusing outputs built with different settings "
                      "would mix stages from two different analyses.")
                print("        Re-run without --resume to redo the pipeline, "
                      "or pass --resume-force if you are certain the change "
                      "does not affect the finished steps.")
                sys.exit(1)
            print("[WARN] --resume-force given; reusing the existing outputs "
                  "anyway. The result will mix settings from both runs.")
    elif manifest_path is None:
        print("[INFO] --resume: no previous run manifest in the output "
              "directory, so parameters could not be compared.")

    done, detail = detect_completed_steps(dirs, tumour, normal, args)
    if not done:
        print("[INFO] --resume: no completed steps found; running everything.")
        return set()

    print(f"\n[INFO] --resume: reusing {len(done)} completed step(s):")
    for step in sorted(done):
        print(f"          {step:9s} {detail[step]}")
    print("       Delete the corresponding output, or drop --resume, to "
          "force a step to re-run.")
    return done


# =============================================================================
# SECTION 14: ARGUMENT PARSING
# =============================================================================

def build_parser():
    """Build the command-line argument parser with all options."""
    parser = argparse.ArgumentParser(
        description="Comprehensive somatic variant calling pipeline for "
                    "cancer DNA (Illumina NextSeq 1000).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Input / Output ---
    io = parser.add_argument_group("input / output")
    io.add_argument("-i", "--input-dir",
                    help="Directory containing raw FASTQ files (required in "
                         "explicit and --auto-discover modes; ignored when "
                         "--manifest is provided).")
    io.add_argument("-o", "--output-dir", required=True,
                    help="Directory for all pipeline outputs.")
    io.add_argument("-r", "--reference", required=True,
                    help="Reference genome FASTA or known name (e.g. 'hg38').")
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
                         "When provided, cleaned FASTQ paths, sample names, "
                         "and lane-merged status are read from it, allowing "
                         "the QC step to be skipped and reusing the previous "
                         "QC output. This enables seamless chaining: "
                         "run fastq_qc_clean.py, then pass its manifest here.")

    # --- Sample specification ---
    # NOTE: not required=True because --manifest / --auto-discover modes
    # resolve sample names from the manifest or filenames. The explicit
    # mode requires them (enforced in validate_args).
    sample = parser.add_argument_group("sample specification")
    sample.add_argument("--tumour-sample", default=None,
                        help="Sample name for the tumour. Required in the "
                             "explicit-input mode; optional with "
                             "--manifest/--auto-discover (auto-picked).")
    sample.add_argument("--tumour-r1", default=None,
                        help="Read 1 FASTQ for the tumour. Required in the "
                             "explicit-input mode only.")
    sample.add_argument("--tumour-r2", default=None,
                        help="Read 2 FASTQ for the tumour. Required in the "
                             "explicit-input mode only.")
    sample.add_argument("--normal-sample", default=None,
                        help="Sample name for the matched normal.")
    sample.add_argument("--normal-r1", default=None,
                        help="Read 1 FASTQ for the matched normal.")
    sample.add_argument("--normal-r2", default=None,
                        help="Read 2 FASTQ for the matched normal.")

    # --- FASTQ patterns ---
    fq = parser.add_argument_group("FASTQ patterns")
    fq.add_argument("--r1-pattern", default="*_R1_*.fastq*",
                    help="Glob pattern for Read 1 files.")
    fq.add_argument("--r2-pattern", default="*_R2_*.fastq*",
                    help="Glob pattern for Read 2 files.")
    fq.add_argument("--recursive", action="store_true",
                    help="Search subdirectories recursively for FASTQs.")
    fq.add_argument("--merge-lanes", action="store_true",
                    help="Accepted for interface compatibility with "
                         "fastq_qc_clean.py. In auto-discover mode "
                         "lane-split samples (_L001.._L004) are ALWAYS "
                         "concatenated (there is no option to disable "
                         "merging), so this flag has no effect.")
    fq.add_argument("--auto-discover", action="store_true",
                    help="Auto-discover FASTQ pairs from input-dir "
                         "(overrides --tumour-r1/r2 and --normal-r1/r2). "
                         "ASSUMES ONE PATIENT PER DIRECTORY: the first pair "
                         "found is the tumour, the second is treated as its "
                         "MATCHED NORMAL, and any others are ignored -- all "
                         "silently. Two patients in one directory yields "
                         "somatic calls on the difference between two "
                         "people. Use one directory per patient.")

    # --- Panel / assay profile ---
    # Everything below this line used to have to be set flag by flag, from
    # memory, per kit. A profile is the same settings named once: see
    # panel_profiles.py for what each one carries and why. The profile only
    # ever fills in what the command line did not state, so adding --panel
    # to a command that already works cannot change it.
    panel = parser.add_argument_group(
        "panel / assay profile",
        "Configure the run by naming the assay instead of setting nine "
        "flags. Explicit flags always win over the profile.")
    panel.add_argument("--panel", default=None, metavar="ID",
                       help="Panel profile to configure this run from, e.g. "
                            "'illumina-tso500', 'thermo-oncomine-cav3', "
                            "'qiagen-qiaseq-dna', or one of the generic "
                            "shapes ('generic-capture', 'generic-amplicon', "
                            "'generic-hotspot', 'ctdna-capture', 'wes', "
                            "'wgs'). Run --list-panels for the full set. The "
                            "profile chooses interval padding, the VAF and "
                            "depth floors, whether duplicate marking and "
                            "BQSR are meaningful for that chemistry, the UMI "
                            "layout and which report statistics the "
                            "footprint can support -- each of which is "
                            "silent when it is wrong.")
    panel.add_argument("--panel-file", action="append", default=[],
                       metavar="JSON",
                       help="Additional profile file or directory to load, "
                            "repeatable. This is how a kit that is not built "
                            "in gets added: write the JSON once and every "
                            "run can name it. Profiles are also picked up "
                            "automatically from ~/.config/cancer_pipeline/"
                            "panels and $CANCER_PIPELINE_PANEL_DIR.")
    panel.add_argument("--panel-bed", default=None, metavar="BED",
                       help="The kit's target BED, used for BOTH --intervals "
                            "and --coverage-bed. They are almost always the "
                            "same file, and setting only one of them is the "
                            "commonest way to end up with a run that either "
                            "calls genome-wide or makes no coverage "
                            "statement. Either flag, given explicitly, still "
                            "wins over this. No profile ships a BED -- "
                            "vendor BEDs are licensed, version-specific "
                            "content that must come from your kit.")
    panel.add_argument("--list-panels", action="store_true",
                       help="List every known panel profile and exit.")
    panel.add_argument("--describe-panel", default=None, metavar="ID",
                       help="Print everything a profile sets, with its "
                            "caveats, and exit.")
    panel.add_argument("--new-panel-template", default=None, metavar="JSON",
                       help="Write a starting-point profile for a kit that "
                            "is not built in, and exit. Combine with --panel "
                            "to base it on an existing profile.")
    panel.add_argument("--save-panel-as", default=None, metavar="JSON",
                       help="After the settings are resolved, save them as a "
                            "reusable profile. The intended workflow is to "
                            "tune a run with ordinary flags under --dry-run "
                            "until it is right, then save it so the next "
                            "operator names the assay instead of "
                            "reconstructing it. Target BED paths are "
                            "deliberately not stored.")
    panel.add_argument("--panel-id", default=None, metavar="ID",
                       help="Id to give the profile written by "
                            "--save-panel-as (default: the file's stem).")

    # --- Resources ---
    res = parser.add_argument_group("annotation resources")
    res.add_argument("--dbsnp", default=None,
                     help="dbSNP VCF for BQSR and Mutect2.")
    res.add_argument("--germline-resource", default=None,
                     help="gnomAD VCF for germline filtering (tumour-only).")
    res.add_argument("--panel-of-normals", default=None,
                     help="Panel of Normals VCF.")
    res.add_argument("--contamination-resource", default=None,
                     help="Common biallelic germline SNP VCF (with AF in "
                          "INFO) used to estimate cross-sample "
                          "contamination. Supplying it enables STEP 7, so "
                          "FilterMutectCalls tests against a measured "
                          "fraction instead of an assumed zero. Broad's "
                          "small_exac_common_3.hg38.vcf.gz is the standard "
                          "choice; it must use the same contig naming as "
                          "the reference.")
    res.add_argument("--known-indels", nargs="*", default=[],
                     help="Known indel VCFs for BQSR.")
    res.add_argument("--intervals", default=None,
                     help="Target intervals (.bed / .interval_list / .list) "
                          "restricting BQSR and Mutect2. On capture panels "
                          "this is what separates real calls from off-target "
                          "noise: without it Mutect2 walks the whole genome "
                          "and calls on the stray reads every capture "
                          "protocol leaves behind. Use the panel "
                          "manufacturer's BED, in the same contig naming as "
                          "the reference.")
    res.add_argument("--coverage-bed", default=None,
                     help="Target regions to check for coverage before the "
                          "clinical report is written (BED). Every region is "
                          "measured against the depth a call must reach, and "
                          "the stretches that fall short are named in a "
                          "per-sample HTML report under <output-dir>/coverage. "
                          "Without it the run makes no coverage statement at "
                          "all: an unsequenced exon and a wild-type exon look "
                          "identical in the VCF and identical in the report. "
                          "Optional -- the run proceeds without the step when "
                          "no BED can be provided.")
    res.add_argument("--coverage-min-depth", type=int, default=None,
                     help="Depth a base must reach to count as covered "
                          "(default: --min-depth, the depth below which a "
                          "call would be filtered out anyway).")
    res.add_argument("--coverage-min-mapq", type=int, default=20,
                     help="Ignore reads below this mapping quality when "
                          "measuring coverage (default: 20).")
    res.add_argument("--coverage-min-baseq", type=int, default=20,
                     help="Ignore bases below this base quality when "
                          "measuring coverage (default: 20).")
    res.add_argument("--interval-padding", type=int, default=0,
                     help="Bases to extend each interval by on both sides "
                          "(GATK --interval-padding). 50-100 is usual for "
                          "capture panels: reads run past the target edges, "
                          "so a variant at the first or last base needs the "
                          "flanking reads to be callable. Ignored without "
                          "--intervals.")
    res.add_argument("--msi-models", default=None, metavar="DIR",
                     help="MSIsensor2 model directory for this genome build "
                          "(e.g. ~/data/msisensor2/models_hg38). Supplying it "
                          "enables STEP 8, which scores microsatellite "
                          "instability from the tumour BAM alone. PCGR will "
                          "not do this on a panel or a tumour-only query -- "
                          "it restricts MSI to WGS/WES tumour-control runs "
                          "and omits the section without comment. Upstream "
                          "calls MSI-H at >=20%% for tumour-only data; the "
                          "3.5%% figure belongs to the older paired caller "
                          "and would call almost anything unstable.")
    res.add_argument("--foldback-min-match", type=float,
                     default=FOLDBACK_MIN_MATCH,
                     help="Tag an insertion 'foldback' when this fraction of "
                          "it reverse-complement matches the flanking "
                          "reference. Those are library-prep artefacts: a "
                          "fragment end folded back on itself and was "
                          "extended, so the read carries a sequence followed "
                          "by its own reverse complement. Nothing else here "
                          "removes them -- they are not strand-biased, not "
                          "low-depth, not in tandem repeats, and reach VAF "
                          "above 0.3. ON by default (%(default)s) because "
                          "the failure mode is silent: on the panel this was "
                          "written for they produced 206 of 270 PASS "
                          "frameshift/stop calls, including BRCA1, BRCA2 and "
                          "PALB2 frameshifts a clinical reporter tiered as "
                          "oncogenic. Tagged, never removed; INFO/FBMATCH "
                          "records the evidence. 0 disables.")
    res.add_argument("--min-allele-fraction", type=float, default=0.0,
                     help="Minimum tumour allele fraction; FilterMutectCalls "
                          "tags anything below it low_allele_frac. GATK's "
                          "own default is 0.0 -- no VAF floor at all. On a "
                          "panel this is the main handle on PCR-slippage "
                          "indels, which cluster at 3-10%% VAF and otherwise "
                          "PASS in numbers that swamp the real calls. Set it "
                          "from the assay's validated limit of detection, "
                          "NOT from what makes the output look tidy: too "
                          "high and genuine subclonal variants disappear. "
                          "0 disables.")
    res.add_argument("--min-depth", type=int, default=0,
                     help="Tag calls with tumour FORMAT/DP below this as "
                          "'low_depth' after FilterMutectCalls, which has no "
                          "depth filter of its own and will PASS a variant "
                          "supported by two reads. Records are tagged, never "
                          "removed. 0 disables it.")
    res.add_argument("--mutect2-extra-args", default=None,
                     help="Extra arguments forwarded verbatim to Mutect2, as "
                          "one quoted string. The escape hatch for a kit "
                          "whose chemistry needs a caller setting this "
                          "script has no flag of its own for -- e.g. "
                          "'--dont-use-soft-clipped-bases true' on an "
                          "amplicon panel whose primers were not trimmed. "
                          "Nothing here validates it: GATK does, at the "
                          "start of the calling step.")
    res.add_argument("--cosmic", default=None,
                     help="COSMIC VCF for known somatic cancer variants. "
                          "Download from https://cancer.sanger.ac.uk/cosmic")

    # --- QC options ---
    qc = parser.add_argument_group("QC / trimming options")
    qc.add_argument("--min-read-length", "--min-length", type=int, default=50,
                    help="Minimum read length after trimming (fastp). Same "
                         "flag name as fastq_qc_clean.py "
                         "(--min-length still accepted).")
    qc.add_argument("--min-base-quality", "--qualified-quality", type=int,
                    default=15,
                    help="Minimum base Phred quality for qualified bases. "
                         "Same flag name as fastq_qc_clean.py "
                         "(--qualified-quality still accepted).")
    qc.add_argument("--max-unqualified-pct", "--unqualified-percent",
                    type=int, default=40,
                    help="Max percentage of unqualified bases per read. Same "
                         "flag name as fastq_qc_clean.py "
                         "(--unqualified-percent still accepted).")
    qc.add_argument("--average-qual", type=int, default=20,
                    help="Minimum mean read quality (0 = disabled).")
    qc.add_argument("--max-n-bases", "--n-base-limit", type=int, default=5,
                    help="Maximum N bases per read. Same flag name as "
                         "fastq_qc_clean.py (--n-base-limit still accepted).")
    qc.add_argument("--adapter-r1", default=None,
                    help="Explicit Read 1 adapter sequence.")
    qc.add_argument("--adapter-r2", default=None,
                    help="Explicit Read 2 adapter sequence.")
    qc.add_argument("--no-detect-adapter-pe", action="store_true",
                    help="Disable PE adapter auto-detection.")
    qc.add_argument("--no-poly-g", action="store_true",
                    help="Disable poly-G trimming (NextSeq 1000 specific).")
    qc.add_argument("--poly-g-min-len", type=int, default=10,
                    help="Minimum poly-G run length to trim.")
    qc.add_argument("--cut-right", action="store_true",
                    help="Enable 3' sliding-window quality trimming. OFF by "
                         "default: it can distort allele fractions at read "
                         "ends in low-VAF somatic calling.")
    qc.add_argument("--cut-right-mean-quality", type=int, default=20,
                    help="Sliding-window mean quality, if --cut-right is set.")
    qc.add_argument("--correction", action="store_true",
                    help="Enable overlap-based base correction. OFF by "
                         "default: it interferes with UMI consensus calling.")
    qc.add_argument("--umi-loc", default=None,
                    choices=["index1", "index2", "read1", "read2",
                             "per_index", "per_read"],
                    help="UMI location for duplex/consensus panels.")
    qc.add_argument("--umi-len", type=int, default=None,
                    help="UMI length in bases.")
    qc.add_argument("--umi-skip", type=int, default=None,
                    help="Bases to skip after the UMI (e.g. a linker).")
    qc.add_argument("--failed-out", action="store_true",
                    help="Write discarded reads to a per-sample file for "
                         "audit.")
    qc.add_argument("--no-overrepresentation", action="store_true",
                    help="Disable fastp overrepresented-sequence analysis.")

    # --- Alignment options ---
    align = parser.add_argument_group("alignment options")
    align.add_argument("--two-pass", action="store_true",
                       help="Enable BWA-MEM2 two-pass alignment.")
    align.add_argument("--skip-bqsr", action="store_true",
                       help="Skip Base Quality Score Recalibration.")
    align.add_argument("--aligned-dir", default=None,
                       help="Directory containing pre-aligned, sorted+indexed "
                            "BAMs named <sample>.sorted.bam. Used together "
                            "with '--skip-steps align' to reuse BAMs produced "
                            "by a separate align_reads.py run instead of "
                            "re-aligning. Defaults to <output-dir>/aligned.")

    # --- PCGR clinical reporting ---
    # PCGR runs its own VEP annotation and is entirely optional: when the
    # executable or the reference bundle is absent the step reports
    # "skipped" and the pipeline continues, like SnpEff and bcftools.
    pcgr = parser.add_argument_group("PCGR clinical report")
    pcgr.add_argument("--pcgr-refdata-dir", default=None,
                      help="PCGR reference data bundle. Supplying it "
                           "enables the PCGR step (actionability tiers, "
                           "TMB, MSI, mutational signatures). Downloaded "
                           "separately; must match the installed PCGR "
                           "version and genome build.")
    pcgr.add_argument("--vep-dir", default=None,
                      help="Ensembl VEP cache directory used by PCGR.")
    pcgr.add_argument("--pcgr-estimate-tmb", action="store_true",
                      help="Ask PCGR to compute tumour mutational burden. "
                           "PCGR leaves it OFF by default and says nothing "
                           "when it is off -- the report just has no TMB "
                           "section, which reads as 'nothing found' rather "
                           "than 'never calculated'. Pair with "
                           "--pcgr-target-size-mb on a panel.")
    pcgr.add_argument("--pcgr-estimate-msi", action="store_true",
                      help="Ask PCGR to predict microsatellite instability. "
                           "OFF by default, and PCGR honours it ONLY for "
                           "WGS/WES tumour-control runs: on a TARGETED or "
                           "tumour-only query it warns and skips, so this "
                           "flag does nothing on a panel.")
    pcgr.add_argument("--pcgr-estimate-signatures", action="store_true",
                      help="Ask PCGR to fit COSMIC mutational signatures "
                           "(SBS3/HRD, MMR, APOBEC). OFF by default. Wants "
                           "a few hundred SNVs to be trustworthy.")
    pcgr.add_argument("--pcgr-target-size-mb", type=float, default=None,
                      help="Megabases actually sequenced, used as PCGR's TMB "
                           "denominator. SET THIS FOR ANY PANEL: PCGR "
                           "defaults to 34 Mb for a TARGETED assay, which is "
                           "exome-sized, so a 4.5 Mb panel reports a TMB "
                           "around 7.5x too low. Use the panel's target "
                           "footprint -- the same BED you pass to "
                           "--intervals.")
    pcgr.add_argument("--pcgr-assay", default="TARGETED",
                      choices=["WGS", "WES", "TARGETED"],
                      help="Assay type reported to PCGR; drives TMB "
                           "scaling.")
    pcgr.add_argument("--pcgr-tumour-site", type=int, default=None,
                      help="PCGR tumour-site code (0 = unspecified). "
                           "Site-specific actionability depends on it.")
    pcgr.add_argument("--pcgr-lift-tags", action="store_true",
                      help="Before running PCGR, copy the tumour's (and "
                           "normal's) FORMAT DP/AF into INFO tags and point "
                           "PCGR at them. Mutect2 writes those as FORMAT "
                           "fields only, so without this PCGR has nothing to "
                           "filter on and its TMB is computed over "
                           "unfiltered calls. Writes a new VCF into the pcgr "
                           "directory; the called VCF is left untouched.")
    pcgr.add_argument("--pcgr-tumor-dp-tag", default=None,
                      help="INFO tag holding tumour depth, for PCGR's TMB "
                           "filtering. NOTE: Mutect2 writes DP/AD/AF as "
                           "FORMAT fields, not INFO, so this is normally "
                           "unset and TMB is then computed unfiltered.")
    pcgr.add_argument("--pcgr-tumor-af-tag", default=None,
                      help="INFO tag holding tumour allele fraction. Same "
                           "caveat as --pcgr-tumor-dp-tag.")
    pcgr.add_argument("--pcgr-legacy-v1", action="store_true",
                      help="Emit the PCGR 1.x bundle flag (--pcgr_dir) "
                           "instead of the 2.x --refdata_dir.")
    pcgr.add_argument("--pcgr-extra-args", default=None,
                      help="Extra arguments forwarded verbatim to pcgr, as "
                           "one quoted string.")

    # --- Runtime options ---
    runtime = parser.add_argument_group("runtime options")
    runtime.add_argument("--threads", type=int, default=8,
                         help="Number of CPU threads.")
    runtime.add_argument("--snpeff-memory", default="8g",
                         help="JVM max heap for SnpEff, as a -Xmx value "
                              "(e.g. '8g', '16g'). The bioconda launcher "
                              "defaults to 1g, which is too small for hg38 "
                              "and fails partway through with an "
                              "OutOfMemoryError, so a value is always "
                              "passed explicitly.")
    runtime.add_argument("--reads-to-process", type=int, default=None,
                         help="Only process this many reads from each FASTQ. "
                              "Useful for quick test runs / debugging.")
    runtime.add_argument("--skip-steps", nargs="*", default=[],
                         help="Skip pipeline stages: qc, align, dedup, bqsr, "
                              "mutect2, contamination, msi, filter, cosmic, "
                              "annotate, coverage, pcgr.")
    runtime.add_argument("--resume", action="store_true",
                         help="Reuse pipeline steps whose output already "
                              "exists in --output-dir, instead of redoing "
                              "them. Opt-in (unlike fastq_qc_clean.py and "
                              "align_reads.py, which skip finished work by "
                              "default) because a stage-3 artefact depends "
                              "on the reference and the known-sites "
                              "resources: reusing one built with different "
                              "settings would silently mix two analyses. "
                              "Refuses to run if the previous run's "
                              "parameters differ.")
    runtime.add_argument("--resume-force", action="store_true",
                         help="Resume even though the previous run's "
                              "parameters differ. The result will mix "
                              "settings from both runs; you are asserting "
                              "that is safe.")
    runtime.add_argument("--dry-run", action="store_true",
                         help="Show commands without executing.")
    return parser


def load_qc_manifest(manifest_path):
    """
    Load a fastq_qc_clean.py run manifest and extract cleaned FASTQ paths.

    THIS IS THE KEY INTERFACING FUNCTION: it lets this script consume the
    output of fastq_qc_clean.py directly, so you can run QC once (with its
    superior QC options) and then feed the cleaned reads here without
    re-running fastp.

    The fastq_qc_clean.py manifest structure per sample record:
        {
            "sample": "...",            # sample name (stem)
            "outputs": {
                "r1": "path/to/{sample}_R1.clean.fastq.gz",
                "r2": "path/to/{sample}_R2.clean.fastq.gz",
            },
            "lanes_merged": bool,
            "inputs": [{"r1": ..., "r2": ...}],
            "status": "ok" | "skipped_existing" | "failed",
            "metrics": {...},
        }

    NOTE on "skipped_existing": when fastq_qc_clean.py found the cleaned
    outputs already present it sets status="skipped_existing" but does NOT
    add an "outputs" key. Those files DO exist though, at the standard
    location <qc_output>/cleaned_fastq/{sample}_R1.clean.fastq.gz. We
    reconstruct that path from the manifest's own output_dir rather than
    dropping an otherwise-good sample.

    RETURNS:
        dict: {sample_name: {"r1": path, "r2": path, "lanes_merged": bool}}
    """
    with open(manifest_path, encoding="utf-8") as fh:
        data = json.load(fh)

    qc_output_dir = data.get("output_dir") or os.path.dirname(manifest_path)
    samples = {}
    for record in data.get("samples", []):
        if record.get("status") == "failed":
            # A failed sample must not be fed forward.
            print(f"[WARN] Sample '{record.get('sample')}' failed in QC; "
                  f"excluding from downstream.")
            continue
        name = record.get("sample")
        outputs = record.get("outputs", {})
        if outputs.get("r1") and outputs.get("r2"):
            r1, r2 = outputs["r1"], outputs["r2"]
        else:
            # "skipped_existing" (or a non-ok, non-failed status): the
            # cleaned files exist at the standard location.
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


# =============================================================================
# SECTION 13c: PANEL PROFILES
# =============================================================================
# Two things happen here, and they happen at different moments in the run.
#
#   1. configure_from_panel() runs immediately after parsing, because a
#      profile changes what validate_args() is validating.
#   2. finalise_panel_targets() runs after the reference is resolved,
#      because both of the things it does -- measuring the BED and checking
#      its contig naming -- are questions about the BED AND the genome
#      together, and the genome may still have been a name at parse time.
#
# panel_profiles.py is imported lazily and its absence is tolerated, the
# same way pcgr_report.py and coverage_report.py are: this file must keep
# working when only it was copied somewhere. The difference is that asking
# for --panel when the module is missing is FATAL rather than skipped --
# the operator asked for a configuration and would otherwise get a
# differently configured run that said nothing about it.


def _load_panel_module(required, parser=None):
    """Import panel_profiles.py from beside this script, or None."""
    try:
        import panel_profiles
        return panel_profiles
    except ImportError as exc:
        if required:
            message = (f"--panel needs panel_profiles.py beside this script, "
                       f"and it could not be imported ({exc}). Copy it into "
                       f"the same directory, or configure the run with "
                       f"explicit flags.")
            if parser:
                parser.error(message)
            print(f"[ERROR] {message}")
            sys.exit(1)
        return None


def configure_from_panel(args, parser, argv):
    """
    Apply --panel-bed and --panel to the parsed arguments, in place.

    Returns (record, explicit): the record goes into the run manifest so a
    later reader can see which assay the numbers belong to, and `explicit`
    is the set of dests the command line named, which the caller needs
    again later for the same precedence rule.

    ORDER MATTERS. --panel-bed is applied BEFORE the profile, so that a
    profile which sets coverage thresholds is reasoning about a run that
    already has a BED.
    """
    record = {"panel": None, "applied": [], "overridden": [],
              "panel_bed": None, "notes": []}

    wants_panel = bool(args.panel or args.save_panel_as)
    panel_profiles = _load_panel_module(required=wants_panel, parser=parser)
    explicit = set()
    if panel_profiles is not None:
        explicit = set(panel_profiles.explicit_dests(parser, argv))

    # --- one BED, two flags -------------------------------------------
    # --intervals and --coverage-bed are separate because they answer
    # different questions, and a laboratory occasionally reports coverage
    # over a narrower subset than it calls over. In every other case they
    # are the same file, and requiring it twice means one of them gets
    # forgotten -- which is either a genome-wide call set or a run that
    # makes no coverage statement, both silent.
    if args.panel_bed:
        if not os.path.isfile(args.panel_bed):
            parser.error(f"--panel-bed not found: {args.panel_bed}")
        record["panel_bed"] = os.path.abspath(args.panel_bed)
        for dest in ("intervals", "coverage_bed"):
            # Both default to None, so a value already here can only have
            # come from the command line. Testing the value rather than
            # consulting `explicit` keeps this right even when
            # panel_profiles.py could not be imported and `explicit` is
            # therefore empty.
            if getattr(args, dest):
                record["notes"].append(
                    f"--{dest.replace('_', '-')} was given explicitly, so "
                    f"--panel-bed did not set it.")
                continue
            setattr(args, dest, args.panel_bed)
            # Counts as explicit from here on. A site profile may name a
            # BED of its own -- the schema allows it, for a laboratory
            # whose design file lives at a fixed path -- and without this
            # it would overwrite the one the operator just named on the
            # command line, which is the opposite of the precedence rule
            # everything else here follows.
            explicit.add(dest)

    if not args.panel:
        return record, explicit

    # --- the profile ---------------------------------------------------
    try:
        registry = panel_profiles.available_panels(
            extra_files=args.panel_file,
            warn=lambda msg: print(f"[WARN] {msg}"))
        profile = panel_profiles.resolve_panel(args.panel, registry)
        applied, overridden = panel_profiles.apply_profile(
            args, profile, explicit)
    except panel_profiles.PanelError as exc:
        parser.error(str(exc))
        return record, explicit           # unreachable; keeps linters quiet

    print()
    print(panel_profiles.format_application(profile, applied, overridden))

    # Recorded as the profile's own dict rather than as a name, so the
    # manifest still describes the run after somebody edits or deletes the
    # profile it was configured from. A run record that says only
    # "--panel our-lung-panel" is worthless a year later.
    record["panel"] = dict(profile)
    record["applied"] = [[k, v] for k, v in applied]
    record["overridden"] = [[k, v] for k, v in overridden]
    return record, explicit


def finalise_panel_targets(args, ref, dirs, record, explicit):
    """
    Make the BED usable against this reference, and measure it.

    THE CONTIG PROBLEM. Vendor BEDs arrive in whichever naming their
    designer used: Ensembl writes "1" and "MT", UCSC and the hg38 this
    pipeline installs write "chr1" and "chrM". A mismatched BED makes GATK
    abort, which is at least loud -- but the SAME file handed to the
    coverage check simply finds no reads in any region, and the report then
    states that none of the panel was covered. That is a confident,
    entirely wrong clinical statement, so the mismatch is detected and a
    renamed copy written beside the run.

    THE DENOMINATOR PROBLEM. PCGR divides by an assumed 34 Mb for a
    TARGETED assay unless told otherwise, so a 2 Mb panel reports a TMB
    about 17x too low and says nothing about having done so. The footprint
    is therefore measured from the BED actually being used -- which beats
    any published figure, because it accounts for the kit version, the
    build and whatever the laboratory spiked in.
    """
    panel_profiles = _load_panel_module(required=False)
    if panel_profiles is None:
        if args.intervals or args.coverage_bed:
            print("[WARN] panel_profiles.py not importable: the target BED "
                  "was not checked against the reference's contig naming, "
                  "and no TMB denominator was measured from it.")
        return

    # --- contig naming -------------------------------------------------
    # Cached by source path: --intervals and --coverage-bed are usually the
    # same file, and translating it twice would write the same copy twice
    # and print the same paragraph twice, which reads like two separate
    # problems.
    translated = {}
    for dest in ("intervals", "coverage_bed"):
        bed = getattr(args, dest, None)
        if not bed:
            continue
        if bed not in translated:
            fixed, note = panel_profiles.harmonise_bed_contigs(
                bed, ref, dirs["panel"], dry_run=args.dry_run)
            translated[bed] = fixed
            if note:
                print(f"[INFO] {note}")
                record["notes"].append(note)
        if translated[bed] != bed:
            setattr(args, dest, translated[bed])

    # --- footprint ------------------------------------------------------
    # Measured from --intervals when there is one, because that is the
    # region variants can actually be called in; --coverage-bed is the
    # fallback for a run that reports coverage without restricting the
    # caller.
    bed = args.intervals or args.coverage_bed
    if not bed:
        return
    if not os.path.isfile(bed):
        # A renamed copy that a dry run did not actually write. The base
        # count is a property of the regions, not of what they are called,
        # so fall back to the file the operator supplied rather than
        # reporting no footprint at all -- a dry run is where somebody
        # checks the denominator before committing hours to it.
        bed = record.get("panel_bed") or next(
            (k for k, v in translated.items() if v == bed), bed)
    measured = panel_profiles.footprint_mb(bed)
    if measured is None:
        print(f"[WARN] could not measure a target footprint from "
              f"{os.path.basename(bed)}; TMB, if requested, will use PCGR's "
              f"assumed default.")
        return

    record["target_size_mb_measured"] = measured
    nominal = (record.get("panel") or {}).get("nominal_target_size_mb")
    note = panel_profiles.target_size_note(measured, nominal)
    if note:
        print(f"[INFO] {note}")
        record["notes"].append(note)

    # Only filled in when the command line did not state it. Somebody who
    # passes --pcgr-target-size-mb has a reason -- a callable-bases figure
    # from their validation, say -- and measuring the raw BED would quietly
    # replace a considered number with a cruder one.
    if args.pcgr_target_size_mb is None and "pcgr_target_size_mb" \
            not in explicit:
        args.pcgr_target_size_mb = measured
        print(f"[INFO] TMB denominator set to the measured footprint "
              f"({measured:.3f} Mb). Pass --pcgr-target-size-mb to override "
              f"it, e.g. with a callable-bases figure from your validation.")
        record["notes"].append(
            f"TMB denominator taken from the target BED: {measured:.3f} Mb.")

    advice = panel_profiles.tmb_advice(measured, args.pcgr_estimate_tmb)
    if advice:
        print(f"[WARN] {advice}")
        record["notes"].append(advice)


def save_panel_profile(args, record, explicit=frozenset()):
    """
    Write the run's effective settings out as a reusable profile.

    The other half of customisation: a laboratory bringing up a new kit
    tunes it with ordinary flags under --dry-run, and this turns the
    command line they arrived at into something the next operator can name.
    """
    panel_profiles = _load_panel_module(required=True)
    stem = os.path.splitext(os.path.basename(args.save_panel_as))[0]
    base = (record.get("panel") or {})
    panel_id = args.panel_id or stem
    notes = None
    if base:
        notes = [
            f"Derived from the profile '{base.get('id')}' "
            f"({base.get('name')}), with this run's settings applied on top. "
            f"Review it before using it on another assay.",
            "Target BEDs are not stored in a profile: pass --panel-bed when "
            "you use it.",
        ]
    profile = panel_profiles.profile_from_args(
        args,
        panel_id=panel_id,
        # Named after the id rather than after the profile it was derived
        # from: a file called mylab-spe.json whose name says "QIAseq
        # Targeted DNA Panel" is the sort of thing that gets mistaken for
        # the vendor's own profile a year later.
        name=panel_id,
        manufacturer=base.get("manufacturer") or "local",
        chemistry=base.get("chemistry") or "hybrid-capture",
        notes=notes,
        # Keep the denominator only if the operator stated it themselves;
        # one measured from this run's BED describes that design, not the
        # assay.
        keep_target_size="pcgr_target_size_mb" in explicit,
    )
    try:
        path = panel_profiles.save_profile(profile, args.save_panel_as)
    except (panel_profiles.PanelError, OSError) as exc:
        print(f"[WARN] could not save the panel profile: {exc}")
        return
    print(f"[INFO] Saved this run's settings as panel profile "
          f"'{profile['id']}' in {path}. Use it with "
          f"--panel-file {path} --panel {profile['id']}.")


def validate_args(args, parser):
    """Validate argument combinations that argparse can't catch."""
    if args.threads < 1:
        parser.error("--threads must be >= 1.")

    # Normal requires all three normal args.
    if bool(args.normal_sample) != bool(args.normal_r1 and args.normal_r2):
        parser.error("If --normal-sample is provided, you must also provide "
                     "--normal-r1 and --normal-r2 (and vice versa).")

    # UMI length is only meaningful when the UMI sits INSIDE a read: for
    # index1/index2/per_index the whole index read is the UMI, so fastp
    # needs no length. Requiring it for those would reject a command that
    # fastq_qc_clean.py accepts -- keep this predicate identical to the one
    # in fastq_qc_clean.validate_args().
    if args.umi_loc in ("read1", "read2", "per_read") and not args.umi_len:
        parser.error("--umi-len is required when --umi-loc is "
                     "read1/read2/per_read.")
    if args.umi_skip is not None and not args.umi_loc:
        parser.error("--umi-skip given without --umi-loc.")

    # A PCGR bundle with no VEP cache produces a run that completes and then
    # cannot report: PCGR requires --vep_dir whenever it is handed an input
    # VCF, and it only says so once step 13 starts -- hours in. Refuse the
    # combination here, where it costs nothing.
    if args.pcgr_refdata_dir and not args.vep_dir:
        parser.error("--pcgr-refdata-dir needs --vep-dir: PCGR annotates "
                     "with Ensembl VEP and refuses an input VCF without a "
                     "cache, so the report would fail after the run. The "
                     "installer puts one in ~/data/vep_cache. Omit "
                     "--pcgr-refdata-dir to run without a clinical report.")

    # Manifest and auto-discover are mutually exclusive input modes.
    if args.manifest and args.auto_discover:
        parser.error("--manifest and --auto-discover are mutually exclusive. "
                     "Use one input mode: manifest, auto-discover, or explicit "
                     "sample args.")

    # If a manifest is given but would require sample args, warn.
    if args.manifest:
        if not os.path.isfile(args.manifest):
            parser.error(f"Manifest not found: {args.manifest}")
        # Explicit sample args are allowed but overridden if --tumour-sample
        # matches a sample in the manifest.

    # Explicit-input mode: all three tumour args are required.
    if not args.manifest and not args.auto_discover:
        if not (args.tumour_sample and args.tumour_r1 and args.tumour_r2):
            parser.error("Please provide --tumour-sample, --tumour-r1 and "
                         "--tumour-r2, or use --auto-discover, or reuse a "
                         "previous QC run via --manifest.")

    # --input-dir is required in every mode except --manifest.
    if not args.manifest:
        if not args.input_dir:
            parser.error("--input-dir is required (unless --manifest is "
                         "given).")

    # Resolve FASTQ paths (try input-dir if not found).
    if not args.auto_discover and not args.manifest:
        for role in ("tumour_r1", "tumour_r2", "normal_r1", "normal_r2"):
            val = getattr(args, role)
            if val is not None and not os.path.isfile(val):
                candidate = os.path.join(args.input_dir, os.path.basename(val))
                if os.path.isfile(candidate):
                    setattr(args, role, candidate)
                else:
                    parser.error(f"FASTQ not found: {val} (also tried {candidate})")

    # Validate COSMIC path.
    if args.cosmic and not os.path.isfile(args.cosmic):
        parser.error(f"COSMIC VCF not found: {args.cosmic}. "
                     "Download from https://cancer.sanger.ac.uk/cosmic")


# =============================================================================
# SECTION 15: MAIN PIPELINE ORCHESTRATION
# =============================================================================

def main():
    """
    Main entry point: parse arguments, validate inputs, run the
    10-step comprehensive variant calling pipeline.
    """
    # --list-panels, --describe-panel and --new-panel-template answer a
    # question that has nothing to do with a run, so they must work without
    # --output-dir and --reference. argparse cannot express "these required
    # arguments are not required today", so they are handled from sys.argv
    # before the real parser runs. They are still declared on the parser, so
    # they appear in --help and a misspelling is still rejected.
    panel_module = _load_panel_module(required=False)
    if panel_module is not None:
        panel_module.handle_panel_queries()

    parser = build_parser()
    args = parser.parse_args()

    # Before validate_args: a profile changes the values being validated.
    panel_record, panel_explicit = configure_from_panel(
        args, parser, sys.argv[1:])

    validate_args(args, parser)

    started = datetime.now(timezone.utc)
    skip_steps = set(args.skip_steps)

    # ---- Pre-flight tool checks ----
    # fastp is only required when QC actually runs (it is auto-skipped in
    # manifest mode). Tools are not required at all for --dry-run, which
    # should preview commands even on machines without them installed.
    required_tools = {
        "bwa-mem2": "bwa-mem2",
        "samtools": "samtools",
        "gatk": "gatk",
    }
    if not args.manifest:
        required_tools["fastp"] = "fastp"
    if not args.dry_run:
        missing = [label for label, exe in required_tools.items()
                   if shutil.which(exe) is None]
        if missing:
            print(f"[ERROR] Required tools not found: {', '.join(missing)}")
            sys.exit(1)

    # ---- Create output directories ----
    output_dir = os.path.abspath(args.output_dir)
    dirs = {
        "cleaned": os.path.join(output_dir, "cleaned_fastq"),
        "align": os.path.abspath(args.aligned_dir) if args.aligned_dir
                 else os.path.join(output_dir, "aligned"),
        "dedup": os.path.join(output_dir, "dedup"),
        "bqsr": os.path.join(output_dir, "bqsr"),
        "mutect2": os.path.join(output_dir, "mutect2"),
        "annotated": os.path.join(output_dir, "annotated"),
        "pcgr": os.path.join(output_dir, "pcgr"),
        "coverage": os.path.join(output_dir, "coverage"),
        # Anything derived from the panel BED -- currently a contig-renamed
        # copy. Kept inside the run so the run is self-describing about
        # what it actually called over, and so the laboratory's own copy of
        # a controlled vendor document is never edited.
        "panel": os.path.join(output_dir, "panel"),
        "metrics": os.path.join(output_dir, "metrics"),
        "stats": os.path.join(output_dir, "stats"),
        "logs": os.path.join(output_dir, "logs"),
        "reports": os.path.join(output_dir, "fastp_reports"),
        "tmp": os.path.join(output_dir, "tmp"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    # ---- Reference genome ----
    ref = ensure_reference(args)
    ensure_indices(ref, args.threads, dry_run=args.dry_run)

    # Now that the genome is a real indexed FASTA, the target BED can be
    # checked against its contig naming and measured for the TMB
    # denominator. Both need the reference, which is why neither happened
    # at parse time.
    finalise_panel_targets(args, ref, dirs, panel_record, panel_explicit)
    if args.save_panel_as:
        save_panel_profile(args, panel_record, panel_explicit)

    # ---- Discover or use specified FASTQs ----
    # THREE INPUT MODES (priority order):
    #   1. --manifest: reuse fastq_qc_clean.py output (cleaned FASTQs).
    #   2. --auto-discover: find FASTQ pairs by filename pattern.
    #   3. Explicit --tumour-r1/r2 (and optional --normal-r1/r2) args.
    # Concatenated lane files created below (auto-discover mode only). They
    # are copies of the input FASTQs, so they can be very large; we delete
    # them once the run finishes. Tracked here so the cleanup at the end of
    # main() can find them.
    lane_merge_files = []

    manifest_samples = None
    if args.manifest:
        manifest_samples = load_qc_manifest(args.manifest)
        if not manifest_samples:
            print("[ERROR] No usable samples found in manifest. All samples "
                  "either failed QC or have missing outputs.")
            sys.exit(1)

        # Match tumour by user-provided sample name; otherwise take the
        # first manifest entry as tumour and the second as normal.
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
                # Take the next sample as the normal (skip the tumour's key).
                candidates = [k for k in sorted(manifest_samples.keys())
                              if k != tumour["name"]]
                if candidates:
                    nkey = candidates[0]
                    normal = {"name": nkey,
                              "r1": manifest_samples[nkey]["r1"],
                              "r2": manifest_samples[nkey]["r2"]}

        # IMPORTANT: QC was already run by fastq_qc_clean.py. The manifest
        # points at the cleaned reads, so the "qc" step MUST be skipped.
        skip_steps.add("qc")
        print(f"[INFO] Loading {len(manifest_samples)} cleaned sample(s) from "
              f"manifest: {args.manifest}")
        print("[INFO] 'qc' step auto-skipped because cleaned reads are "
              "provided via manifest.")
    elif args.auto_discover:
        pairs, warnings = find_paired_fastqs(
            args.input_dir, args.r1_pattern, args.r2_pattern,
            recursive=args.recursive)
        for w in warnings:
            print(f"[WARN] {w}")
        if not pairs:
            print("[ERROR] No FASTQ pairs found with auto-discover.")
            sys.exit(1)

        def resolve_lanes(pair, label):
            """Return r1/r2 for a pair, concatenating lanes if lane-split.
            Lane-split data always contributes ALL lanes (never only the
            first); --merge-lanes merely suppresses the informational note.

            DRY RUN: the concatenation is skipped -- on real lane-split WGS
            input it would copy tens of GB, which a preview must not do. We
            still return the paths the merged files WOULD occupy so the
            commands printed downstream read correctly."""
            if pair.get("lanes") and len(pair["lanes"]) > 1:
                print(f"[INFO] {len(pair['lanes'])} lanes found for "
                      f"{label} sample '{pair['sample']}'; concatenating "
                      f"all lanes into one sample.")
                r1_out = os.path.join(
                    dirs["tmp"], f"{pair['sample']}_{label}_merged_R1.fastq.gz")
                r2_out = os.path.join(
                    dirs["tmp"], f"{pair['sample']}_{label}_merged_R2.fastq.gz")

                if args.dry_run:
                    print(f"[DRY RUN] Would merge {len(pair['lanes'])} lanes "
                          f"-> {r1_out} / {r2_out}")
                    return r1_out, r2_out

                os.makedirs(dirs["tmp"], exist_ok=True)
                r1_merge = concatenate([ln["r1"] for ln in pair["lanes"]],
                                       r1_out)
                r2_merge = concatenate([ln["r2"] for ln in pair["lanes"]],
                                       r2_out)
                # Remember them so they get cleaned up at the end of the run.
                lane_merge_files.extend([r1_merge, r2_merge])
                return r1_merge, r2_merge
            return pair["r1"], pair["r2"]

        # Assign the first pair as tumour, second as normal.
        t_r1, t_r2 = resolve_lanes(pairs[0], "tumour")
        tumour = {"name": pairs[0]["sample"], "r1": t_r1, "r2": t_r2}
        normal = None
        if len(pairs) > 1:
            n_r1, n_r2 = resolve_lanes(pairs[1], "normal")
            normal = {"name": pairs[1]["sample"], "r1": n_r1, "r2": n_r2}
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

    # ---- Resume: fold already-finished steps into skip_steps ----
    # Done here, once the samples are known but before any step runs, so
    # the rest of main() needs no resume-specific logic: it simply sees
    # more steps in skip_steps. `resumed_steps` is kept separately because
    # a step skipped by --resume must fall back to ITS OWN output, whereas
    # a step named in --skip-steps falls back to the previous stage's.
    resumed_steps = set()
    if args.resume:
        resumed_steps = plan_resume(dirs, tumour, normal, args, output_dir)
        skip_steps |= resumed_steps

    # ---- Print summary ----
    print(f"\n{'=' * 72}")
    print("COMPREHENSIVE CANCER VARIANT CALLING PIPELINE")
    print(f"{'=' * 72}")
    print(f"Tumour  : {tumour['name']}")
    if normal:
        print(f"Normal  : {normal['name']}")
    else:
        print("Mode    : tumour-only (no matched normal)")
    if panel_record.get("panel"):
        print(f"Panel   : {panel_record['panel']['name']} "
              f"[{panel_record['panel']['id']}]")
    elif args.intervals or args.coverage_bed:
        print("Panel   : none named (target BED given directly)")
    else:
        # Worth one line: the commonest configuration mistake on a panel
        # run is to leave the target out entirely, and it is invisible
        # until the VCF turns out to hold tens of thousands of calls.
        print("Panel   : none -- no target restriction (whole genome)")
    print(f"COSMIC  : {args.cosmic or 'not provided'}")
    print(f"Ref     : {ref}")
    print(f"Threads : {args.threads}")
    print(f"Output  : {output_dir}")
    print(f"Skip    : {', '.join(sorted(skip_steps)) or 'none'}")
    if resumed_steps:
        print(f"Resumed : {', '.join(sorted(resumed_steps))}")
    print(f"{'=' * 72}")

    # ---- Build known sites for BQSR ----
    known_sites = []
    if args.dbsnp:
        known_sites.append(os.path.abspath(args.dbsnp))
    for indel in (args.known_indels or []):
        known_sites.append(os.path.abspath(indel))

    # ---- Run manifest (full pipeline record) ----
    manifest = {
        "script": os.path.basename(__file__),
        "started_utc": started.isoformat(),
        "host": os.uname().nodename if hasattr(os, "uname") else None,
        "python": sys.version.split()[0],
        "command_line": shlex.join(sys.argv),
        "parameters": vars(args),
        "tool_versions": {n: tool_version(e) for n, e in required_tools.items()},
        "reference": ref,
        "tumour": tumour["name"],
        "normal": normal["name"] if normal else None,
        "cosmic": args.cosmic,
        # The assay, as data rather than as a name: the profile's own
        # contents are copied in, so this manifest still describes the run
        # after the profile it came from is edited or deleted.
        "panel_profile": panel_record,
        "steps_completed": [],
    }

    def should_run(step):
        """Check if a pipeline step should run."""
        return step not in skip_steps

    # ================================================================
    # STEP 1: QC & CLEANING
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 1/14: QC & Cleaning (fastp)")
    print(f"{'#' * 72}")

    tumour_r1_clean = os.path.join(dirs["cleaned"],
                                   f"{tumour['name']}_R1.clean.fastq.gz")
    tumour_r2_clean = os.path.join(dirs["cleaned"],
                                   f"{tumour['name']}_R2.clean.fastq.gz")
    normal_r1_clean = None
    normal_r2_clean = None

    if should_run("qc"):
        # Clean tumour reads.
        code = run_fastp(
            tumour["name"], tumour["r1"], tumour["r2"],
            tumour_r1_clean, tumour_r2_clean,
            dirs["reports"], args,
            os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if code != 0:
            print("[FATAL] Tumour QC failed.")
            sys.exit(1)
        manifest["steps_completed"].append("qc_tumour")

        # Clean normal reads if provided.
        if normal:
            normal_r1_clean = os.path.join(
                dirs["cleaned"], f"{normal['name']}_R1.clean.fastq.gz")
            normal_r2_clean = os.path.join(
                dirs["cleaned"], f"{normal['name']}_R2.clean.fastq.gz")
            code = run_fastp(
                normal["name"], normal["r1"], normal["r2"],
                normal_r1_clean, normal_r2_clean,
                dirs["reports"], args,
                os.path.join(dirs["logs"], f"{normal['name']}.log"),
                dry_run=args.dry_run,
            )
            if code != 0:
                print("[FATAL] Normal QC failed.")
                sys.exit(1)
            manifest["steps_completed"].append("qc_normal")

        # Record fastp metrics.
        for name, r1_clean in [(tumour["name"], tumour_r1_clean),
                                (normal["name"] if normal else None, normal_r1_clean)]:
            if name and r1_clean:
                json_path = os.path.join(dirs["reports"], f"{name}.fastp.json")
                if os.path.exists(json_path):
                    manifest[f"fastp_{name}"] = summarise_fastp_json(json_path)
                    # A file that is not FASTQ, an empty one, or a mismatched
                    # pair leaves fastp exiting 0 with nothing to show for
                    # it. Every later step then "succeeds" over no data and
                    # the run reports success with an empty VCF -- a result
                    # indistinguishable from a genuine absence of variants.
                    # Nothing downstream can recover from this, so stop.
                    kept = manifest[f"fastp_{name}"].get("reads_after")
                    if kept == 0 and not args.dry_run:
                        print(f"[FATAL] {name}: fastp kept 0 reads. Nothing "
                              f"downstream can be computed from an empty "
                              f"read set, and a run that continued would "
                              f"report success with an empty VCF.")
                        print(f"        Check the input files are real, "
                              f"paired FASTQ: {tumour['r1'] if name == tumour['name'] else (normal or {}).get('r1')}")
                        sys.exit(1)
    elif "qc" in resumed_steps:
        # Resumed: the cleaned FASTQs this step would have written are
        # already on disk, so keep pointing at THEM (the dirs["cleaned"]
        # paths set above) rather than at the raw inputs.
        print("[SKIP] QC skipped (resumed: cleaned FASTQs already present).")
        if normal:
            normal_r1_clean = os.path.join(
                dirs["cleaned"], f"{normal['name']}_R1.clean.fastq.gz")
            normal_r2_clean = os.path.join(
                dirs["cleaned"], f"{normal['name']}_R2.clean.fastq.gz")
    else:
        print("[SKIP] QC skipped.")
        # Assume the given FASTQs are already clean.
        tumour_r1_clean = tumour["r1"]
        tumour_r2_clean = tumour["r2"]
        if normal:
            normal_r1_clean = normal["r1"]
            normal_r2_clean = normal["r2"]

    # ================================================================
    # STEP 2: ALIGNMENT
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 2/14: Alignment (BWA-MEM2)")
    print(f"{'#' * 72}")

    tumour_bam = os.path.join(dirs["align"], f"{tumour['name']}.sorted.bam")
    normal_bam = None

    if should_run("align"):
        result = align_sample(
            tumour["name"], tumour_r1_clean, tumour_r2_clean,
            ref, tumour_bam, args.threads, args.two_pass,
            os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if result is None:
            print("[FATAL] Tumour alignment failed.")
            sys.exit(1)
        manifest["steps_completed"].append("align_tumour")

        if normal and normal_r1_clean:
            normal_bam = os.path.join(dirs["align"],
                                      f"{normal['name']}.sorted.bam")
            result = align_sample(
                normal["name"], normal_r1_clean, normal_r2_clean,
                ref, normal_bam, args.threads, args.two_pass,
                os.path.join(dirs["logs"], f"{normal['name']}.log"),
                dry_run=args.dry_run,
            )
            if result is None:
                print("[FATAL] Normal alignment failed.")
                sys.exit(1)
            manifest["steps_completed"].append("align_normal")

        # Collect alignment stats (writes real flagstat/idxstats/stats files
        # into dirs["stats"] and parses real metrics -- see collect_stats()).
        for name, bam in [(tumour["name"], tumour_bam),
                          (normal["name"] if normal else None, normal_bam)]:
            if name and bam and os.path.exists(bam):
                alignment_stats = collect_stats(
                    bam, name, dirs["stats"],
                    os.path.join(dirs["logs"], f"{name}.log"),
                    dry_run=args.dry_run,
                )
                manifest[f"alignment_stats_{name}"] = alignment_stats
                # The same trap as an empty QC result, one step later and
                # reachable even when QC was skipped or a manifest supplied
                # the reads: nothing aligned means every later step runs
                # over an empty BAM and the run reports success with an
                # empty VCF. flagstat has already been parsed, so this
                # costs nothing.
                mapped = alignment_stats.get("mapped_reads")
                if not args.dry_run and mapped is not None and \
                        str(mapped).strip() in ("0", "0.0"):
                    print(f"[FATAL] {name}: 0 reads aligned to the "
                          f"reference. Every later step would run over an "
                          f"empty BAM and the run would report success with "
                          f"an empty VCF.")
                    print("        Check that the reads and the reference "
                          "belong together (species, build and contig "
                          "naming), and that the input really is FASTQ.")
                    sys.exit(1)
    else:
        print(f"[SKIP] Alignment skipped; expecting an existing BAM at "
              f"{tumour_bam}")

        def require_existing_bam(path):
            """Abort unless a usable pre-aligned BAM is sitting at `path`.

            Both the BAM and its .bai must be present: every GATK tool
            downstream (MarkDuplicates, BaseRecalibrator, Mutect2) needs
            random access, and a BAM without an index fails deep inside
            GATK with a far less obvious error than this one.
            align_reads.py applies the same bam+bai rule when it decides
            whether a sample can be skipped."""
            if args.dry_run:
                return  # Nothing exists yet in a preview; that is expected.
            if not os.path.exists(path):
                print(f"[FATAL] '--skip-steps align' requires an existing "
                      f"BAM at {path} (use --aligned-dir to point "
                      f"elsewhere).")
                sys.exit(1)
            if not os.path.exists(path + ".bai"):
                print(f"[FATAL] BAM index missing: {path}.bai -- GATK needs "
                      f"it. Run: samtools index {path}")
                sys.exit(1)

        require_existing_bam(tumour_bam)
        if normal:
            normal_bam = os.path.join(dirs["align"],
                                      f"{normal['name']}.sorted.bam")
            require_existing_bam(normal_bam)

    # ================================================================
    # STEP 3: DUPLICATE MARKING
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 3/14: Duplicate Marking (GATK)")
    print(f"{'#' * 72}")

    tumour_dedup = os.path.join(dirs["dedup"],
                                f"{tumour['name']}.dedup.bam")
    normal_dedup = None

    if should_run("dedup"):
        ok = mark_duplicates(
            tumour_bam, tumour_dedup,
            os.path.join(dirs["metrics"],
                         f"{tumour['name']}.dup_metrics.txt"),
            dirs["tmp"], True, args.threads,
            os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if not ok:
            print("[FATAL] Tumour dedup failed.")
            sys.exit(1)
        manifest["steps_completed"].append("dedup_tumour")

        if normal and normal_bam:
            normal_dedup = os.path.join(dirs["dedup"],
                                        f"{normal['name']}.dedup.bam")
            ok = mark_duplicates(
                normal_bam, normal_dedup,
                os.path.join(dirs["metrics"],
                             f"{normal['name']}.dup_metrics.txt"),
                dirs["tmp"], True, args.threads,
                os.path.join(dirs["logs"], f"{normal['name']}.log"),
                dry_run=args.dry_run,
            )
            if not ok:
                print("[FATAL] Normal dedup failed.")
                sys.exit(1)
            manifest["steps_completed"].append("dedup_normal")
    elif "dedup" in resumed_steps:
        # Resumed: use the duplicate-marked BAMs already on disk.
        print("[SKIP] Dedup skipped (resumed: dedup BAM already present).")
        if normal:
            normal_dedup = os.path.join(dirs["dedup"],
                                        f"{normal['name']}.dedup.bam")
    else:
        print("[SKIP] Dedup skipped.")
        tumour_dedup = tumour_bam
        if normal_bam:
            normal_dedup = normal_bam

    # ================================================================
    # STEP 4: BQSR
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 4/14: Base Quality Score Recalibration (BQSR)")
    print(f"{'#' * 72}")

    tumour_bqsr = os.path.join(dirs["bqsr"],
                               f"{tumour['name']}.bqsr.bam")
    normal_bqsr = None

    if should_run("bqsr") and not args.skip_bqsr and known_sites:
        ok = run_bqsr(
            tumour_dedup, ref, known_sites, tumour_bqsr,
            os.path.join(dirs["metrics"], f"{tumour['name']}.bqsr.table"),
            dirs["tmp"],
            os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
            intervals=args.intervals,
            interval_padding=args.interval_padding,
        )
        if not ok:
            print("[FATAL] Tumour BQSR failed.")
            sys.exit(1)
        manifest["steps_completed"].append("bqsr_tumour")

        if normal and normal_dedup:
            normal_bqsr = os.path.join(dirs["bqsr"],
                                       f"{normal['name']}.bqsr.bam")
            ok = run_bqsr(
                normal_dedup, ref, known_sites, normal_bqsr,
                os.path.join(dirs["metrics"],
                             f"{normal['name']}.bqsr.table"),
                dirs["tmp"],
                os.path.join(dirs["logs"], f"{normal['name']}.log"),
                dry_run=args.dry_run,
                intervals=args.intervals,
                interval_padding=args.interval_padding,
            )
            if not ok:
                print("[FATAL] Normal BQSR failed.")
                sys.exit(1)
            manifest["steps_completed"].append("bqsr_normal")
    elif "bqsr" in resumed_steps:
        # Resumed: use the recalibrated BAMs already on disk.
        print("[SKIP] BQSR skipped (resumed: recalibrated BAM already "
              "present).")
        if normal:
            normal_bqsr = os.path.join(dirs["bqsr"],
                                       f"{normal['name']}.bqsr.bam")
    else:
        if args.skip_bqsr:
            print("[SKIP] BQSR skipped (--skip-bqsr).")
        elif not known_sites:
            print("[SKIP] BQSR skipped (no --dbsnp or --known-indels).")
        else:
            print("[SKIP] BQSR skipped.")
        tumour_bqsr = tumour_dedup
        if normal_dedup:
            normal_bqsr = normal_dedup

    # ================================================================
    # STEP 5: MUTECT2 VARIANT CALLING
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 5/14: Somatic Variant Calling (GATK Mutect2)")
    print(f"{'#' * 72}")

    mutect2_vcf = os.path.join(dirs["mutect2"],
                               f"{tumour['name']}.mutect2.vcf.gz")
    f1r2_dir = os.path.join(dirs["mutect2"], "f1r2")

    if should_run("mutect2"):
        ok = run_mutect2(
            tumour_bqsr, tumour["name"], ref, mutect2_vcf,
            normal_bam=normal_bqsr,
            normal_name=normal["name"] if normal else None,
            germline_resource=args.germline_resource,
            panel_of_normals=args.panel_of_normals,
            f1r2_dir=f1r2_dir,
            tmp_dir=dirs["tmp"],
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
            intervals=args.intervals,
            interval_padding=args.interval_padding,
            extra_args=(shlex.split(args.mutect2_extra_args)
                        if args.mutect2_extra_args else None),
        )
        if not ok:
            print("[FATAL] Mutect2 failed.")
            sys.exit(1)
        manifest["steps_completed"].append("mutect2")
        manifest["mutect2_vcf"] = mutect2_vcf
    elif "mutect2" in resumed_steps:
        print("[SKIP] Mutect2 skipped (resumed: raw VCF already present).")
        manifest["mutect2_vcf"] = mutect2_vcf
    else:
        print("[SKIP] Mutect2 skipped.")

    # ================================================================
    # STEP 6: STRAND BIAS MODELLING
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 6/14: Strand Bias Modelling (LearnReadOrientation)")
    print(f"{'#' * 72}")

    orientation_model = os.path.join(dirs["mutect2"],
                                     "read_orientation_model.tar.gz")

    if should_run("filter"):
        f1r2_tar = os.path.join(f1r2_dir,
                                f"{tumour['name']}.f1r2.tar.gz")
        if os.path.exists(f1r2_tar) or args.dry_run:
            run_learn_read_orientation(
                f1r2_tar, orientation_model,
                log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
                dry_run=args.dry_run,
            )
        else:
            print("[WARN] F1R2 tarball not found; skipping orientation model.")

    # ================================================================
    # STEP 7: CONTAMINATION ESTIMATION
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 7/14: Contamination Estimation (GATK)")
    print(f"{'#' * 72}")

    # Enrichment, not a dependency: contamination_table stays None unless
    # the step both runs and succeeds, and filter_mutect2() omits the
    # arguments when it is None. Every failure path below therefore lands
    # on the pre-existing behaviour rather than stopping the run.
    contamination_table = None
    segmentation_table = None
    contam_paths = (
        os.path.join(dirs["metrics"], f"{tumour['name']}.pileups.table"),
        os.path.join(dirs["metrics"], f"{tumour['name']}.contamination.table"),
        os.path.join(dirs["metrics"], f"{tumour['name']}.segments.table"),
    )

    if should_run("contamination") and args.contamination_resource:
        contamination_table, segmentation_table = estimate_contamination(
            tumour_bqsr,
            args.contamination_resource,
            contam_paths[0], contam_paths[1], contam_paths[2],
            tmp_dir=dirs["tmp"],
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if contamination_table:
            manifest["steps_completed"].append("contamination")
            manifest["contamination_table"] = contamination_table
            manifest["tumour_segmentation"] = segmentation_table
    elif "contamination" in resumed_steps:
        # Resumed: reuse this step's OWN artefacts (see --resume vs
        # --skip-steps in the module docstring), but only if both are
        # really on disk -- a half-written pair must not reach GATK.
        if os.path.exists(contam_paths[1]) and os.path.exists(contam_paths[2]):
            contamination_table, segmentation_table = contam_paths[1], contam_paths[2]
            print("[SKIP] Contamination estimation skipped (resumed: "
                  "tables already present).")
            manifest["contamination_table"] = contamination_table
            manifest["tumour_segmentation"] = segmentation_table
        else:
            print("[WARN] Contamination tables missing despite resume; "
                  "filtering without them.")
    elif not args.contamination_resource:
        print("[SKIP] Contamination estimation skipped "
              "(--contamination-resource not provided).")
        print("       FilterMutectCalls will assume zero contamination "
              "rather than measure it, so real contamination in this "
              "sample would go undetected.")
    else:
        print("[SKIP] Contamination estimation skipped "
              "(--skip-steps contamination).")

    # ================================================================
    # STEP 8: MICROSATELLITE INSTABILITY
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 8/14: Microsatellite Instability (MSIsensor2)")
    print(f"{'#' * 72}")

    msi_result = None
    if should_run("msi") and args.msi_models:
        msi_result = run_msisensor2(
            tumour_bqsr, args.msi_models,
            os.path.join(dirs["metrics"], f"{tumour['name']}.msi"),
            intervals=args.intervals, threads=args.threads,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if msi_result:
            manifest["steps_completed"].append("msi")
            manifest["msi"] = msi_result
    elif not args.msi_models:
        print("[SKIP] MSI skipped (--msi-models not provided).")
        print("       PCGR does not cover this: it restricts MSI to WGS/WES")
        print("       tumour-control runs and omits the section otherwise.")
    else:
        print("[SKIP] MSI skipped (--skip-steps msi).")

    # ================================================================
    # STEP 9: MUTECT2 CALL FILTERING
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 9/14: Filtering Mutect2 Calls")
    print(f"{'#' * 72}")

    filtered_vcf = os.path.join(dirs["mutect2"],
                                f"{tumour['name']}.mutect2.filtered.vcf.gz")

    if should_run("filter"):
        ok = filter_mutect2(
            mutect2_vcf, ref, filtered_vcf,
            orientation_model=orientation_model
                if os.path.exists(orientation_model) else None,
            contamination_table=contamination_table,
            segmentation_table=segmentation_table,
            min_allele_fraction=args.min_allele_fraction,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if not ok:
            print("[FATAL] FilterMutectCalls failed.")
            sys.exit(1)
        manifest["steps_completed"].append("filter")
        manifest["filtered_vcf"] = filtered_vcf

        # Part of filtering, not a step of its own: it refines the same
        # artefact in place, so it must not be separately skippable or
        # the FILTER column would depend on which flags a resume used.
        apply_depth_floor(
            filtered_vcf, tumour["name"], args.min_depth,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        manifest["min_depth"] = args.min_depth

        # Also part of filtering: a sequence-structure test that no
        # likelihood model can make, since the reads supporting a foldback
        # are real and map cleanly.
        tag_foldback_insertions(
            filtered_vcf, ref, args.foldback_min_match,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        manifest["foldback_min_match"] = args.foldback_min_match
    elif "filter" in resumed_steps:
        # Resumed: the filtered VCF exists, so keep pointing at it rather
        # than falling back to the unfiltered Mutect2 output.
        print("[SKIP] Mutect2 call filtering skipped (resumed: filtered VCF "
              "already present).")
    else:
        print("[SKIP] Mutect2 call filtering skipped.")
        filtered_vcf = mutect2_vcf

    # ================================================================
    # STEP 10: COSMIC ANNOTATION
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 10/14: COSMIC Database Annotation")
    print(f"{'#' * 72}")

    cosmic_vcf = os.path.join(dirs["mutect2"],
                              f"{tumour['name']}.mutect2.cosmic.vcf.gz")
    vcf_for_annotation = filtered_vcf

    if should_run("cosmic") and args.cosmic:
        ok = annotate_cosmic(
            filtered_vcf, cosmic_vcf, args.cosmic,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if ok is None:
            print("[SKIP] COSMIC annotation skipped "
                  "(bcftools unavailable or COSMIC VCF missing).")
        elif ok:
            manifest["steps_completed"].append("cosmic")
            manifest["cosmic_vcf"] = cosmic_vcf
            vcf_for_annotation = cosmic_vcf

            # Count COSMIC overlaps for the manifest.
            cosmic_counts = count_cosmic_overlaps(
                cosmic_vcf, args.cosmic, dry_run=args.dry_run)
            manifest["cosmic_overlap"] = cosmic_counts
        else:
            print("[WARN] COSMIC annotation failed.")
    elif "cosmic" in resumed_steps:
        # Resumed: hand the already-annotated VCF to SnpEff/PCGR, exactly
        # as a successful COSMIC run would have.
        print("[SKIP] COSMIC skipped (resumed: annotated VCF already "
              "present).")
        vcf_for_annotation = cosmic_vcf
        manifest["cosmic_vcf"] = cosmic_vcf
    else:
        # Several different reasons to land here -- say which, so a silent
        # step never leaves the reader guessing.
        if not args.cosmic:
            print("[SKIP] COSMIC skipped (--cosmic not provided).")
        else:
            print("[SKIP] COSMIC skipped (--skip-steps cosmic).")

    # ================================================================
    # STEP 11: SnpEff ANNOTATION
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 11/14: Gene Annotation (SnpEff)")
    print(f"{'#' * 72}")

    annotated_vcf = os.path.join(dirs["annotated"],
                                 f"{tumour['name']}.annotated.vcf")

    # Auto-detect the genome build from the reference filename. Computed
    # here, OUTSIDE the should_run() guard, because STEP 12 (PCGR) needs it
    # too -- deriving it inside the branch would raise NameError whenever
    # 'annotate' is skipped but PCGR still runs.
    snpeff_genome = "hg38"
    if "hg19" in ref.lower() or "grch37" in ref.lower():
        snpeff_genome = "hg19"

    if should_run("annotate"):
        ok = run_snpeff(
            vcf_for_annotation, annotated_vcf, genome=snpeff_genome,
            memory=args.snpeff_memory,
            log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
            dry_run=args.dry_run,
        )
        if ok is None:
            print("[SKIP] SnpEff annotation skipped (snpEff not installed).")
        elif ok:
            manifest["steps_completed"].append("annotate")
            manifest["annotated_vcf"] = annotated_vcf
        else:
            print("[WARN] SnpEff annotation failed.")
    else:
        print("[SKIP] SnpEff annotation skipped.")

    # ================================================================
    # STEP 12: TARGET COVERAGE CHECK
    # ================================================================
    # Deliberately BEFORE the report: a region the sequencing never covered
    # produces no variant, and a VCF cannot distinguish that from a region
    # that is genuinely wild type. Both reach PCGR as silence. This is the
    # only place in the run where that difference is written down.
    print(f"\n{'#' * 72}")
    print("# STEP 12/14: Target Coverage Check")
    print(f"{'#' * 72}")

    coverage_reports = []
    if not args.coverage_bed:
        print("[SKIP] No --coverage-bed given, so no coverage statement is "
              "made. Absent variants in this run cannot be told apart from "
              "unsequenced regions.")
    elif not should_run("coverage"):
        print("[SKIP] Coverage check skipped (--skip-steps coverage).")
    else:
        try:
            from coverage_report import run_coverage_check
        except ImportError as exc:
            print(f"[WARN] coverage_report.py not importable ({exc}); "
                  f"skipping the coverage check.")
            run_coverage_check = None

        if run_coverage_check is not None:
            cov_dir = dirs["coverage"]
            depth_floor = args.coverage_min_depth or args.min_depth
            # Both samples are checked. A gap in the NORMAL is just as
            # capable of hiding a somatic call: with nothing to compare
            # against, Mutect2 has no evidence the site is somatic.
            targets = [(tumour["name"], tumour_bqsr)]
            if normal and normal_bqsr:
                targets.append((normal["name"], normal_bqsr))

            for name, bam in targets:
                print(f"[INFO] Checking {args.coverage_bed} against {name} "
                      f"at {depth_floor}x ...")
                result = run_coverage_check(
                    bam, args.coverage_bed, name, cov_dir,
                    min_depth=depth_floor,
                    min_mapq=args.coverage_min_mapq,
                    min_baseq=args.coverage_min_baseq,
                    dry_run=args.dry_run)
                if args.dry_run:
                    # measure() has already printed the samtools command it
                    # would run; nothing was measured, so there is nothing
                    # to report and nothing to warn about.
                    print(f"[DRY RUN] Would write {cov_dir}/"
                          f"{name}.coverage.html")
                    continue
                if not result["ok"]:
                    # Never fatal. The calls are already on disk and cost
                    # hours; a missing coverage statement is a caveat on
                    # them, not a reason to throw them away.
                    print(f"[WARN] Coverage report for {name} not produced: "
                          f"{result['error']}")
                    continue
                summary = result["summary"]
                coverage_reports.append(result)
                print(f"[OK] {name}: "
                      f"{summary['regions_fully_covered']}/"
                      f"{summary['regions']} regions fully covered at "
                      f"{depth_floor}x -> {result['html']}")
                if not summary["all_covered"]:
                    print(f"[WARN] {name}: "
                          f"{summary['bases_below_threshold']:,} target "
                          f"bases below {depth_floor}x "
                          f"({summary['regions_absent']} regions with no "
                          f"usable coverage). A variant in those stretches "
                          f"would be absent from the VCF whether or not it "
                          f"is there.")

            if coverage_reports:
                manifest["steps_completed"].append("coverage")
                manifest["coverage_bed"] = os.path.abspath(args.coverage_bed)
                manifest["coverage_reports"] = [
                    {"sample": r["sample"], "html": r["html"],
                     "json": r["json"], "tsv": r["tsv"],
                     "summary": r["summary"]}
                    for r in coverage_reports]

    # ================================================================
    # STEP 13: PCGR CLINICAL REPORT
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 13/14: Clinical Interpretation Report (PCGR)")
    print(f"{'#' * 72}")

    if should_run("pcgr") and args.pcgr_refdata_dir:
        # pcgr_report.py lives beside this script. Import it lazily and
        # tolerate its absence, so this file keeps working standalone if
        # only it was copied somewhere.
        try:
            from pcgr_report import run_pcgr, lift_format_to_info
        except ImportError:
            print("[SKIP] pcgr_report.py not found next to this script; "
                  "skipping the PCGR report.")
            run_pcgr = None

        if run_pcgr is not None:
            # Deliberately NOT the SnpEff output: PCGR runs its own VEP
            # annotation and validates the INFO fields of what it is given,
            # so it gets the same call set SnpEff consumed (COSMIC-annotated
            # when COSMIC ran, otherwise the filtered VCF).
            pcgr_input = vcf_for_annotation
            pcgr_dp_tag = args.pcgr_tumor_dp_tag
            pcgr_af_tag = args.pcgr_tumor_af_tag
            pcgr_control_dp = None
            pcgr_control_af = None

            # Optional FORMAT -> INFO lift, so PCGR's depth/AF filters (and
            # therefore its TMB) actually have something to read. Writes a
            # separate file; the called VCF is never modified.
            if args.pcgr_lift_tags:
                lifted = os.path.join(
                    dirs["pcgr"], f"{tumour['name']}.pcgr_input.vcf")
                lift_stats = lift_format_to_info(
                    pcgr_input, lifted,
                    tumour_sample=tumour["name"],
                    normal_sample=normal["name"] if normal else None,
                    dry_run=args.dry_run,
                )
                if lift_stats is None:
                    # Refused (ambiguous samples, already lifted, unreadable).
                    # Not fatal: PCGR still reports, just without filters.
                    print("[WARN] FORMAT->INFO transform did not run; "
                          "continuing with the original VCF and no depth/AF "
                          "tags.")
                else:
                    pcgr_input = lifted
                    pcgr_dp_tag = pcgr_dp_tag or "TDP"
                    pcgr_af_tag = pcgr_af_tag or "TAF"
                    if lift_stats.get("normal_sample"):
                        pcgr_control_dp = "NDP"
                        pcgr_control_af = "NAF"

            # PCGR's genome vocabulary is lowercase and differs from
            # SnpEff's; derive it the same way the SnpEff genome is derived.
            pcgr_genome = "grch37" if snpeff_genome == "hg19" else "grch38"

            ok = run_pcgr(
                input_vcf=pcgr_input,
                output_dir=dirs["pcgr"],
                sample_id=tumour["name"],
                refdata_dir=args.pcgr_refdata_dir,
                vep_dir=args.vep_dir,
                genome_assembly=pcgr_genome,
                assay=args.pcgr_assay,
                effective_target_size_mb=args.pcgr_target_size_mb,
                estimate_tmb=args.pcgr_estimate_tmb,
                estimate_msi=args.pcgr_estimate_msi,
                estimate_signatures=args.pcgr_estimate_signatures,
                tumour_site=args.pcgr_tumour_site,
                # No matched normal means PCGR must apply its germline
                # filtering, or the report fills with inherited variants.
                tumour_only=normal is None,
                tumour_dp_tag=pcgr_dp_tag,
                tumour_af_tag=pcgr_af_tag,
                control_dp_tag=pcgr_control_dp,
                control_af_tag=pcgr_control_af,
                legacy_v1=args.pcgr_legacy_v1,
                extra_args=(shlex.split(args.pcgr_extra_args)
                            if args.pcgr_extra_args else None),
                log_path=os.path.join(dirs["logs"], f"{tumour['name']}.log"),
                dry_run=args.dry_run,
            )
            if ok is None:
                print("[SKIP] PCGR report skipped (pcgr or its reference "
                      "bundle unavailable).")
            elif ok:
                manifest["steps_completed"].append("pcgr")
                manifest["pcgr_dir"] = dirs["pcgr"]
                manifest["pcgr_input_vcf"] = pcgr_input
            else:
                print("[WARN] PCGR report failed.")
    else:
        # Two distinct reasons; say which.
        if not args.pcgr_refdata_dir:
            print("[SKIP] PCGR skipped (--pcgr-refdata-dir not provided).")
        else:
            print("[SKIP] PCGR skipped (--skip-steps pcgr).")

    # ================================================================
    # STEP 13: SUMMARY STATISTICS
    # ================================================================
    print(f"\n{'#' * 72}")
    print("# STEP 14/14: Summary Statistics")
    print(f"{'#' * 72}")

    if os.path.exists(filtered_vcf):
        manifest["vcf_stats"] = vcf_stats(filtered_vcf, dry_run=args.dry_run)

    # ================================================================
    # CLEAN UP LANE-MERGE TEMPORARIES
    # ================================================================
    # These are concatenated copies of the input FASTQs (auto-discover mode
    # only). Alignment has consumed them by now, and on real data they are
    # as large as the inputs themselves, so leaving them behind would
    # quietly double the run's disk footprint. Note this only runs on the
    # success path -- if the pipeline aborted early they are deliberately
    # left in place for post-mortem debugging.
    for stale in lane_merge_files:
        try:
            os.remove(stale)
        except OSError:
            pass  # Already gone, or not ours to delete -- not worth failing.

    # ================================================================
    # WRITE MANIFEST AND PRINT SUMMARY
    # ================================================================
    finished = datetime.now(timezone.utc)
    manifest["finished_utc"] = finished.isoformat()
    manifest["duration_seconds"] = round(
        (finished - started).total_seconds(), 1)
    manifest["result"] = {
        "tumour": tumour["name"],
        "normal": normal["name"] if normal else None,
        "steps_completed": manifest["steps_completed"],
    }

    manifest_path = os.path.join(
        output_dir,
        f"pipeline_manifest_{started.strftime('%Y%m%dT%H%M%SZ')}.json")
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    print(f"\n{'=' * 72}")
    print("PIPELINE COMPLETE")
    print(f"{'=' * 72}")
    print(f"Duration    : {manifest['duration_seconds']}s")
    print(f"Steps done  : {', '.join(manifest['steps_completed'])}")
    print(f"Tumour BAM  : {tumour_bqsr}")
    if normal:
        print(f"Normal BAM  : {normal_bqsr}")
    print(f"Filtered VCF: {filtered_vcf}")
    if os.path.exists(cosmic_vcf):
        print(f"COSMIC VCF  : {cosmic_vcf}")
    if os.path.exists(annotated_vcf):
        print(f"Annotated   : {annotated_vcf}")
    if "pcgr" in manifest["steps_completed"]:
        print(f"PCGR report : {dirs['pcgr']}")
    print(f"Metrics     : {dirs['metrics']}")
    print(f"Logs        : {dirs['logs']}")
    print(f"Manifest    : {manifest_path}")
    print(f"{'=' * 72}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
