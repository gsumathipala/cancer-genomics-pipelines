#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
align_rna.py
============
Stage 2 of the RNA branch: splice-aware alignment with STAR, configured so
that a fusion caller can use the result.

WHY NOT align_reads.py
----------------------
  bwa-mem2 cannot align RNA. A read that crosses an exon-exon junction has
  no contiguous match anywhere in the genome, so bwa-mem2 either soft-clips
  it back to one exon or throws it away. Both outcomes are silent: the BAM
  looks ordinary, the mapping rate looks only slightly low, and every
  junction-spanning read -- which is to say every read that could have
  evidenced a fusion -- is gone.

  STAR solves the same problem the transcript does: it allows an alignment
  to be split across a gap, and it uses the annotation (the GTF) to know
  where the plausible gaps are.

WHY THE CHIMERIC SETTINGS BELOW ARE NOT OPTIONAL
------------------------------------------------
  A fusion is a read whose two halves land on two different genes. STAR
  calls that a CHIMERIC alignment and, by default, does not look for them
  at all. Run STAR with its defaults and hand the BAM to Arriba and you get
  a clean, empty, entirely wrong fusion table.

  So this module sets the chimeric-detection parameters explicitly, as a
  block, following Arriba's documented invocation. They are written out
  rather than hidden behind a "--fusion-mode" flag precisely because the
  failure they prevent is invisible -- somebody reading the command has to
  be able to see that chimeric detection is on.

WHAT THIS MODULE DOES NOT DO
----------------------------
  * No duplicate marking. On RNA it is wrong, not merely unhelpful: two
    reads at one position in a highly expressed gene are two observations
    of an abundant transcript, not one molecule counted twice. There is no
    flag to enable it, because there is no correct value.
  * No BQSR, no variant calling. See RNA_SCOPE.md.
  * No expression quantification. Also scoped out, also deliberate.

THE INDEX
---------
  STAR's index is built from the genome FASTA *and* the GTF *and* a read
  length (--sjdbOverhang, conventionally read length - 1). All three are
  baked in. An index built for 100 bp reads works perfectly well on 150 bp
  reads and quietly loses junction sensitivity, so this module records what
  an index was built from and refuses to use one whose provenance does not
  match, rather than discovering the mismatch as a slightly disappointing
  fusion table months later.

Usage
-----
    # Build an index once (about an hour, ~32 GB RAM, ~30 GB disk):
    python align_rna.py --build-index \\
        --reference ~/data/references/hg38/hg38.fa \\
        --gtf ~/data/references/gencode/gencode.v44.annotation.gtf \\
        --read-length 151 --star-index ~/data/references/star_hg38_150

    # Align:
    python align_rna.py --manifest qc_out/run_manifest_*.json \\
        --star-index ~/data/references/star_hg38_150 \\
        -o rna_out/ --threads 16
"""

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

# Where a shared STAR index lives by default, beside the shared genome the
# DNA branch already uses. Sharing the FASTA is the point: an RNA run costs
# an index, not another 20 GB of reference.
DEFAULT_REFERENCE_DIR = os.environ.get(
    "PIPELINE_REFERENCE_DIR", os.path.expanduser("~/data/references"))

# Written into the index directory when this script builds one, and checked
# before this script uses one. STAR itself records some of this in
# genomeParameters.txt, but not the GTF's identity -- and the GTF is the
# half that decides which junctions exist.
INDEX_RECORD = "pipeline_index.json"

# Files STAR writes into an index directory. Their presence is how an
# existing index is recognised.
INDEX_SENTINELS = ("SA", "SAindex", "Genome", "genomeParameters.txt")


def run_command(cmd, tag, log_path=None, dry_run=False):
    """
    Run one command, echoing it first.

    Same contract as the DNA scripts': the command is printed before it
    runs so a log records what was actually executed, and --dry-run prints
    without executing.
    """
    printable = " ".join(shlex.quote(str(c)) for c in cmd)
    print(f"\n[{tag}] $ {printable}")
    if dry_run:
        print(f"[{tag}] (dry run - not executed)")
        return 0
    handle = open(log_path, "a", encoding="utf-8") if log_path else None
    try:
        if handle:
            handle.write(f"\n=== {tag} ===\n$ {printable}\n")
            handle.flush()
        proc = subprocess.run(cmd, stdout=handle or None,
                              stderr=subprocess.STDOUT if handle else None)
        return proc.returncode
    except FileNotFoundError:
        print(f"[{tag}] [ERROR] executable not found: {cmd[0]}")
        return 127
    finally:
        if handle:
            handle.close()


def tool_version(executable):
    """Version string for the run record, or None."""
    if shutil.which(executable) is None:
        return None
    for flag in ("--version", "-v", "--help"):
        try:
            out = subprocess.run([executable, flag], capture_output=True,
                                 text=True, timeout=30)
            text = (out.stdout or out.stderr).strip().splitlines()
            if text:
                return text[0][:120]
        except (OSError, subprocess.SubprocessError):
            continue
    return None


# =============================================================================
# SECTION 1: THE INDEX
# =============================================================================

def index_record_path(star_index):
    return os.path.join(star_index, INDEX_RECORD)


def read_index_record(star_index):
    """What an index was built from, or {} when it was not built by us."""
    try:
        with open(index_record_path(star_index), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def index_is_present(star_index):
    """True when a directory holds a complete-looking STAR index."""
    return all(os.path.exists(os.path.join(star_index, name))
               for name in INDEX_SENTINELS)


def check_index_provenance(star_index, gtf, read_length):
    """
    Compare an index against what this run needs, returning warnings.

    NOT fatal, and deliberately so. An index built elsewhere -- by a core
    facility, by a previous version of this script -- carries no record,
    and refusing to use it would be worse than using it. But an index whose
    record says it was built from a different GTF release, or for a
    materially different read length, is a real problem that is otherwise
    completely silent, so it is named.
    """
    warnings = []
    record = read_index_record(star_index)
    if not record:
        warnings.append(
            f"{os.path.basename(star_index)} carries no build record, so "
            f"the GTF and read length it was built from cannot be checked. "
            f"If it came from elsewhere, confirm it was built with the same "
            f"annotation you are about to hand the fusion caller: a "
            f"mismatch costs junction sensitivity and says nothing.")
        return warnings

    built_gtf = record.get("gtf")
    if gtf and built_gtf and os.path.basename(built_gtf) != \
            os.path.basename(gtf):
        warnings.append(
            f"this index was built from {os.path.basename(built_gtf)} but "
            f"the run is using {os.path.basename(gtf)}. The caller will "
            f"annotate breakpoints against one annotation and STAR found "
            f"junctions using another; the two disagreeing is exactly how a "
            f"real fusion ends up unannotated.")

    built_len = record.get("read_length")
    if read_length and built_len and abs(int(built_len) - int(read_length)) > 10:
        warnings.append(
            f"this index was built for {built_len} bp reads and this "
            f"library is {read_length} bp. STAR will run and lose junction "
            f"sensitivity in proportion to the difference -- rebuild the "
            f"index for this read length if fusion detection matters.")
    return warnings


def build_star_index(star_index, reference, gtf, read_length, threads,
                     log_path=None, dry_run=False, extra_args=None):
    """
    Build a STAR index from the shared genome and a GENCODE GTF.

    Expensive: roughly an hour on 16 threads, ~32 GB of RAM and ~30 GB of
    disk. It is built ONCE per (genome, annotation, read length) and shared
    by every run, which is the same bargain the bwa-mem2 index makes on the
    DNA side -- and the reason both live under --reference-dir rather than
    inside an output directory.

    RETURNS:
        bool: True on success.
    """
    os.makedirs(star_index, exist_ok=True)
    # STAR's own convention: the overhang is the length of the sequence
    # either side of a junction that a read could span, so read length - 1.
    overhang = max(int(read_length) - 1, 1)
    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--runThreadN", str(threads),
        "--genomeDir", star_index,
        "--genomeFastaFiles", reference,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", str(overhang),
        # STAR writes its scratch into the working directory otherwise,
        # which on a shared machine is somebody else's problem.
        "--outFileNamePrefix", os.path.join(star_index, "_build."),
    ]
    if extra_args:
        cmd += list(extra_args)

    print(f"[INFO] Building the STAR index in {star_index}")
    print(f"[INFO] Genome {os.path.basename(reference)}, annotation "
          f"{os.path.basename(gtf)}, sjdbOverhang {overhang} "
          f"(read length {read_length}).")
    print("[INFO] This takes about an hour and wants ~32 GB of RAM.")
    code = run_command(cmd, tag="STAR genomeGenerate", log_path=log_path,
                       dry_run=dry_run)
    if code != 0:
        print("[ERROR] STAR index build failed.")
        return False
    if dry_run:
        return True

    # The build record. Without it, the next run cannot tell this index
    # from one built for a different annotation -- see
    # check_index_provenance().
    record = {
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference": os.path.abspath(reference),
        "gtf": os.path.abspath(gtf),
        "read_length": int(read_length),
        "sjdb_overhang": overhang,
        "star_version": tool_version("STAR"),
    }
    try:
        with open(index_record_path(star_index), "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=2)
            fh.write("\n")
    except OSError as exc:
        print(f"[WARN] index built, but its build record could not be "
              f"written ({exc}). Later runs will not be able to check what "
              f"it was built from.")
    return True


# =============================================================================
# SECTION 2: ALIGNMENT
# =============================================================================

def chimeric_args():
    """
    STAR's chimeric-detection block, as Arriba documents it.

    Written as its own function so it appears verbatim in one place and can
    be read. Each line earns its keep:

      --chimSegmentMin 10        look for chimeric alignments at all, with
                                 each segment at least 10 bp. STAR's
                                 default is 0, meaning OFF -- this single
                                 value is the difference between a fusion
                                 panel and an empty table.
      --chimOutType WithinBAM
                       SoftClip  write chimeric alignments into the main
                                 BAM rather than a side file, which is what
                                 lets Arriba read one input.
      --chimJunctionOverhangMin  how far a read must extend past the
                                 breakpoint to count as evidence.
      --chimScoreDropMax 30      tolerate the alignment score penalty a
                                 genuine fusion incurs.
      --chimScoreJunctionNonGTAG 0
                                 a fusion breakpoint is a genomic
                                 rearrangement, not a splice site, so it
                                 must not be penalised for lacking the
                                 canonical GT/AG motif. Default is -1.
      --chimScoreSeparation 1
      --chimSegmentReadGapMax 3  tolerate short unaligned stretches at the
                                 junction.
      --chimMultimapNmax 50      fusions frequently involve repetitive
                                 partners; discarding multimappers discards
                                 them.
      --outFilterMultimapNmax 50 the same argument for the linear
                                 alignments.
      --peOverlapNbasesMin 10    use mate overlap, which on the short
                                 fragments FFPE RNA produces is most of the
                                 library.
      --alignSJstitchMismatchNmax 5 -1 5 5
                                 allow mismatches when stitching junctions,
                                 as Arriba's documented invocation does.
    """
    return [
        "--outSAMtype", "BAM", "Unsorted",
        "--outSAMunmapped", "Within",
        # Arriba reads the BAM once and sorts nothing, so compression is
        # pure cost here.
        "--outBAMcompression", "0",
        "--outFilterMultimapNmax", "50",
        "--peOverlapNbasesMin", "10",
        "--alignSplicedMateMapLminOverLmate", "0.5",
        "--alignSJstitchMismatchNmax", "5", "-1", "5", "5",
        "--chimSegmentMin", "10",
        "--chimOutType", "WithinBAM", "SoftClip",
        "--chimJunctionOverhangMin", "10",
        "--chimScoreDropMax", "30",
        "--chimScoreJunctionNonGTAG", "0",
        "--chimScoreSeparation", "1",
        "--chimSegmentReadGapMax", "3",
        "--chimMultimapNmax", "50",
    ]


def star_align(sample, r1, r2, star_index, out_dir, threads,
               log_path=None, dry_run=False, extra_args=None,
               read_files_command=None):
    """
    Align one sample with STAR, with chimeric detection on.

    ARGS:
        sample:   sample name; used for the output prefix.
        r1, r2:   FASTQs. r2 may be None for a single-end library (Ion
                  Torrent), which STAR handles and which every paired-end
                  assumption elsewhere in this bundle does not.
        out_dir:  where STAR's outputs go.
        extra_args: forwarded verbatim, appended last so they override.

    RETURNS:
        dict with the BAM path and STAR's log path, or None on failure.
    """
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, f"{sample}.")

    reads = [r1] + ([r2] if r2 else [])
    if read_files_command is None:
        # STAR does not detect compression. Handed a .gz without this it
        # reads the compressed bytes as sequence and reports ~0% mapping,
        # which looks like a failed library rather than a missing flag.
        read_files_command = "zcat" if str(r1).endswith(".gz") else "cat"

    cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", star_index,
        # NoSharedMemory: the webapp runs one job at a time and a shared
        # genome segment surviving a cancelled run is a 30 GB leak that
        # needs a manual STAR --genomeLoad Remove to clear.
        "--genomeLoad", "NoSharedMemory",
        "--readFilesIn"] + reads + [
        "--readFilesCommand", read_files_command,
        "--outFileNamePrefix", prefix,
    ]
    cmd += chimeric_args()
    if extra_args:
        cmd += list(extra_args)

    code = run_command(cmd, tag=f"STAR [{sample}]", log_path=log_path,
                       dry_run=dry_run)
    if code != 0:
        print(f"[ERROR] STAR failed for {sample}.")
        return None

    bam = prefix + "Aligned.out.bam"
    result = {
        "sample": sample,
        "bam": bam,
        "star_log": prefix + "Log.final.out",
        "sj_tab": prefix + "SJ.out.tab",
        "prefix": prefix,
    }
    if dry_run:
        return result
    if not os.path.exists(bam):
        print(f"[ERROR] STAR reported success but {bam} is not there.")
        return None
    return result


def sort_and_index(bam, output_bam, threads, log_path=None, dry_run=False):
    """
    Coordinate-sort and index an aligned BAM.

    Arriba wants the UNSORTED BAM -- it reads mates together, and sorting
    separates them -- so this is kept separate from star_align() and is for
    everything else: QC, IGV, and any later step that needs random access.
    """
    code = run_command(
        ["samtools", "sort", "-@", str(max(threads - 1, 1)),
         "-o", output_bam, bam],
        tag="samtools sort", log_path=log_path, dry_run=dry_run)
    if code != 0:
        print("[ERROR] samtools sort failed.")
        return False
    code = run_command(["samtools", "index", "-@", str(max(threads - 1, 1)),
                        output_bam],
                       tag="samtools index", log_path=log_path,
                       dry_run=dry_run)
    return code == 0


def parse_star_log(path):
    """
    STAR's Log.final.out as a dict of floats and ints.

    The file is a human-readable two-column report with '|' separators and
    percentages written with a trailing '%'. Parsed here rather than in the
    QC report so that both the aligner and the QC step read it the same
    way, and so a format change breaks in one place.

    Returns {} when the file cannot be read -- every caller treats these
    numbers as advisory.
    """
    stats = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if "|" not in line:
                    continue
                label, _sep, value = line.partition("|")
                label = label.strip().rstrip(":").strip()
                value = value.strip()
                if not value:
                    continue
                try:
                    if value.endswith("%"):
                        stats[label] = float(value[:-1])
                    elif "." in value:
                        stats[label] = float(value)
                    else:
                        stats[label] = int(value)
                except ValueError:
                    stats[label] = value
    except OSError:
        return {}
    return stats


# =============================================================================
# SECTION 3: CLI
# =============================================================================

def load_qc_manifest(manifest_path):
    """
    Cleaned FASTQ paths from a fastq_qc_clean.py manifest.

    IMPORTED, NOT REIMPLEMENTED. The DNA engine already reads this format
    and handles three things a fresh implementation gets wrong, because
    each was found the hard way:

      * the key is outputs.r1, not outputs.r1_clean;
      * a record with status "skipped_existing" carries NO "outputs" key at
        all, and its cleaned files nevertheless exist at the standard
        location -- dropping those samples silently loses every sample of a
        re-run QC directory;
      * the QC output may have been moved since it ran, so the reads are
        looked for beside the manifest as well.

    A second copy of that logic is a second copy that can drift, and the
    failure mode is a chain that silently analyses nothing. Fall back to a
    local read only if the DNA engine is not importable -- which is the
    same tolerance the rest of the bundle shows for a partial copy.
    """
    try:
        from comprehensive_variant_calling import (
            load_qc_manifest as _load)
    except ImportError:
        print("[WARN] comprehensive_variant_calling.py is not beside this "
              "script, so the QC manifest is read by a simplified reader "
              "that cannot recover samples whose QC was skipped because "
              "their output already existed.")
        return _load_qc_manifest_fallback(manifest_path)

    try:
        samples = _load(manifest_path)
    except (OSError, ValueError) as exc:
        print(f"[ERROR] cannot read manifest {manifest_path}: {exc}")
        return {}

    # The DNA engine always pairs r1 with r2 because its aligner requires
    # it. RNA does not: an Ion Torrent library is single-end and STAR
    # handles that, so a missing or absent r2 becomes None here rather than
    # a path that does not exist.
    for record in samples.values():
        r2 = record.get("r2")
        record["r2"] = r2 if r2 and os.path.isfile(r2) else None
    return samples


def _load_qc_manifest_fallback(manifest_path):
    """Minimal manifest reader, for a bundle missing the DNA engine."""
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError) as exc:
        print(f"[ERROR] cannot read manifest {manifest_path}: {exc}")
        return {}

    samples = {}
    for record in data.get("samples", []):
        if record.get("status") == "failed":
            continue
        name = record.get("sample")
        outputs = record.get("outputs") or {}
        r1, r2 = outputs.get("r1"), outputs.get("r2")
        if name and r1 and os.path.isfile(r1):
            samples[name] = {"r1": r1,
                             "r2": r2 if r2 and os.path.isfile(r2) else None}
    return samples


def build_parser():
    parser = argparse.ArgumentParser(
        description="Splice-aware RNA alignment (STAR) configured for "
                    "fusion detection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="See RNA.md for the whole RNA branch, and RNA_SCOPE.md for "
               "what it deliberately does not do.")

    io = parser.add_argument_group("input / output")
    io.add_argument("-o", "--output-dir", default=None,
                    help="Directory for the aligned BAMs and STAR logs. "
                         "Not needed with --build-index.")
    io.add_argument("--manifest", default=None,
                    help="fastq_qc_clean.py run manifest; cleaned FASTQ "
                         "paths and sample names are read from it.")
    io.add_argument("--sample", default=None,
                    help="Sample name, for the explicit input mode.")
    io.add_argument("--r1", default=None, help="Read 1 FASTQ.")
    io.add_argument("--r2", default=None,
                    help="Read 2 FASTQ. Omit for a single-end library.")

    ref = parser.add_argument_group("reference / index")
    ref.add_argument("-r", "--reference", default=None,
                     help="Genome FASTA. Needed only to BUILD an index; "
                          "alignment reads the genome out of the index.")
    ref.add_argument("--gtf", default=None,
                     help="GENCODE annotation GTF. Required to build an "
                          "index, and recorded so a later run can check "
                          "that the annotation has not changed under it.")
    ref.add_argument("--star-index", default=None,
                     help="STAR index directory. Shared between runs like "
                          "the bwa-mem2 index is; keep it under "
                          "--reference-dir rather than in an output "
                          "directory.")
    ref.add_argument("--reference-dir", default=DEFAULT_REFERENCE_DIR,
                     help="Where shared references live, used to guess "
                          "--star-index when it is not given.")
    ref.add_argument("--read-length", type=int, default=None,
                     help="Library read length. Sets STAR's --sjdbOverhang "
                          "(read length - 1) when building, and is checked "
                          "against the index's build record when aligning. "
                          "An index built for the wrong read length works "
                          "and quietly loses junction sensitivity.")
    ref.add_argument("--build-index", action="store_true",
                     help="Build the index and exit. About an hour, ~32 GB "
                          "RAM, ~30 GB disk, once per genome+annotation+"
                          "read length.")

    run = parser.add_argument_group("runtime")
    run.add_argument("--threads", type=int, default=8,
                     help="CPU threads.")
    run.add_argument("--sorted-bam", action="store_true",
                     help="Also write a coordinate-sorted, indexed copy. "
                          "The fusion caller wants the unsorted BAM (it "
                          "reads mates together), so this is for QC, IGV "
                          "and anything needing random access.")
    run.add_argument("--star-extra-args", default=None,
                     help="Extra arguments forwarded verbatim to STAR, as "
                          "one quoted string, appended last so they "
                          "override what this script builds.")
    run.add_argument("--dry-run", action="store_true",
                     help="Print the commands without running them.")

    panel = parser.add_argument_group("panel / assay profile")
    panel.add_argument("--panel", default=None, metavar="ID",
                       help="RNA panel profile to take the read length, "
                            "annotation and STAR settings from.")
    panel.add_argument("--panel-file", action="append", default=[],
                       metavar="JSON", help="Additional profile file.")
    panel.add_argument("--list-panels", action="store_true",
                       help="List every known panel profile and exit.")
    panel.add_argument("--describe-panel", default=None, metavar="ID",
                       help="Print everything a profile sets, and exit.")
    return parser


def apply_panel(args, parser, argv):
    """Take the RNA settings this script implements from --panel."""
    if not args.panel:
        return
    try:
        import panel_profiles
    except ImportError as exc:
        parser.error(f"--panel needs panel_profiles.py beside this script "
                     f"({exc}).")
        return
    try:
        registry = panel_profiles.available_panels(
            extra_files=args.panel_file,
            warn=lambda msg: print(f"[WARN] {msg}"))
        profile = panel_profiles.resolve_panel(args.panel, registry)
        if panel_profiles.assay_type(profile) != "rna":
            parser.error(
                f"'{profile['id']}' is a DNA profile ({profile['chemistry']}) "
                f"and this is the RNA aligner. RNA profiles are listed by "
                f"--list-panels with [RNA] beside them.")
        applied, overridden = panel_profiles.apply_profile(
            args, profile, panel_profiles.explicit_dests(parser, argv),
            only=panel_profiles.RNA_SETTINGS)
    except panel_profiles.PanelError as exc:
        parser.error(str(exc))
        return
    print()
    print(panel_profiles.format_application(profile, applied, overridden))


def installed_star_indexes(reference_dir):
    """
    Every complete STAR index under the shared reference directory.

    Returns [(path, read_length_or_None)], sorted by read length.
    """
    found = []
    try:
        entries = sorted(os.listdir(os.path.expanduser(reference_dir)))
    except OSError:
        return found
    for name in entries:
        if not name.startswith("star_"):
            continue
        path = os.path.join(os.path.expanduser(reference_dir), name)
        if not os.path.exists(os.path.join(path, "SAindex")):
            continue
        tail = name.rsplit("_", 1)[-1]
        found.append((path, int(tail) if tail.isdigit() else None))
    return sorted(found, key=lambda pair: (pair[1] is None, pair[1] or 0))


def default_star_index(reference_dir, read_length):
    """
    Which index to use when --star-index was not given.

    Named by read length, because that is what makes two indexes different
    and otherwise indistinguishable on disk -- so with a read length this
    is a straight construction.

    WITHOUT one it DISCOVERS rather than guesses. Constructing
    "star_hg38" was the obvious thing and the wrong thing: the installer
    never creates that name (it always appends the read length), so a run
    that omitted --read-length was told there was no index while a
    perfectly good star_hg38_150 sat beside the path it had invented.
    """
    reference_dir = os.path.expanduser(reference_dir)
    if read_length:
        return os.path.join(reference_dir, f"star_hg38_{int(read_length)}")

    found = installed_star_indexes(reference_dir)
    if len(found) == 1:
        path, length = found[0]
        print(f"[INFO] --read-length was not given; using the one installed "
              f"STAR index, {os.path.basename(path)}"
              + (f" (built for {length} bp reads)." if length else "."))
        return path
    if len(found) > 1:
        names = ", ".join(os.path.basename(p) for p, _ in found)
        print(f"[WARN] several STAR indexes are installed ({names}) and "
              f"--read-length was not given, so the right one cannot be "
              f"chosen. Using {os.path.basename(found[0][0])}; pass "
              f"--read-length or --star-index to be explicit. An index "
              f"built for the wrong read length works and quietly loses "
              f"junction sensitivity.")
        return found[0][0]
    return os.path.join(reference_dir, "star_hg38")


def main():
    try:
        import panel_profiles
        panel_profiles.handle_panel_queries()
    except ImportError:
        pass

    parser = build_parser()
    args = parser.parse_args()
    apply_panel(args, parser, sys.argv[1:])

    if not args.star_index:
        args.star_index = default_star_index(args.reference_dir,
                                             args.read_length)
        print(f"[INFO] --star-index not given; using {args.star_index}")

    if shutil.which("STAR") is None and not args.dry_run:
        print("[ERROR] STAR not found on PATH. It lives in the cancer_rna "
              "environment: conda activate cancer_rna")
        return 1

    # ---- build mode ----
    if args.build_index:
        missing = [flag for flag, value in
                   (("--reference", args.reference), ("--gtf", args.gtf),
                    ("--read-length", args.read_length)) if not value]
        if missing:
            parser.error(f"--build-index needs {', '.join(missing)}.")
        if index_is_present(args.star_index):
            print(f"[SKIP] a STAR index is already present in "
                  f"{args.star_index}. Delete it to rebuild.")
            return 0
        ok = build_star_index(
            args.star_index, args.reference, args.gtf, args.read_length,
            args.threads, dry_run=args.dry_run,
            extra_args=(shlex.split(args.star_extra_args)
                        if args.star_extra_args else None))
        return 0 if ok else 1

    # ---- align mode ----
    if not args.output_dir:
        parser.error("-o/--output-dir is required when aligning.")
    if not index_is_present(args.star_index) and not args.dry_run:
        print(f"[ERROR] no STAR index in {args.star_index}. Build one:")
        print(f"        python align_rna.py --build-index --reference "
              f"<hg38.fa> \\\n"
              f"            --gtf <gencode.gtf> --read-length "
              f"{args.read_length or '<len>'} \\\n"
              f"            --star-index {args.star_index}")
        return 1

    for warning in check_index_provenance(args.star_index, args.gtf,
                                          args.read_length):
        print(f"[WARN] {warning}")

    samples = {}
    if args.manifest:
        samples = load_qc_manifest(args.manifest)
        if not samples:
            print("[ERROR] no usable samples in the manifest.")
            return 1
    elif args.sample and args.r1:
        samples = {args.sample: {"r1": args.r1, "r2": args.r2}}
    else:
        parser.error("Provide --manifest, or --sample with --r1 "
                     "(and --r2 for a paired library).")

    output_dir = os.path.abspath(args.output_dir)
    align_dir = os.path.join(output_dir, "rna_aligned")
    log_dir = os.path.join(output_dir, "logs")
    for directory in (align_dir, log_dir):
        os.makedirs(directory, exist_ok=True)

    results = []
    for name in sorted(samples):
        reads = samples[name]
        print(f"\n{'=' * 72}\nSTAR: {name}\n{'=' * 72}")
        result = star_align(
            name, reads["r1"], reads.get("r2"), args.star_index, align_dir,
            args.threads, log_path=os.path.join(log_dir, f"{name}.log"),
            dry_run=args.dry_run,
            extra_args=(shlex.split(args.star_extra_args)
                        if args.star_extra_args else None))
        if result is None:
            return 1
        if args.sorted_bam:
            sorted_bam = os.path.join(align_dir, f"{name}.sorted.bam")
            if sort_and_index(result["bam"], sorted_bam, args.threads,
                              log_path=os.path.join(log_dir, f"{name}.log"),
                              dry_run=args.dry_run):
                result["sorted_bam"] = sorted_bam
        stats = parse_star_log(result["star_log"]) if not args.dry_run else {}
        result["star_stats"] = stats
        if stats:
            print(f"[OK] {name}: "
                  f"{stats.get('Uniquely mapped reads %', '?')}% uniquely "
                  f"mapped, "
                  f"{stats.get('Number of chimeric reads', '?')} chimeric "
                  f"reads.")
        results.append(result)

    manifest = {
        "script": os.path.basename(__file__),
        "finished_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "command_line": shlex.join(sys.argv),
        "star_index": os.path.abspath(args.star_index),
        "index_record": read_index_record(args.star_index),
        "star_version": tool_version("STAR"),
        "samples": results,
    }
    path = os.path.join(output_dir, "rna_align_manifest.json")
    if not args.dry_run:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"\n[OK] Alignment manifest: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
