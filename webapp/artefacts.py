#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
artefacts.py
============
Every file a run produced, described in plain English and ordered by the
question it answers.

WHO THIS IS FOR
---------------
  A molecular pathologist who has just been handed this interface, knows
  the biology cold, and has never seen a VCF directory. Before this module
  the outputs of a run were real, correct, and effectively invisible: a
  PCGR link, a coverage link, and thirty other files sitting in an output
  directory nobody could reach without a terminal.

  The failure that produces is not "an inconvenience". It is a report read
  without its coverage statement, because the coverage statement was three
  directory levels away and the PCGR page was one click away.

THE ORDERING IS THE ARGUMENT
----------------------------
  Groups are deliberately NOT ordered by file type or by pipeline stage.
  They are ordered the way the result should be read:

      1. Read this first      the clinical interpretation
      2. Can I trust it?      coverage / library adequacy / read QC
      3. The underlying data  VCFs, tables
      4. Technical record     manifests, logs, alignment statistics

  Group 2 comes before group 3 for the same reason target coverage
  precedes the variant table in the PDF, and library adequacy precedes the
  fusion table in the RNA one: a reader who sees the findings first has
  already drawn a conclusion by the time they reach the caveat.

  Within each group, what a novice most needs comes first.

EVERY ENTRY SAYS WHAT IT ANSWERS
--------------------------------
  Not what it is. "coverage.html" is a filename; "which parts of the panel
  were actually sequenced deeply enough to call a variant" is the question
  somebody has. The second is what gets shown.

MISSING FILES ARE REPORTED, NOT OMITTED
---------------------------------------
  A run without a coverage report must SAY it has no coverage report and
  why, because a missing section reads as "nothing to report" rather than
  "never measured". absent_notes() supplies those lines.

SECURITY
--------
  Discovery globs inside the run's own two directories and nowhere else,
  and the serving route addresses a file by its INDEX in this list -- never
  by a path from the request. resolve() re-checks containment before
  anything is served, so a job record edited on disk cannot turn into a
  file-read primitive.
"""

import glob
import hashlib
import os

# How a file is offered. HTML and PDF open in a tab; everything else
# downloads, because a browser showing 40 MB of VCF inline helps nobody.
VIEW = "view"
DOWNLOAD = "download"

# Group ids, in the order they are shown. See the module docstring: this
# sequence is a clinical argument, not a filing convention.
GROUPS = [
    ("report", "Read this first",
     "The interpretation. Everything else on this page supports or "
     "qualifies it."),
    ("trust", "Can this result be trusted?",
     "Whether the sequencing was good enough for the findings above to "
     "mean what they appear to mean. Read these BEFORE reporting a "
     "negative result."),
    ("data", "The underlying data",
     "The variant and fusion tables themselves, for review or for loading "
     "into another tool."),
    ("record", "Technical record",
     "What was run, with which settings and which reference data. Needed "
     "to reproduce or audit a result, rarely needed to read one."),
]


def artefact_id(path, run_dir, output_dir):
    """
    A stable id for one file: derived from WHERE it is, never from
    where it happens to sit in a list.

    THIS IS NOT A MICRO-OPTIMISATION. These reports are meant to be read
    while the run is still going -- that is the point of the coverage and
    library-QC steps. So the set of files GROWS under the page: a batch run
    writing SAMPLE_A's coverage report after the page was rendered shifted
    every position after it, and a link labelled "Read quality" then served
    the previous sample's coverage report instead. A 200, the wrong file,
    and nothing anywhere to say so.

    The id is a digest of the path relative to the run, so it survives new
    files appearing, files being removed, and the run directory being
    moved.
    """
    for root in (output_dir, run_dir):
        if root:
            try:
                rel = os.path.relpath(path, root)
            except ValueError:            # different drives, on Windows
                continue
            if not rel.startswith(".."):
                return hashlib.sha1(rel.encode("utf-8")).hexdigest()[:12]
    return hashlib.sha1(os.path.abspath(path).encode("utf-8")).hexdigest()[:12]


def _entry(path, label, what, group, action=VIEW, caution=None):
    """One artefact. `what` is the question it answers, in plain English."""
    return {"path": path, "label": label, "what": what, "group": group,
            "action": action, "caution": caution,
            "size": os.path.getsize(path) if os.path.exists(path) else 0}


def _stem(path, suffix):
    """Sample name out of a report filename."""
    name = os.path.basename(path)
    return name[:-len(suffix)] if name.endswith(suffix) else name


def _sorted_glob(*parts):
    return sorted(glob.glob(os.path.join(*parts)))


# =============================================================================
# DISCOVERY
# =============================================================================

def collect(job_snapshot, run_dir, output_dir):
    """
    Every artefact of one run, ordered for reading.

    Takes the snapshot rather than the Job so this module stays importable
    and testable without the job machinery.
    """
    assay = (job_snapshot or {}).get("assay", "dna")
    items = []
    run_dir = run_dir or ""
    output_dir = output_dir or ""

    # ---------------------------------------------------------------
    # 1. READ THIS FIRST
    # ---------------------------------------------------------------
    for pdf in _sorted_glob(run_dir, "report_*.pdf"):
        items.append(_entry(
            pdf, "Signed-out PDF report",
            "The whole run as one document, with the patient details, what "
            "was found, and the limitations that apply to it. This is the "
            "one to file or to hand on.",
            "report", DOWNLOAD))

    # PCGR's own report, on BOTH branches. PCGR 2.x accepts RNA fusions as
    # a molecular input in their own right, so an RNA run has a real
    # clinical interpretation too -- and when the run was linked to a DNA
    # run, one report covers the specimen's SNVs and its fusions together.
    for path in _sorted_glob(output_dir, "pcgr", "*.html"):
        items.append(_entry(
            path, f"PCGR clinical report — {os.path.basename(path)}",
            "The authoritative interpretation: findings placed into "
            "actionability tiers against the tumour site. On a DNA run that "
            "is variants, TMB and mutational signatures; on an RNA run it "
            "is the fusions; where a run was linked to its partner, it is "
            "both in one document.",
            "report", VIEW))

    if assay == "rna":
        for path in _sorted_glob(output_dir, "fusion_report",
                                 "*.fusions.html"):
            sample = _stem(path, ".fusions.html")
            items.append(_entry(
                path, f"Fusion report — {sample}",
                "The fusions found, ordered by clinical importance rather "
                "than by the caller's confidence, with the actionable genes "
                "marked and same-gene splice events (such as MET exon 14 "
                "skipping) called out separately.",
                "report", VIEW))

    # ---------------------------------------------------------------
    # 2. CAN THIS RESULT BE TRUSTED?
    # ---------------------------------------------------------------
    if assay == "rna":
        for path in _sorted_glob(output_dir, "rna_qc", "*.rna_qc.html"):
            sample = _stem(path, ".rna_qc.html")
            items.append(_entry(
                path, f"Library adequacy — {sample}",
                "Whether an empty or sparse fusion result from this sample "
                "may be reported as a negative at all. A degraded RNA "
                "extraction produces a clean, empty fusion table that looks "
                "exactly like a true negative.",
                "trust", VIEW,
                caution="Read this before reporting 'no fusion detected'."))
    else:
        for path in _sorted_glob(output_dir, "coverage", "*.coverage.html"):
            sample = _stem(path, ".coverage.html")
            items.append(_entry(
                path, f"Target coverage — {sample}",
                "Which parts of the panel were sequenced deeply enough for "
                "a variant to have been called there. A region nobody "
                "sequenced and a region that is wild type look identical in "
                "the variant table; only this tells them apart.",
                "trust", VIEW,
                caution="Read this before reporting any negative finding."))

    for path in _sorted_glob(output_dir, "fastp_reports", "*.fastp.html"):
        sample = _stem(path, ".fastp.html")
        items.append(_entry(
            path, f"Read quality — {sample}",
            "The sequencing reads before and after trimming: quality "
            "scores, adapter content, duplication and length. A failed "
            "library usually shows here first.",
            "trust", VIEW))

    for path in _sorted_glob(output_dir, "fastqc_raw", "*.html") + \
            _sorted_glob(output_dir, "fastqc_cleaned", "*.html"):
        items.append(_entry(
            path, f"FastQC — {os.path.basename(path)[:-5]}",
            "The standard per-sample read-quality report, if FastQC was "
            "run.", "trust", VIEW))

    if assay == "rna":
        for path in _sorted_glob(output_dir, "rna_aligned", "*Log.final.out"):
            sample = os.path.basename(path).split(".")[0]
            items.append(_entry(
                path, f"Alignment statistics — {sample}",
                "STAR's own summary: how many reads aligned uniquely, how "
                "many were lost as too short (the signature of a degraded "
                "library), and how many crossed a chimeric junction.",
                "trust", DOWNLOAD))
    else:
        for path in _sorted_glob(output_dir, "stats", "*flagstat*") + \
                _sorted_glob(output_dir, "stats", "*idxstats*"):
            items.append(_entry(
                path, f"Alignment statistics — {os.path.basename(path)}",
                "How many reads aligned, and where. A low mapping rate "
                "undermines everything downstream.",
                "trust", DOWNLOAD))
        for path in _sorted_glob(output_dir, "metrics", "*dup_metrics.txt"):
            items.append(_entry(
                path, f"Duplicate metrics — {os.path.basename(path)}",
                "What fraction of reads were PCR or optical duplicates. "
                "Very high duplication means the library was sequenced "
                "beyond its real complexity.",
                "trust", DOWNLOAD))

    # ---------------------------------------------------------------
    # 3. THE UNDERLYING DATA
    # ---------------------------------------------------------------
    if assay == "rna":
        for path in _sorted_glob(output_dir, "fusion_report",
                                 "*.fusions.annotated.tsv"):
            items.append(_entry(
                path, f"Fusion table — {_stem(path, '.fusions.annotated.tsv')}",
                "The same fusions as the report above, as a spreadsheet: "
                "gene pair, confidence, supporting reads, breakpoints.",
                "data", DOWNLOAD))
        for path in _sorted_glob(output_dir, "fusion_report",
                                 "*.pcgr_fusions.tsv"):
            items.append(_entry(
                path, f"PCGR fusion input — {_stem(path, '.pcgr_fusions.tsv')}",
                "The same fusions in the format PCGR reads. Hand this to "
                "pcgr_report.py --input-rna-fusion, with the DNA run's VCF "
                "alongside it, to build one combined report for the "
                "specimen.",
                "data", DOWNLOAD))
        for path in _sorted_glob(output_dir, "fusions", "*.fusions.tsv"):
            items.append(_entry(
                path, f"Caller output (raw) — {_stem(path, '.fusions.tsv')}",
                "Arriba's untouched output, before this pipeline reordered "
                "or annotated anything.",
                "data", DOWNLOAD))
        for path in _sorted_glob(output_dir, "fusions",
                                 "*.fusions.discarded.tsv"):
            items.append(_entry(
                path, f"Discarded candidates — {_stem(path, '.fusions.discarded.tsv')}",
                "Fusion candidates the caller rejected, with its reason. "
                "Worth checking only when a clinically expected fusion is "
                "absent from the report.",
                "data", DOWNLOAD,
                caution="Rejected calls. Not findings."))
    else:
        for path in _sorted_glob(output_dir, "pcgr", "*.tsv") + \
                _sorted_glob(output_dir, "pcgr", "*.tsv.gz"):
            if path.endswith(".pcgr_input.vcf"):
                continue
            items.append(_entry(
                path, f"PCGR variant table — {os.path.basename(path)}",
                "Every variant PCGR considered, with its tier, consequence "
                "and supporting evidence, as a spreadsheet.",
                "data", DOWNLOAD))
        for pattern, label, what in (
            ("annotated/*.vcf.gz", "Annotated variants",
             "The final variant calls with gene and consequence "
             "annotation. This is the file to load into a genome browser "
             "or hand to another tool."),
            ("mutect2/*.mutect2.vcf.gz", "Raw caller output",
             "Everything Mutect2 proposed, before filtering. Most of it is "
             "not real; consult it only when a expected variant is missing "
             "from the annotated set."),
        ):
            for path in _sorted_glob(output_dir, *pattern.split("/")):
                items.append(_entry(
                    path, f"{label} — {os.path.basename(path)}", what,
                    "data", DOWNLOAD,
                    caution=("Unfiltered. Not a result set."
                             if "Raw" in label else None)))

    # ---------------------------------------------------------------
    # 4. TECHNICAL RECORD
    # ---------------------------------------------------------------
    for pattern in ("pipeline_manifest_*.json", "rna_manifest_*.json"):
        for path in _sorted_glob(output_dir, pattern):
            items.append(_entry(
                path, "Run manifest",
                "Exactly what was run: every setting, every reference file, "
                "every tool version, and which steps completed. This is the "
                "record that makes a result reproducible.",
                "record", DOWNLOAD))

    log = os.path.join(run_dir, "webapp_run.log")
    if os.path.exists(log):
        items.append(_entry(
            log, "Full run log",
            "Everything the pipeline printed, in order. The place to look "
            "when a step reports that it was skipped and you want to know "
            "why.",
            "record", DOWNLOAD))

    for path in _sorted_glob(output_dir, "logs", "*.log"):
        items.append(_entry(
            path, f"Tool log — {os.path.basename(path)}",
            "The raw output of the external tools for one sample.",
            "record", DOWNLOAD))

    # Ids derived from each file's own location, NOT from its position in
    # this list -- see artefact_id(). Display order is still the order
    # above; only the identity is independent of it.
    for item in items:
        item["id"] = artefact_id(item["path"], run_dir, output_dir)
    return items


def grouped(items):
    """[(group_id, title, blurb, [items])] in display order, empties dropped."""
    out = []
    for group_id, title, blurb in GROUPS:
        members = [i for i in items if i["group"] == group_id]
        if members:
            out.append((group_id, title, blurb, members))
    return out


def absent_notes(job_snapshot, items, meta):
    """
    What this run did NOT produce, and what that means.

    A missing report has to be stated. Left out, its absence reads as
    "nothing to report" rather than "this was never measured" -- which is
    the whole failure mode the coverage and library-QC steps exist to
    prevent, reappearing one level up in the interface.
    """
    assay = (job_snapshot or {}).get("assay", "dna")
    have = {i["group"] for i in items}
    labels = " ".join(i["label"] for i in items)
    notes = []

    if assay == "rna":
        if "Library adequacy" not in labels:
            notes.append(
                "No library adequacy report. Nothing here says whether an "
                "empty or sparse fusion result reflects the specimen or an "
                "unusable library, so it cannot be reported as a negative.")
        if "PCGR clinical report" not in labels:
            notes.append(
                "No PCGR interpretation for this run. Either no fusion "
                "cleared the confidence bar (PCGR rejects an empty fusion "
                "file, so 'nothing found' produces no report), or no "
                "reference bundle was configured. The fusion report below "
                "is this pipeline's own ranking, not a tiered "
                "interpretation.")
        if "Fusion report" not in labels:
            notes.append(
                "No fusion report. Either fusion calling did not complete, "
                "or the reporting step was skipped -- the run log says "
                "which.")
    else:
        if "Target coverage" not in labels:
            notes.append(
                "No coverage report, because this run was submitted without "
                "a target BED. Nothing distinguishes a region that was "
                "sequenced and is wild type from one the sequencing never "
                "reached: both are simply absent from the variant table.")
        if "PCGR clinical report" not in labels:
            notes.append(
                "No PCGR report, so there are no actionability tiers, no "
                "TMB and no mutational signatures. What exists is an "
                "annotated variant table, not an interpretation.")

    if "Read quality" not in labels:
        notes.append(
            "No read-quality report: QC was skipped, or the reads were "
            "cleaned by an earlier run.")
    if "report" not in have:
        notes.append(
            "No PDF has been generated yet. Use the button above; it can be "
            "regenerated at any time.")
    return notes


def resolve(items, artefact, run_dir, output_dir):
    """
    The path for one artefact id, re-checked for containment.

    Defence in depth. The id is a digest and cannot express a path, but a
    job record is a file on disk and a file on disk can be edited. This
    confirms the resolved path really does sit inside the run's own two
    directories before anything is served.

    Returns (path, item) or (None, None).
    """
    match = next((i for i in items if i["id"] == artefact), None)
    if match is None:
        return None, None
    path = os.path.realpath(match["path"])
    for root in (run_dir, output_dir):
        if not root:
            continue
        root = os.path.realpath(root)
        if path == root or path.startswith(root + os.sep):
            return path, match
    return None, None
