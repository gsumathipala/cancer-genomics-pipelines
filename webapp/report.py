#!/usr/bin/env python3
"""
report.py
=========
Assemble a PDF report from a finished pipeline run plus the patient details
captured in the web form.

WHAT THIS REPORT IS, AND IS NOT
-------------------------------
  It is a run summary: who the specimen belongs to, exactly how it was
  processed, what the pipeline observed, and where the full outputs live.
  Every number in it is read back from the artefacts the pipeline itself
  wrote (pipeline_manifest_*.json, the PCGR TSV) rather than recomputed
  here, so the PDF cannot disagree with the run.

  It is NOT a diagnostic report and does not interpret variants. PCGR's own
  HTML report is the interpretation, and this PDF points at it. The
  disclaimer at the end of the document says so explicitly, because a PDF
  with a patient name at the top will be treated as a clinical document by
  whoever finds it later.

DEFENSIVE READING
-----------------
  Everything parsed here comes from files that may be absent, truncated by
  a crash, or written by a different version of PCGR than the one this was
  written against. Every reader below therefore degrades to "not available"
  instead of raising -- a report that is missing a section is useful, one
  that fails to generate is not.
"""

import csv
import datetime
import glob
import json
import os

from pdfwriter import PDFDocument, unsupported_characters

# Fields captured by the patient form, in the order they are printed.
PATIENT_FIELDS = [
    ("patient_id", "Patient / MRN"),
    ("name", "Name"),
    ("dob", "Date of birth"),
    ("sex", "Sex"),
    ("diagnosis", "Clinical diagnosis"),
    ("specimen_id", "Specimen ID"),
    ("specimen_type", "Specimen type"),
    ("collection_date", "Collection date"),
    ("referring_clinician", "Referring clinician"),
    ("notes", "Notes"),
]


# ---------------------------------------------------------------------------
# Readers for the pipeline's own artefacts
# ---------------------------------------------------------------------------
def load_pipeline_manifest(output_dir):
    """Return the newest pipeline_manifest_*.json as a dict, or {}."""
    hits = sorted(glob.glob(os.path.join(output_dir,
                                         "pipeline_manifest_*.json")))
    if not hits:
        return {}
    try:
        with open(hits[-1], encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}


def find_pcgr_outputs(output_dir):
    """
    Locate PCGR's outputs under <output_dir>/pcgr.

    Filenames have shifted between PCGR releases, so glob rather than
    reconstruct exact names.
    """
    pcgr_dir = os.path.join(output_dir, "pcgr")
    found = {"dir": pcgr_dir, "html": [], "tsv": []}
    if not os.path.isdir(pcgr_dir):
        return found
    found["html"] = sorted(glob.glob(os.path.join(pcgr_dir, "*.html")))
    found["tsv"] = sorted(
        p for p in glob.glob(os.path.join(pcgr_dir, "*.tsv*"))
        # The lifted input VCF is written into this directory too; never
        # mistake our own input for a PCGR result.
        if not p.endswith(".pcgr_input.vcf")
    )
    return found


def read_pcgr_variants(tsv_paths, limit=40):
    """
    Pull a variant table out of a PCGR TSV, if one is recognisable.

    PCGR's column names vary by release, so match on a set of plausible
    aliases per column and skip anything unrecognised. Returns
    (headers, rows, source_path) or (None, [], None).

    NOTE: this has never been run against real PCGR output -- see the
    project's engineering notes. It is written to fail soft: if the columns
    do not match, the PDF simply omits the variant table and still points
    at PCGR's own HTML report.
    """
    aliases = {
        "Gene": ("SYMBOL", "GENE_SYMBOL", "SYMBOL_ENTREZ", "GENE"),
        "Variant": ("HGVSp_short", "HGVSP_SHORT", "PROTEIN_CHANGE",
                    "HGVSc", "CONSEQUENCE"),
        "Consequence": ("CONSEQUENCE", "VARIANT_CLASS", "MUTATION_TYPE"),
        "VAF": ("TVAF", "AF_TUMOR", "TAF", "VAF"),
        "Tier": ("TIER", "ACTIONABILITY_TIER", "CLINICAL_SIGNIFICANCE"),
    }
    for path in tsv_paths:
        try:
            with open(path, newline="", encoding="utf-8",
                      errors="replace") as fh:
                # PCGR TSVs sometimes carry leading '#' comment lines.
                lines = [l for l in fh if not l.startswith("##")]
            if not lines:
                continue
            reader = csv.DictReader(lines, delimiter="\t")
            cols = reader.fieldnames or []
            chosen = {}
            for label, options in aliases.items():
                for opt in options:
                    if opt in cols:
                        chosen[label] = opt
                        break
            # Require at least a gene column, else this is not a variant TSV.
            if "Gene" not in chosen:
                continue
            headers = list(chosen)
            rows = []
            for rec in reader:
                rows.append([(rec.get(chosen[h]) or "").strip()
                             for h in headers])
                if len(rows) >= limit:
                    break
            if rows:
                return headers, rows, path
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
    return None, [], None


def summarise_run(manifest):
    """Condense the pipeline manifest into printable key/value rows."""
    stats = manifest.get("vcf_stats") or {}
    overlap = manifest.get("cosmic_overlap") or {}

    def stat(*names):
        for n in names:
            if n in stats:
                return stats[n]
        return None

    rows = [
        ("Variants (records)", stat("number of records", "number of SNPs")),
        ("SNVs", stat("number of SNPs")),
        ("Indels", stat("number of indels")),
        ("Multiallelic sites", stat("number of multiallelic sites")),
        ("Ti/Tv ratio", stat("ts/tv")),
    ]
    if overlap:
        rows += [
            ("Variants in COSMIC", overlap.get("cosmic_overlapping")),
            ("COSMIC overlap", (f"{overlap['cosmic_overlap_pct']}%"
                                if overlap.get("cosmic_overlap_pct") is not None
                                else None)),
        ]
    return [(k, v) for k, v in rows if v not in (None, "")]


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------
def patient_warnings(patient):
    """
    Characters in the patient details that the PDF cannot render.

    Surfaced in the UI so nobody issues a report with a mangled name.
    """
    bad = set()
    for key, _label in PATIENT_FIELDS:
        bad |= unsupported_characters(patient.get(key, ""))
    return sorted(bad)


def build_report(patient, run, output_dir, pdf_path):
    """
    Write the PDF. `run` is the job record (parameters, timings, status).

    RETURNS:
        dict describing what went into the report, for display in the UI.
    """
    manifest = load_pipeline_manifest(output_dir)
    pcgr = find_pcgr_outputs(output_dir)
    headers, variants, tsv_source = read_pcgr_variants(pcgr["tsv"])

    label = patient.get("patient_id") or patient.get("specimen_id") or "unidentified"
    doc = PDFDocument(
        title=f"Somatic variant report - {label}",
        footer="Research use only - not a validated diagnostic report",
    )

    doc.heading("Somatic Variant Analysis Report", 1)
    doc.paragraph(
        "Generated " +
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M") +
        " by the cancer DNA pipeline web interface.",
        size=9, grey=0.4)

    # --- Patient -------------------------------------------------------
    doc.heading("Patient and specimen", 2)
    doc.key_values([(lbl, patient.get(key, ""))
                    for key, lbl in PATIENT_FIELDS])

    # --- How it was processed ------------------------------------------
    doc.heading("Analysis provenance", 2)
    params = manifest.get("parameters") or {}
    provenance = [
        ("Tumour sample", manifest.get("tumour") or run.get("tumour_sample")),
        ("Normal sample", manifest.get("normal") or run.get("normal_sample")
         or "none (tumour-only)"),
        ("Reference", manifest.get("reference") or params.get("reference")),
        ("COSMIC", params.get("cosmic")),
        ("dbSNP", params.get("dbsnp")),
        ("Germline resource", params.get("germline_resource")),
        ("Panel of normals", params.get("panel_of_normals")),
        ("Contamination resource", params.get("contamination_resource")),
        ("Target intervals", params.get("intervals")),
        ("Minimum depth", params.get("min_depth")),
        ("Panel size for TMB (Mb)", params.get("pcgr_target_size_mb")),
        ("Pipeline started", manifest.get("started_utc") or run.get("started")),
        ("Pipeline finished", manifest.get("finished_utc")
         or run.get("finished")),
        ("Duration (s)", manifest.get("duration_seconds")),
        ("Steps completed", ", ".join(manifest.get("steps_completed") or [])),
    ]
    doc.key_values([(k, v) for k, v in provenance if v not in (None, "")])

    versions = manifest.get("tool_versions") or {}
    if versions:
        doc.heading("Tool versions", 3)
        doc.table(["Tool", "Version"],
                  [[k, v] for k, v in sorted(versions.items())],
                  widths=[1, 3])

    # --- What the run observed -----------------------------------------
    summary = summarise_run(manifest)
    doc.heading("Variant summary", 2)
    if summary:
        doc.key_values(summary)
    else:
        doc.paragraph(
            "No VCF statistics were recorded in the pipeline manifest. This "
            "is expected when bcftools was unavailable or the calling step "
            "did not complete.", size=9, grey=0.35)

    # --- PCGR ------------------------------------------------------------
    doc.heading("Clinical interpretation (PCGR)", 2)
    if variants:
        doc.paragraph(
            f"Extracted from {os.path.basename(tsv_source)}. This is a "
            f"convenience extract; PCGR's own HTML report is the "
            f"authoritative interpretation.", size=9, grey=0.35)
        doc.spacer(4)
        doc.table(headers, variants,
                  widths=[1.2, 2.2, 1.6, 0.8, 0.8][:len(headers)],
                  max_rows=40)
    elif pcgr["html"]:
        doc.paragraph(
            "A PCGR report was produced but no variant table could be "
            "extracted from its TSV output (the column names may differ in "
            "this PCGR release). See the HTML report for the full "
            "interpretation.", size=9, grey=0.35)
    else:
        doc.paragraph(
            "PCGR did not run for this analysis, so no actionability tiers, "
            "TMB or mutational-signature results are available. Re-run with "
            "--pcgr-refdata-dir to enable it.", size=9, grey=0.35)

    if pcgr["html"]:
        doc.heading("Full report files", 3)
        doc.key_values([("PCGR report", os.path.basename(pcgr["html"][0])),
                        ("Location", pcgr["dir"])])

    # --- Where everything lives -----------------------------------------
    doc.heading("Output files", 2)
    files = [
        ("Filtered VCF", manifest.get("filtered_vcf")),
        ("COSMIC-annotated VCF", manifest.get("cosmic_vcf")),
        ("SnpEff-annotated VCF", manifest.get("annotated_vcf")),
        ("Output directory", output_dir),
    ]
    doc.key_values([(k, v) for k, v in files if v])

    # --- The bit that stops this being mistaken for a diagnostic report --
    doc.heading("Limitations and intended use", 2)
    for line in (
        "This report summarises an automated research pipeline. It is not a "
        "validated diagnostic test and must not be used as the sole basis "
        "for a clinical decision.",
        "Variant calls are somatic SNVs and short indels only. Copy-number "
        "changes, structural variants and gene fusions are not assessed.",
        "Any therapeutic associations shown originate from PCGR's knowledge "
        "bases and require review by a qualified molecular pathologist "
        "against the current literature.",
        "Tumour mutational burden is only meaningful when depth and allele "
        "fraction were available to PCGR; see the pipeline's --pcgr-lift-tags "
        "option and the run log.",
    ):
        doc.paragraph("- " + line, size=9)

    doc.save(pdf_path)
    return {
        "pdf": pdf_path,
        "variants_in_table": len(variants),
        "pcgr_html": pcgr["html"][0] if pcgr["html"] else None,
        "tsv_source": tsv_source,
        "manifest_found": bool(manifest),
    }
