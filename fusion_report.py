#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
fusion_report.py
================
Turn a fusion caller's table into something a human can act on.

WHY THIS EXISTS AT ALL
----------------------
  PCGR does not do fusions. It is a somatic SNV/indel reporter, it takes a
  VCF, and there is no equivalent of it for rearrangements -- so the RNA
  branch has no clinical interpretation layer unless one is written. This
  is that layer, and it is deliberately modest: it sorts, groups and
  annotates what the caller found. It does not tier evidence, it does not
  match therapies, and it does not consult a knowledge base.

WHAT IT ADDS TO THE CALLER'S TSV
--------------------------------
  1. THE ACTIONABLE LIST IS SEPARATED OUT. Arriba's output is ordered by
     its own confidence, which is a statement about the EVIDENCE, not about
     the clinical meaning. A high-confidence read-through between two
     housekeeping genes sorts above a medium-confidence EML4-ALK. So
     anything involving a gene on the list below is surfaced first,
     regardless of rank.

  2. SAME-GENE EVENTS ARE LABELLED SEPARATELY. MET exon 14 skipping is one
     of the most actionable findings a lung panel produces and it is not a
     fusion between two genes -- it is a splice event within one, and in a
     table sorted by gene pair it reads as "MET--MET" and is trivially
     missed. EGFRvIII has the same shape.

  3. THE CONFIDENCE BAR IS STATED, NOT APPLIED SILENTLY. Everything the
     caller reported appears in the output; the calls below the bar are
     marked rather than dropped, because a laboratory investigating a
     clinically expected fusion needs to see that it was called at low
     confidence, not to be told nothing was found.

THE GENE LIST IS A CONVENIENCE, NOT A KNOWLEDGE BASE
----------------------------------------------------
  The list below is the set of genes whose fusions have well-established
  targeted therapies or diagnostic significance in solid tumours and
  haematological malignancy. It is a highlighting aid so that a reader does
  not scroll past ALK. It is NOT:

    * a claim that a listed fusion is actionable in this tumour type -- the
      same fusion is predictive in one tissue and incidental in another;
    * a claim that an unlisted fusion is not important;
    * a substitute for a curated database or for a molecular pathologist.

  Keep it in mind when the report marks something "actionable gene".

Usage
-----
    python fusion_report.py --fusions rna_out/fusions/S1.fusions.tsv \\
        --sample-id S1 --output-dir rna_out/fusion_report
"""

import argparse
import csv
import html
import json
import os
import sys
from datetime import datetime, timezone

# Genes whose rearrangements carry established therapeutic or diagnostic
# weight in oncology. Grouped by why they are here, so the list can be
# maintained by somebody who was not there when it was written.
ACTIONABLE_GENES = {
    # Receptor tyrosine kinases with approved fusion-directed inhibitors.
    "ALK", "ROS1", "RET", "NTRK1", "NTRK2", "NTRK3", "MET", "FGFR1",
    "FGFR2", "FGFR3", "NRG1", "EGFR", "ERBB2", "PDGFRA", "PDGFRB",
    "BRAF", "RAF1", "ABL1", "ABL2", "KIT", "FLT3", "JAK2", "CSF1R",
    "LTK", "NUTM1", "AXL",
    # Recurrent diagnostic fusion partners -- sarcoma, leukaemia and
    # lymphoma, where the fusion defines the entity rather than the
    # treatment.
    "EWSR1", "FUS", "SS18", "PAX3", "PAX7", "FOXO1", "DDIT3", "TFE3",
    "TFEB", "MYB", "MYBL1", "NUTM2A", "BCR", "KMT2A", "RUNX1", "ETV6",
    "PML", "RARA", "CBFB", "MYH11", "NPM1", "TCF3", "PBX1", "CIC",
    "BCOR", "STAT6", "NAB2", "USP6", "GLI1", "HMGA2", "PLAG1", "MDM2",
}

# Same-gene events that are clinically meaningful in their own right, and
# that a gene-pair table hides. The report names these explicitly.
SPLICE_EVENTS = {
    "MET": "MET exon 14 skipping -- an established indication for MET "
           "inhibitors in NSCLC. Confirm the breakpoints fall either side "
           "of exon 14 before reporting it.",
    "EGFR": "EGFR internal rearrangement (EGFRvIII and similar). Check the "
            "exon boundaries; a same-gene EGFR call is not automatically "
            "vIII.",
    "AR": "AR splice variant (AR-V7 and similar), relevant in prostate "
          "cancer.",
}

CONFIDENCE_ORDER = {"high": 0, "medium": 1, "low": 2}


def _split_genes(value):
    """
    Gene names out of one of Arriba's gene columns.

    Arriba writes the gene at a breakpoint as a name, and where the
    breakpoint is intergenic it writes the flanking genes separated by a
    comma or a parenthesised distance -- "GENE1(1234),GENE2". Everything
    before the first '(' of each comma-separated part is the symbol.
    """
    genes = []
    for part in str(value or "").split(","):
        name = part.split("(")[0].strip()
        if name and name not in (".", "-"):
            genes.append(name)
    return genes


def load_arriba(path):
    """
    Read an Arriba fusions.tsv into a list of dicts.

    Arriba's header line begins with '#gene1' in some releases and 'gene1'
    in others; both are handled rather than matched exactly, because a
    header that fails to parse would produce an empty report that looks
    like a negative.
    """
    rows = []
    try:
        with open(path, encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = None
            for fields in reader:
                if not fields:
                    continue
                if header is None:
                    header = [f.lstrip("#").strip() for f in fields]
                    continue
                if len(fields) < len(header):
                    fields += [""] * (len(header) - len(fields))
                rows.append(dict(zip(header, fields)))
    except OSError:
        return None
    return rows


def classify(row, min_confidence="medium"):
    """
    Everything this report adds to one caller row.

    Returns the row with the added keys, rather than a parallel structure,
    so the HTML and the JSON cannot disagree about what was decided.
    """
    gene1 = _split_genes(row.get("gene1"))
    gene2 = _split_genes(row.get("gene2"))
    genes = gene1 + gene2
    confidence = str(row.get("confidence", "")).strip().lower()

    hit = sorted({g for g in genes if g in ACTIONABLE_GENES})
    same_gene = bool(gene1 and gene2 and gene1[0] == gene2[0])
    splice_note = SPLICE_EVENTS.get(gene1[0]) if same_gene and gene1 else None

    bar = CONFIDENCE_ORDER.get(str(min_confidence).lower(), 1)
    rank = CONFIDENCE_ORDER.get(confidence, 3)

    row["_genes"] = genes
    row["_pair"] = f"{gene1[0] if gene1 else '?'}--{gene2[0] if gene2 else '?'}"
    row["_actionable_genes"] = hit
    row["_same_gene_event"] = same_gene
    row["_splice_note"] = splice_note
    row["_confidence"] = confidence or "unknown"
    row["_meets_confidence"] = rank <= bar
    # Supporting reads: split reads either side of the breakpoint PLUS
    # discordant mates -- total evidence for the event, which is what a
    # reader scanning a table wants as one number.
    #
    # NOTE the deliberate difference from write_pcgr_fusion_tsv(), which
    # sums only the split reads. PCGR's column is named SplitReads and its
    # --fusion_min_split_reads threshold is written against split reads, so
    # feeding it this total would silently inflate every fusion past that
    # threshold. Two questions, two numbers; the column headings say which
    # is which.
    total = 0
    for column in ("split_reads1", "split_reads2", "discordant_mates"):
        try:
            total += int(str(row.get(column, "0")).strip() or 0)
        except ValueError:
            pass
    row["_supporting_reads"] = total
    return row


def order_key(row):
    """
    Sort order: actionable first, then confidence, then read support.

    This is the whole point of the report. The caller's own order is a
    statement about evidence quality; a clinical reader needs the
    actionable calls at the top even when the evidence for them is
    middling, because that is the call they will act on or investigate.
    """
    return (
        0 if row["_actionable_genes"] else 1,
        0 if row["_splice_note"] else 1,
        CONFIDENCE_ORDER.get(row["_confidence"], 3),
        -row["_supporting_reads"],
    )


def summarise(rows, min_confidence):
    """Counts for the top of the report and for the run manifest."""
    reportable = [r for r in rows if r["_meets_confidence"]]
    actionable = [r for r in reportable if r["_actionable_genes"]]
    splice = [r for r in reportable if r["_splice_note"]]
    by_confidence = {}
    for row in rows:
        by_confidence[row["_confidence"]] = \
            by_confidence.get(row["_confidence"], 0) + 1
    return {
        "total_called": len(rows),
        "at_or_above_confidence": len(reportable),
        "min_confidence": min_confidence,
        "involving_actionable_gene": len(actionable),
        "same_gene_splice_events": len(splice),
        "by_confidence": by_confidence,
        "actionable_pairs": sorted({r["_pair"] for r in actionable}),
    }


def write_html(path, sample_id, rows, summary, context):
    """One self-contained page; no external CSS, for the usual reasons."""
    def esc(value):
        return html.escape(str(value))

    def row_html(row):
        tint = {"high": "#1a7f37", "medium": "#b7791f",
                "low": "#57606a"}.get(row["_confidence"], "#57606a")
        tags = []
        if row["_actionable_genes"]:
            tags.append(
                f"<span style='background:#e7f5ec;color:#1a7f37;"
                f"padding:1px 6px;border-radius:10px;font-size:11px'>"
                f"actionable gene: "
                f"{esc(', '.join(row['_actionable_genes']))}</span>")
        if row["_splice_note"]:
            tags.append(
                "<span style='background:#fdf3e0;color:#b7791f;"
                "padding:1px 6px;border-radius:10px;font-size:11px'>"
                "same-gene splice event</span>")
        if not row["_meets_confidence"]:
            tags.append(
                "<span style='background:#f3f4f6;color:#57606a;"
                "padding:1px 6px;border-radius:10px;font-size:11px'>"
                "below the confidence bar</span>")
        note = (f"<div style='color:#b7791f;font-size:12px;margin-top:4px'>"
                f"{esc(row['_splice_note'])}</div>"
                if row["_splice_note"] else "")
        return (
            f"<tr>"
            f"<td><strong>{esc(row['_pair'])}</strong>"
            f"<div style='margin-top:3px'>{' '.join(tags)}</div>{note}</td>"
            f"<td style='color:{tint}'>{esc(row['_confidence'])}</td>"
            f"<td>{esc(row['_supporting_reads'])}</td>"
            f"<td class='mono'>{esc(row.get('breakpoint1', ''))}<br>"
            f"{esc(row.get('breakpoint2', ''))}</td>"
            f"<td>{esc(row.get('type', ''))}<br>"
            f"<span class='muted'>{esc(row.get('reading_frame', ''))}</span>"
            f"</td></tr>")

    body = "".join(row_html(r) for r in rows)
    pairs = ", ".join(summary["actionable_pairs"]) or "none"

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Fusion report — {esc(sample_id)}</title>
<style>
 body {{ font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
        margin: 0; padding: 24px; color: #1f2328; background: #f6f8fa; }}
 main {{ max-width: 1000px; margin: 0 auto; }}
 h1 {{ font-size: 20px; margin: 0 0 4px; }}
 h2 {{ font-size: 14px; margin: 0 0 10px; }}
 .card {{ background: #fff; border: 1px solid #d0d7de; border-radius: 8px;
          padding: 16px 18px; margin-bottom: 16px; }}
 table {{ width: 100%; border-collapse: collapse; }}
 td, th {{ border-bottom: 1px solid #eaeef2; padding: 8px;
           text-align: left; vertical-align: top; font-size: 13px; }}
 th {{ color: #57606a; font-size: 12px; text-transform: uppercase;
       letter-spacing: .03em; }}
 .muted {{ color: #57606a; font-size: 12.5px; }}
 .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 12px; }}
</style></head><body><main>
<h1>Fusion report — {esc(sample_id)}</h1>
<p class="muted">{esc(context.get('generated', ''))} ·
 {esc(context.get('panel', 'no panel named'))} ·
 caller: {esc(context.get('caller', 'unknown'))}</p>

<div class="card">
 <h2>What was found</h2>
 <table><tbody>
  <tr><td>Calls reported by the caller</td>
      <td>{summary['total_called']}</td></tr>
  <tr><td>At or above the '{esc(summary['min_confidence'])}' bar</td>
      <td>{summary['at_or_above_confidence']}</td></tr>
  <tr><td>Involving a gene on the actionable list</td>
      <td>{summary['involving_actionable_gene']} ({esc(pairs)})</td></tr>
  <tr><td>Same-gene splice events</td>
      <td>{summary['same_gene_splice_events']}</td></tr>
 </tbody></table>
 <p class="muted" style="margin-top:10px">An empty table is only a negative
  result if the library was adequate. Read the RNA QC report beside this
  one before reporting one.</p>
</div>

<div class="card">
 <h2>Calls</h2>
 <table>
  <thead><tr><th>Fusion</th><th>Confidence</th>
   <th title="Split reads either side of the breakpoint plus discordant mates">Supporting reads</th>
   <th>Breakpoints</th><th>Type / frame</th></tr></thead>
  <tbody>{body or
   '<tr><td colspan="5" class="muted">No fusions were called.</td></tr>'}
  </tbody>
 </table>
</div>

<div class="card">
 <h2>How to read this</h2>
 <p class="muted">Calls are ordered by clinical salience, not by the
  caller's own ranking: anything involving a gene with established
  fusion-directed therapy or diagnostic significance comes first, then
  same-gene splice events, then confidence, then read support. Everything
  the caller reported is shown; calls below the confidence bar are marked
  rather than hidden, because a laboratory investigating a clinically
  expected fusion needs to see that it was called weakly rather than to be
  told nothing was found.</p>
 <p class="muted"><strong>The actionable-gene list is a highlighting aid,
  not a knowledge base.</strong> It does not claim that a marked fusion is
  actionable in this tumour type -- the same fusion can be predictive in
  one tissue and incidental in another -- and it does not claim that an
  unmarked fusion is unimportant. Interpretation belongs to a qualified
  molecular pathologist against current evidence.</p>
 <p class="muted"><strong>Supporting reads</strong> above counts split
  reads on both sides of the breakpoint plus discordant mates. The file
  handed to PCGR counts only the split reads, because PCGR's own
  <code>SplitReads</code> threshold is defined that way \u2014 so the two
  numbers for one fusion differ, and deliberately.</p>
 <p class="muted">Fusion calls from RNA are evidence of a transcript, not
  proof of a genomic rearrangement, and callers disagree substantially on
  real data. Orthogonal confirmation (FISH, RT-PCR, or a second caller) is
  expected before clinical reporting.</p>
</div>
</main></body></html>
""")
    return path


def write_tsv(path, rows):
    """The annotated table, for anyone who would rather have a spreadsheet."""
    columns = ["pair", "confidence", "meets_confidence",
               "supporting_reads_incl_discordant",
               "actionable_genes", "same_gene_event", "breakpoint1",
               "breakpoint2", "type", "reading_frame", "site1", "site2"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([
                row["_pair"], row["_confidence"],
                "yes" if row["_meets_confidence"] else "no",
                row["_supporting_reads"],
                ",".join(row["_actionable_genes"]),
                "yes" if row["_same_gene_event"] else "no",
                row.get("breakpoint1", ""), row.get("breakpoint2", ""),
                row.get("type", ""), row.get("reading_frame", ""),
                row.get("site1", ""), row.get("site2", ""),
            ])
    return path


# =============================================================================
# HANDOFF TO PCGR
# =============================================================================
# PCGR 2.x accepts RNA fusions as a first-class molecular input alongside
# (or instead of) a somatic VCF, which means the fusions from an RNA
# library and the SNVs from the DNA library of the same specimen can be
# interpreted in ONE clinical report rather than two that a reader has to
# reconcile.
#
# Its schema, taken from pcgr.validate.is_valid_rna_fusion() in the
# installed version rather than from documentation:
#
#     FusionGene        "GENE1--GENE2" or "GENE1::GENE2", exactly two parts,
#                       each a gene symbol and not all digits
#     LeftBreakpoint    "<chrom>:<pos>", with or without the chr prefix
#     RightBreakpoint   same
#     SplitReads        non-negative integer
#     Score             optional, float
#
# Arriba's own columns map onto it cleanly, with two judgement calls
# recorded below.

PCGR_FUSION_COLUMNS = ["FusionGene", "LeftBreakpoint", "RightBreakpoint",
                       "SplitReads"]


def write_pcgr_fusion_tsv(rows, path, min_confidence=None):
    """
    Convert classified Arriba rows into PCGR's fusion input format.

    TWO JUDGEMENT CALLS, both deliberate:

    SplitReads is split_reads1 + split_reads2. Arriba counts the reads
    crossing the breakpoint on each side separately; PCGR wants one number
    for "reads supporting this fusion", and the sum is that number.
    Discordant mates are NOT added: they support the pairing without
    crossing the junction, and PCGR's own --fusion_min_split_reads
    threshold is written against split reads.

    Score is left OUT. It is optional and typed float, and Arriba grades
    confidence as high/medium/low -- mapping those onto 1.0/0.5/0.1 would
    invent a number with no meaning and PCGR would then rank on it.

    Rows below `min_confidence` are excluded when it is given: unlike the
    HTML report, which marks them, a weak call handed to PCGR becomes an
    entry in a clinical interpretation, and the report has no way to show
    that it was weak.

    RETURNS:
        (path, written, skipped) -- or (None, 0, 0) if nothing qualified.
        PCGR REJECTS AN EMPTY FUSION FILE outright, so "no fusions" must
        produce no file rather than an empty one.
    """
    bar = CONFIDENCE_ORDER.get(str(min_confidence).lower(), None) \
        if min_confidence else None

    written = []
    skipped = 0
    for row in rows:
        if bar is not None and \
                CONFIDENCE_ORDER.get(row.get("_confidence"), 3) > bar:
            skipped += 1
            continue
        genes = _split_genes(row.get("gene1")), _split_genes(row.get("gene2"))
        if not genes[0] or not genes[1]:
            skipped += 1
            continue
        left, right = genes[0][0], genes[1][0]
        # PCGR rejects an all-digit gene name and requires two parts; an
        # intergenic breakpoint that yielded neither is dropped rather than
        # passed through to fail validation hours later.
        if not left or not right or left.isdigit() or right.isdigit():
            skipped += 1
            continue
        bp1 = str(row.get("breakpoint1", "")).strip()
        bp2 = str(row.get("breakpoint2", "")).strip()
        if ":" not in bp1 or ":" not in bp2:
            skipped += 1
            continue
        support = 0
        for column in ("split_reads1", "split_reads2"):
            try:
                support += int(str(row.get(column, "0")).strip() or 0)
            except ValueError:
                pass
        written.append([f"{left}--{right}", bp1, bp2, str(max(support, 0))])

    if not written:
        return None, 0, skipped

    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(PCGR_FUSION_COLUMNS)
        writer.writerows(written)
    return path, len(written), skipped


def run_fusion_report(fusions_tsv, sample_id, output_dir,
                      min_confidence="medium", panel_label=None,
                      caller="arriba", dry_run=False):
    """
    Build the report. Importable, so fusion_calling.py runs it as a step.

    Returns a dict and never raises for a file it could not read: a missing
    report is a caveat on a run, not a reason to discard the calls that are
    already on disk.
    """
    os.makedirs(output_dir, exist_ok=True)
    html_path = os.path.join(output_dir, f"{sample_id}.fusions.html")
    json_path = os.path.join(output_dir, f"{sample_id}.fusions.json")
    tsv_path = os.path.join(output_dir, f"{sample_id}.fusions.annotated.tsv")

    if dry_run:
        return {"ok": True, "html": html_path, "json": json_path,
                "tsv": tsv_path,
                "pcgr_tsv": os.path.join(output_dir,
                                         f"{sample_id}.pcgr_fusions.tsv"),
                "summary": {}, "error": None}

    raw = load_arriba(fusions_tsv)
    if raw is None:
        return {"ok": False, "html": None, "json": None, "tsv": None,
                "summary": {},
                "error": f"could not read {fusions_tsv}"}

    rows = sorted((classify(r, min_confidence) for r in raw), key=order_key)
    summary = summarise(rows, min_confidence)
    context = {
        "generated": datetime.now(timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"),
        "panel": panel_label or "no panel named",
        "caller": caller,
        "source": os.path.abspath(fusions_tsv),
    }
    pcgr_path = os.path.join(output_dir, f"{sample_id}.pcgr_fusions.tsv")
    try:
        write_html(html_path, sample_id, rows, summary, context)
        write_tsv(tsv_path, rows)
        # The PCGR handoff, written whether or not PCGR will be run: it
        # costs nothing, and a laboratory that decides later to produce a
        # combined DNA+RNA report should not have to re-run the pipeline
        # to get the input file for it.
        pcgr_path, pcgr_written, pcgr_skipped = write_pcgr_fusion_tsv(
            rows, pcgr_path, min_confidence=min_confidence)
        summary["pcgr_fusions_written"] = pcgr_written
        summary["pcgr_fusions_skipped"] = pcgr_skipped
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump({"sample": sample_id, "context": context,
                       "summary": summary, "fusions": rows}, fh, indent=2)
    except OSError as exc:
        return {"ok": False, "html": None, "json": None, "tsv": None,
                "pcgr_tsv": None, "summary": summary, "error": str(exc)}
    return {"ok": True, "html": html_path, "json": json_path,
            "tsv": tsv_path, "pcgr_tsv": pcgr_path, "summary": summary,
            "error": None}


def main():
    parser = argparse.ArgumentParser(
        description="Annotate and rank a fusion caller's output for a "
                    "clinical reader.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--fusions", required=True,
                        help="Arriba fusions.tsv.")
    parser.add_argument("--sample-id", required=True,
                        help="Sample name; used for the output filenames.")
    parser.add_argument("--output-dir", required=True,
                        help="Where the HTML/JSON/TSV reports are written.")
    parser.add_argument("--min-confidence", default="medium",
                        choices=["low", "medium", "high"],
                        help="The bar at which a call counts as reportable. "
                             "Calls below it are MARKED, never dropped.")
    parser.add_argument("--panel-label", default=None,
                        help="Assay name, recorded on the report.")
    parser.add_argument("--caller", default="arriba",
                        help="Which caller produced the input, for the "
                             "record.")
    args = parser.parse_args()

    result = run_fusion_report(
        args.fusions, args.sample_id, args.output_dir,
        min_confidence=args.min_confidence, panel_label=args.panel_label,
        caller=args.caller)
    if not result["ok"]:
        print(f"[WARN] Fusion report not produced: {result['error']}")
        return 1
    summary = result["summary"]
    print(f"[OK] Fusion report: {result['html']}")
    print(f"[OK] {summary['total_called']} call(s); "
          f"{summary['at_or_above_confidence']} at or above "
          f"'{summary['min_confidence']}'; "
          f"{summary['involving_actionable_gene']} involving an actionable "
          f"gene.")
    if result.get("pcgr_tsv"):
        print(f"[OK] PCGR fusion input: {result['pcgr_tsv']} "
              f"({summary.get('pcgr_fusions_written', 0)} fusion(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
