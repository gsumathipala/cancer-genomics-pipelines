#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
rna_qc_report.py
================
What an RNA run's negative result is entitled to say.

THE PROBLEM THIS EXISTS FOR
---------------------------
  coverage_report.py answers that question on the DNA side: a region the
  sequencing never reached and a region that is wild type both appear as
  "no variant", and only a coverage statement tells them apart.

  RNA has the same hole and it is worse, because the failure is
  library-wide rather than regional. A degraded FFPE extraction, a failed
  library prep, or an RNA input below the assay's validated DV200 produces:

      * a BAM that opens
      * a mapping rate that looks only somewhat low
      * a fusion table that is clean, well-formatted and EMPTY

  which is indistinguishable, in every downstream artefact, from a
  specimen that genuinely carries no fusion. A report that says "no fusions
  detected" on that library is not a negative result. It is no result,
  reported as a negative one.

  So this measures the library and says plainly whether a negative is
  interpretable.

WHAT IT MEASURES, AND FROM WHAT
-------------------------------
  Deliberately cheap, and deliberately dependent on nothing but samtools
  and files the run has already produced -- the same constraint that keeps
  coverage_report.py usable:

    STAR's Log.final.out    input reads, uniquely mapped %, multi-mapping,
                            reads lost to "too short" (the signature of a
                            degraded or contaminated library), splice
                            junction counts, chimeric reads.
    samtools idxstats       mitochondrial fraction. High chrM with low
                            exonic signal is what a degraded FFPE RNA
                            library looks like.
    the GTF                 rRNA gene intervals, derived once, so the rRNA
                            fraction can be measured with samtools. On a
                            poly-A or rRNA-depleted library a high rRNA
                            fraction means the depletion failed and most of
                            the reads are carrying no information.

WHAT IT DOES NOT MEASURE, AND WHY NOT
-------------------------------------
  * DUPLICATION RATE. On RNA this is not a quality metric: a highly
    expressed gene produces identical reads legitimately, and the
    "duplication rate" of a good amplicon RNA library approaches 100%.
    Reporting it would invite somebody to act on it.
  * DV200 / RIN. These are properties of the extracted RNA measured on an
    instrument before sequencing. They cannot be recovered from FASTQ, and
    they are the single best predictor of whether fusion detection will
    work at all -- so they are a field the operator records, not a number
    this computes. The web form asks for it.
  * GENE-BODY COVERAGE and STRANDEDNESS. Both want a transcript model and
    a per-base pass that this deliberately avoids.

THRESHOLDS ARE CONVENTIONAL, NOT VALIDATED
------------------------------------------
  The pass/fail marks come from the panel profile, whose defaults are
  ordinary practice for the chemistry -- not from a validation of this
  pipeline on your assay. They are a prompt to look, not a release
  criterion.

Usage
-----
    python rna_qc_report.py --sample-id S1 \\
        --star-log rna_out/rna_aligned/S1.Log.final.out \\
        --bam rna_out/rna_aligned/S1.sorted.bam \\
        --gtf ~/data/references/gencode/gencode.v50.primary_assembly.annotation.gtf \\
        --output-dir rna_out/rna_qc
"""

import argparse
import html
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

# Above this fraction of reads on the mitochondrial contig, a library is
# dominated by mitochondrial transcripts -- which on a nuclear-gene fusion
# panel means the informative fraction is small. Conventional, not
# validated; see the module docstring.
MT_FRACTION_WARN = 0.30

# Above this fraction of reads in rRNA genes, depletion or selection has
# failed and most of the library carries no information about the panel.
RRNA_FRACTION_WARN = 0.30

# STAR's "% of reads unmapped: too short" is the degraded-library
# signature: fragments too short to align across their own length. On FFPE
# RNA a modest figure is normal and a large one is not.
TOO_SHORT_WARN = 30.0


def run(cmd, dry_run=False):
    """Run a command and return stdout, or None."""
    if dry_run:
        print("$ " + " ".join(str(c) for c in cmd))
        return None
    try:
        out = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        return None
    if out.returncode != 0:
        return None
    return out.stdout


# =============================================================================
# SECTION 1: MEASUREMENTS
# =============================================================================

def rrna_intervals(gtf, destination):
    """
    Write a BED of rRNA gene intervals, derived from the GTF.

    Derived rather than downloaded: the annotation already says which genes
    are rRNA, and a separate interval file is one more thing that can be
    the wrong release. Matches gene_type "rRNA" and "rRNA_pseudogene",
    which is how GENCODE labels them.

    RETURNS:
        (path, count) or (None, 0) when the GTF is unreadable or has none.
    """
    count = 0
    try:
        with open(gtf, encoding="utf-8", errors="replace") as src, \
                open(destination, "w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attributes = fields[8]
                if 'gene_type "rRNA' not in attributes and \
                        'gene_biotype "rRNA' not in attributes:
                    continue
                # GTF is 1-based inclusive; BED is 0-based half-open.
                dst.write(f"{fields[0]}\t{int(fields[3]) - 1}\t{fields[4]}\n")
                count += 1
    except (OSError, ValueError):
        return None, 0
    return (destination, count) if count else (None, 0)


def idxstats_summary(bam, dry_run=False):
    """
    Mapped reads in total and on the mitochondrial contig.

    samtools idxstats is a read of the BAM index rather than of the BAM, so
    this costs nothing even on a large file.
    """
    out = run(["samtools", "idxstats", bam], dry_run=dry_run)
    if not out:
        return {}
    total = 0
    mito = 0
    for line in out.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            mapped = int(fields[2])
        except ValueError:
            continue
        total += mapped
        # Both conventions, because a BAM may have come from either
        # reference naming -- the same trap the panel BEDs have.
        if fields[0] in ("chrM", "MT", "M", "chrMT"):
            mito += mapped
    if not total:
        return {}
    return {"mapped_reads": total, "mito_reads": mito,
            "mito_fraction": round(mito / total, 4)}


def region_read_count(bam, bed, dry_run=False):
    """Reads overlapping a BED, via samtools view -c -L."""
    out = run(["samtools", "view", "-c", "-L", bed, bam], dry_run=dry_run)
    if out is None:
        return None
    try:
        return int(out.strip())
    except ValueError:
        return None


def collect(sample_id, star_log, bam, gtf, output_dir, thresholds,
            dry_run=False):
    """
    Every measurement this report makes, as one dict.

    Nothing here is fatal. Each measurement that cannot be taken is simply
    absent and is reported as "not measured" rather than as a zero -- a
    zero would read as a finding.
    """
    os.makedirs(output_dir, exist_ok=True)
    metrics = {"sample": sample_id, "measured": {}, "not_measured": []}

    # --- STAR's own report ---
    star_stats = {}
    if star_log and os.path.exists(star_log):
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from align_rna import parse_star_log
            star_stats = parse_star_log(star_log)
        except ImportError:
            metrics["not_measured"].append(
                "STAR alignment statistics (align_rna.py not importable)")
    else:
        metrics["not_measured"].append(
            "STAR alignment statistics (no Log.final.out given)")

    if star_stats:
        input_reads = star_stats.get("Number of input reads") or 0
        metrics["measured"].update({
            "input_reads": input_reads,
            "input_reads_millions": round(input_reads / 1e6, 2)
            if input_reads else 0,
            "uniquely_mapped_pct":
                star_stats.get("Uniquely mapped reads %"),
            "multimapping_pct":
                star_stats.get("% of reads mapped to multiple loci"),
            "too_short_pct":
                star_stats.get("% of reads unmapped: too short"),
            "splice_junctions_total":
                star_stats.get("Number of splices: Total"),
            "splice_junctions_annotated":
                star_stats.get("Number of splices: Annotated (sjdb)"),
            "chimeric_reads": star_stats.get("Number of chimeric reads"),
            "chimeric_pct": star_stats.get("% of chimeric reads"),
            "mean_read_length":
                star_stats.get("Average input read length"),
        })

    # --- the BAM ---
    if bam and (dry_run or os.path.exists(bam)):
        metrics["measured"].update(idxstats_summary(bam, dry_run=dry_run))
        if gtf and os.path.exists(gtf):
            bed = os.path.join(output_dir, f"{sample_id}.rrna.bed")
            path, genes = rrna_intervals(gtf, bed)
            if path:
                count = region_read_count(bam, path, dry_run=dry_run)
                mapped = metrics["measured"].get("mapped_reads")
                if count is not None and mapped:
                    metrics["measured"]["rrna_reads"] = count
                    metrics["measured"]["rrna_fraction"] = round(
                        count / mapped, 4)
                    metrics["measured"]["rrna_genes_in_annotation"] = genes
                else:
                    metrics["not_measured"].append(
                        "rRNA fraction (samtools could not count it)")
            else:
                metrics["not_measured"].append(
                    "rRNA fraction (no rRNA genes found in the annotation)")
        else:
            metrics["not_measured"].append("rRNA fraction (no GTF given)")
    else:
        metrics["not_measured"].append(
            "mitochondrial and rRNA fractions (no indexed BAM given)")

    # Never computed, always explained -- see the module docstring.
    metrics["not_measured"].append(
        "duplication rate (not a quality metric on RNA: an abundant "
        "transcript produces identical reads legitimately)")
    metrics["not_measured"].append(
        "DV200 / RIN (an instrument measurement of the extracted RNA; it "
        "cannot be recovered from FASTQ and it is the best predictor of "
        "whether fusion detection will work -- record it with the run)")

    metrics["checks"] = evaluate(metrics["measured"], thresholds)
    metrics["verdict"] = verdict(metrics["checks"])
    return metrics


def evaluate(measured, thresholds):
    """
    Turn the measurements into pass / concern lines.

    Each check states the number, the bar and what a failure MEANS, because
    "uniquely mapped 41%" on its own tells a reader nothing about whether
    they may report a negative.
    """
    checks = []

    def add(name, ok, detail):
        checks.append({"name": name, "state": "pass" if ok else "concern",
                       "detail": detail})

    def unjudged(name, detail):
        # Measured, but with no bar to measure it against. Stated rather
        # than dropped: see the comment below.
        checks.append({"name": name, "state": "unassessed",
                       "detail": detail})

    # Depth and mapping rate are the two checks that most decide whether an
    # empty table means "no fusion" -- and their floors are assay-specific,
    # so they come from the panel profile or the operator, never from a
    # guess here. When neither supplied one, these checks used to be LEFT
    # OUT, and the verdict then read "every measured check cleared its bar
    # ... an empty fusion table is a negative result" for a library of any
    # depth at all: 0.6 million reads was certified. They are now always
    # present, as "not judged" when there is no floor, and the verdict
    # refuses to call a negative on a library whose depth was not judged.
    reads_m = measured.get("input_reads_millions")
    floor_m = thresholds.get("min_reads_millions")
    if reads_m is not None and floor_m:
        add("Library size", reads_m >= floor_m,
            f"{reads_m:.1f} M reads against a floor of {floor_m:.1f} M for "
            f"this assay. Below it, a fusion present at low abundance may "
            f"simply not have been sampled, so an empty result is not a "
            f"negative.")
    elif reads_m is not None:
        unjudged("Library size",
                 f"{reads_m:.1f} M reads, but no floor is set for this "
                 f"assay (choose an RNA panel profile, or set "
                 f"--rna-min-reads-millions), so whether that is enough to "
                 f"have sampled a low-abundance fusion cannot be judged.")

    unique = measured.get("uniquely_mapped_pct")
    floor_u = thresholds.get("min_unique_mapped_pct")
    if unique is not None and floor_u:
        add("Uniquely mapped", unique >= floor_u,
            f"{unique:.1f}% against a floor of {floor_u:.1f}%. A low rate "
            f"on RNA usually means degradation, contamination, or an index "
            f"built from the wrong annotation.")
    elif unique is not None:
        unjudged("Uniquely mapped",
                 f"{unique:.1f}%, but no floor is set for this assay "
                 f"(choose an RNA panel profile, or set "
                 f"--rna-min-unique-mapped-pct), so it cannot be judged.")

    too_short = measured.get("too_short_pct")
    if too_short is not None:
        add("Reads lost as 'too short'", too_short <= TOO_SHORT_WARN,
            f"{too_short:.1f}% (concern above {TOO_SHORT_WARN:.0f}%). This "
            f"is the degraded-library signature: fragments too short to "
            f"align across their own length. Expect some on FFPE.")

    mito = measured.get("mito_fraction")
    if mito is not None:
        add("Mitochondrial fraction", mito <= MT_FRACTION_WARN,
            f"{mito * 100:.1f}% of mapped reads (concern above "
            f"{MT_FRACTION_WARN * 100:.0f}%). A library dominated by "
            f"mitochondrial transcripts has little left for the nuclear "
            f"genes a fusion panel targets.")

    rrna = measured.get("rrna_fraction")
    if rrna is not None:
        add("rRNA fraction", rrna <= RRNA_FRACTION_WARN,
            f"{rrna * 100:.1f}% of mapped reads (concern above "
            f"{RRNA_FRACTION_WARN * 100:.0f}%). A high figure means rRNA "
            f"depletion or poly-A selection did not work, and most of the "
            f"library is carrying no information.")

    junctions = measured.get("splice_junctions_total")
    if junctions is not None:
        # Not a threshold so much as a sanity check: a library with almost
        # no junction-spanning reads cannot evidence a fusion, whatever
        # else its numbers say.
        add("Splice junctions observed", junctions > 0,
            f"{junctions:,} junction-spanning alignments. A fusion is "
            f"evidenced by reads that cross a junction; a library with "
            f"none cannot produce one regardless of depth.")

    chimeric = measured.get("chimeric_reads")
    if chimeric is not None:
        # The condition is `> 0`, not `is not None`. Written the second way
        # it sat inside `if chimeric is not None` and was therefore always
        # true -- so the one failure this check exists to catch, a run with
        # chimeric detection switched off, could never raise it. An empty
        # fusion table then looked like a clean negative.
        add("Chimeric detection active", chimeric > 0,
            f"{chimeric:,} chimeric reads found by STAR. Zero here, with "
            f"an otherwise healthy library, usually means STAR ran without "
            f"chimeric detection -- in which case the fusion table is "
            f"empty for a reason that has nothing to do with the "
            f"specimen.")
    return checks


# How each check state is shown, shared by this report and the web PDF.
CHECK_MARKS = {
    "pass": ("#1a7f37", "pass"),
    "concern": ("#b7791f", "look at this"),
    "unassessed": ("#57606a", "not judged"),
}


def verdict(checks):
    """One line saying whether a negative from this library is reportable."""
    concerns = [c for c in checks if c["state"] == "concern"]
    if not checks:
        return {"state": "unknown",
                "text": "Nothing could be measured, so this run makes no "
                        "statement about library quality. Treat an empty "
                        "fusion table as uninterpreted rather than as a "
                        "negative."}
    unjudged = [c for c in checks if c["state"] == "unassessed"]
    if not concerns and unjudged:
        return {"state": "partial",
                "text": "The checks that could be judged cleared their "
                        "bar, but " + " and ".join(
                            c["name"].lower() for c in unjudged) +
                        " had no floor to be judged against. An empty "
                        "fusion table from this library cannot yet be "
                        "reported as a negative: set the floors for this "
                        "assay (an RNA panel profile carries them)."}
    if not concerns:
        return {"state": "ok",
                "text": "Every measured check cleared its bar. An empty "
                        "fusion table from this library is a negative "
                        "result, within the limits of the assay."}
    return {"state": "concern",
            "text": f"{len(concerns)} of {len(checks)} checks did not clear "
                    f"their bar. An empty fusion table from this library "
                    f"should NOT be reported as a negative without "
                    f"explaining these: " +
                    "; ".join(c["name"] for c in concerns) + "."}


# =============================================================================
# SECTION 2: OUTPUT
# =============================================================================

def write_html(path, metrics, context):
    """
    A self-contained HTML page, in the same shape as the coverage report.

    No external CSS or JavaScript: these get emailed, archived and opened
    on machines that are not this one.
    """
    def esc(value):
        return html.escape(str(value))

    state = metrics["verdict"]["state"]
    colour = {"ok": "#1a7f37", "concern": "#b7791f", "partial": "#57606a",
              "unknown": "#57606a"}.get(state, "#57606a")

    rows = []
    for key, value in metrics["measured"].items():
        if value is None:
            continue
        label = key.replace("_", " ").capitalize()
        shown = f"{value:,}" if isinstance(value, int) else value
        rows.append(f"<tr><td>{esc(label)}</td><td>{esc(shown)}</td></tr>")

    checks = []
    for check in metrics["checks"]:
        tint, mark = CHECK_MARKS.get(check["state"], CHECK_MARKS["concern"])
        checks.append(
            f"<tr><td><strong>{esc(check['name'])}</strong></td>"
            f"<td style='color:{tint};white-space:nowrap'>{esc(mark)}</td>"
            f"<td>{esc(check['detail'])}</td></tr>")

    unmeasured = "".join(f"<li>{esc(item)}</li>"
                         for item in metrics["not_measured"])

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>RNA library QC — {esc(metrics['sample'])}</title>
<style>
 body {{ font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
        margin: 0; padding: 24px; color: #1f2328; background: #f6f8fa; }}
 main {{ max-width: 900px; margin: 0 auto; }}
 h1 {{ font-size: 20px; margin: 0 0 4px; }}
 h2 {{ font-size: 14px; margin: 24px 0 8px; }}
 .card {{ background: #fff; border: 1px solid #d0d7de; border-radius: 8px;
          padding: 16px 18px; margin-bottom: 16px; }}
 .verdict {{ border-left: 4px solid {colour}; }}
 .verdict p {{ margin: 0; }}
 table {{ width: 100%; border-collapse: collapse; }}
 td, th {{ border-bottom: 1px solid #eaeef2; padding: 6px 8px;
           text-align: left; vertical-align: top; font-size: 13px; }}
 .muted {{ color: #57606a; font-size: 12.5px; }}
 ul {{ margin: 6px 0 0 18px; padding: 0; }}
</style></head><body><main>
<h1>RNA library QC — {esc(metrics['sample'])}</h1>
<p class="muted">{esc(context.get('generated', ''))} ·
 {esc(context.get('panel', 'no panel named'))}</p>

<div class="card verdict">
 <h2 style="margin-top:0">Can a negative be reported from this library?</h2>
 <p>{esc(metrics['verdict']['text'])}</p>
</div>

<div class="card">
 <h2 style="margin-top:0">Checks</h2>
 <table><tbody>{''.join(checks) or
   '<tr><td class="muted">Nothing could be checked.</td></tr>'}</tbody></table>
</div>

<div class="card">
 <h2 style="margin-top:0">Measurements</h2>
 <table><tbody>{''.join(rows) or
   '<tr><td class="muted">None.</td></tr>'}</tbody></table>
</div>

<div class="card">
 <h2 style="margin-top:0">Not measured</h2>
 <p class="muted">Named rather than omitted, so their absence is not read
  as a pass.</p>
 <ul>{unmeasured}</ul>
</div>

<p class="muted">Thresholds are conventional practice for the chemistry,
 taken from the panel profile. They are not a validation of this pipeline
 on your assay, and they are a prompt to look rather than a release
 criterion.</p>
</main></body></html>
""")
    return path


def run_rna_qc(sample_id, star_log, bam, gtf, output_dir,
               min_reads_millions=None, min_unique_mapped_pct=None,
               panel_label=None, dry_run=False):
    """
    Measure one RNA library and write the HTML/JSON reports.

    Importable, so fusion_calling.py can call it as a step in the same way
    comprehensive_variant_calling.py calls run_coverage_check(). Returns a
    dict; never raises for a measurement it could not take.
    """
    thresholds = {"min_reads_millions": min_reads_millions,
                  "min_unique_mapped_pct": min_unique_mapped_pct}
    metrics = collect(sample_id, star_log, bam, gtf, output_dir, thresholds,
                      dry_run=dry_run)
    context = {
        "generated": datetime.now(timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"),
        "panel": panel_label or "no panel named",
    }
    metrics["context"] = context
    metrics["thresholds"] = thresholds

    html_path = os.path.join(output_dir, f"{sample_id}.rna_qc.html")
    json_path = os.path.join(output_dir, f"{sample_id}.rna_qc.json")
    if dry_run:
        return {"ok": True, "html": html_path, "json": json_path,
                "metrics": metrics, "error": None}
    try:
        write_html(html_path, metrics, context)
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(metrics, fh, indent=2)
    except OSError as exc:
        return {"ok": False, "html": None, "json": None,
                "metrics": metrics, "error": str(exc)}
    return {"ok": True, "html": html_path, "json": json_path,
            "metrics": metrics, "error": None}


def main():
    parser = argparse.ArgumentParser(
        description="Measure an RNA library and say whether a negative "
                    "fusion result from it is interpretable.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--sample-id", required=True,
                        help="Sample name; used for the output filenames.")
    parser.add_argument("--star-log", default=None,
                        help="STAR's Log.final.out for this sample.")
    parser.add_argument("--bam", default=None,
                        help="Coordinate-sorted, indexed BAM. Needed for "
                             "the mitochondrial and rRNA fractions.")
    parser.add_argument("--gtf", default=None,
                        help="GENCODE GTF, used to derive rRNA intervals.")
    parser.add_argument("--output-dir", required=True,
                        help="Where the HTML/JSON reports are written.")
    parser.add_argument("--min-reads-millions", type=float, default=None,
                        help="Library-size floor for this assay. Comes from "
                             "the panel profile in a normal run.")
    parser.add_argument("--min-unique-mapped-pct", type=float, default=None,
                        help="Unique-mapping floor for this assay.")
    parser.add_argument("--panel-label", default=None,
                        help="Assay name, recorded on the report.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the samtools commands without running "
                             "them.")
    args = parser.parse_args()

    result = run_rna_qc(
        args.sample_id, args.star_log, args.bam, args.gtf, args.output_dir,
        min_reads_millions=args.min_reads_millions,
        min_unique_mapped_pct=args.min_unique_mapped_pct,
        panel_label=args.panel_label, dry_run=args.dry_run)
    if not result["ok"]:
        print(f"[WARN] RNA QC report not produced: {result['error']}")
        return 1
    print(f"[OK] RNA QC report: {result['html']}")
    print(f"[{result['metrics']['verdict']['state'].upper()}] "
          f"{result['metrics']['verdict']['text']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
