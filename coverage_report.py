#!/usr/bin/env python3
"""
coverage_report.py -- did the sequencing actually cover the regions we
report on?

WHY THIS EXISTS
---------------
A variant caller can only call what it can see. If a region of the panel
was not sequenced deeply enough -- capture dropout, a GC-rich exon, a
degraded FFPE block, too little input DNA -- no variant is emitted there,
and the VCF looks exactly like a region that was sequenced beautifully and
happens to be wild type. Both come out as "nothing reported".

The clinical report inherits that ambiguity: PCGR is given a VCF, and a
VCF has no way to say "this exon was never covered". So a false negative
reaches the report indistinguishable from a true negative, and nothing
downstream flags it. THIS is the mechanism that was missing; the pipeline
already had none.

This script answers the question directly, before the report is written:
for every region in a BED file, how much of it reached the depth the
variant filter demands, and which stretches did not. Regions that fall
short are named, with their uncovered intervals, so a reader knows exactly
which part of the panel the report cannot speak for.

WHAT IT PRODUCES
----------------
For each sample, in --output-dir:
  <sample>.coverage.html    Human-readable report, self-contained, opened
                            from the webapp in a new tab.
  <sample>.coverage.json    The same numbers for programmatic use; the
                            pipeline manifest points at this.
  <sample>.coverage.tsv     Per-region table for a spreadsheet.

DEPTH IS MEASURED THE WAY THE CALLER SEES IT
--------------------------------------------
`samtools depth` is run with the same minimum mapping and base qualities
the caller uses, so "covered" here means "covered as far as Mutect2 is
concerned", not "some reads piled up nearby". The default threshold is the
pipeline's --min-depth, because that is the depth below which a call would
have been filtered out anyway.

USAGE
    python coverage_report.py --bam T.bqsr.bam --bed panel.bed \\
        --sample-id TUMOUR_01 --output-dir out/coverage --min-depth 10
"""

import argparse
import bisect
import html
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone


# A region needs this fraction of its bases at or above the depth threshold
# before it is called fully covered. It is 1.0 deliberately: a panel region
# is reported on as a whole, and "99% covered" still means some codon was
# invisible. The summary distinguishes complete from near-complete anyway.
FULL_COVERAGE_FRACTION = 1.0


def parse_bed(path):
    """
    Read a BED file into a list of regions.

    Track/browser lines and comments are skipped. BED is half-open and
    0-based; the report speaks in 1-based inclusive coordinates because
    that is what every genome browser and every clinician's variant
    nomenclature uses, so the conversion happens here, once.
    """
    regions = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{lineno}: need at least chrom/start/end, "
                    f"got {line!r}")
            chrom, start, end = fields[0], fields[1], fields[2]
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise ValueError(
                    f"{path}:{lineno}: start and end must be integers, "
                    f"got {start!r} and {end!r}") from None
            if end <= start:
                raise ValueError(
                    f"{path}:{lineno}: empty or reversed region {line!r}")
            regions.append({
                "chrom": chrom,
                "start": start + 1,          # 1-based inclusive
                "end": end,
                "name": fields[3] if len(fields) > 3 else "",
                "length": end - start,
            })
    if not regions:
        raise ValueError(f"{path}: no usable regions")
    return regions


def index_regions(regions):
    """
    Group regions by chromosome for a bisect lookup of a position.

    samtools emits one line per base, and a panel can be millions of
    bases; scanning the region list for each of them is the difference
    between seconds and hours.
    """
    by_chrom = {}
    for idx, reg in enumerate(regions):
        by_chrom.setdefault(reg["chrom"], []).append((reg["start"],
                                                      reg["end"], idx))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    starts = {c: [r[0] for r in v] for c, v in by_chrom.items()}
    return by_chrom, starts


def depth_command(bam, bed, min_mapq, min_baseq):
    """The samtools invocation, kept in one place so the report can quote it."""
    return ["samtools", "depth", "-a", "-b", bed,
            "-Q", str(min_mapq), "-q", str(min_baseq), bam]


def measure(bam, bed, regions, min_depth, min_mapq, min_baseq,
            dry_run=False):
    """
    Per-base depth over the BED, folded into per-region statistics.

    Returns (regions_with_stats, error). Depth is streamed rather than
    collected: a whole-exome BED at one line per base is gigabytes of
    text, and none of it is needed twice.
    """
    by_chrom, starts = index_regions(regions)
    stats = [{"covered_bases": 0, "total_depth": 0, "min_depth": None,
              "low_runs": [], "_open_run": None} for _ in regions]

    cmd = depth_command(bam, bed, min_mapq, min_baseq)
    if dry_run:
        print(f"[DRY RUN] Would run: {' '.join(cmd)}")
        return None, None

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, bufsize=1)
    except OSError as exc:
        return None, f"could not run samtools depth: {exc}"

    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1]), int(parts[2])
        chrom_regions = by_chrom.get(chrom)
        if not chrom_regions:
            continue
        # Regions may overlap, so walk back from the insertion point over
        # every region that could still contain this position.
        i = bisect.bisect_right(starts[chrom], pos) - 1
        while i >= 0:
            start, end, idx = chrom_regions[i]
            if end >= pos:
                st = stats[idx]
                st["total_depth"] += depth
                if st["min_depth"] is None or depth < st["min_depth"]:
                    st["min_depth"] = depth
                if depth >= min_depth:
                    st["covered_bases"] += 1
                    if st["_open_run"] is not None:
                        st["low_runs"].append(st["_open_run"])
                        st["_open_run"] = None
                else:
                    if st["_open_run"] is None:
                        st["_open_run"] = [pos, pos]
                    else:
                        st["_open_run"][1] = pos
            i -= 1

    stderr = proc.stderr.read()
    proc.wait()
    if proc.returncode != 0:
        return None, (f"samtools depth exited {proc.returncode}: "
                      f"{stderr.strip().splitlines()[-1] if stderr.strip() else 'no message'}")

    out = []
    for reg, st in zip(regions, stats):
        if st["_open_run"] is not None:
            st["low_runs"].append(st["_open_run"])
        length = reg["length"]
        covered = st["covered_bases"]
        # samtools -a reports every position in the BED, but a contig
        # missing from the BAM header yields nothing at all: treat the
        # unseen bases as the zero-depth bases they are.
        unseen = length - (covered + sum(r[1] - r[0] + 1
                                         for r in st["low_runs"]))
        if unseen > 0:
            st["low_runs"].append([reg["start"], reg["end"]])
            st["min_depth"] = 0
        fraction = covered / length if length else 0.0
        out.append(dict(reg, **{
            "covered_bases": covered,
            "uncovered_bases": length - covered,
            "fraction_covered": fraction,
            "mean_depth": (st["total_depth"] / length) if length else 0.0,
            "min_depth": st["min_depth"] if st["min_depth"] is not None else 0,
            "gaps": [{"start": a, "end": b} for a, b in st["low_runs"]],
            "fully_covered": fraction >= FULL_COVERAGE_FRACTION,
        }))
    return out, None


def summarise(measured, min_depth):
    """Run-level counts, so the reader gets the verdict before the table."""
    total = len(measured)
    full = sum(1 for r in measured if r["fully_covered"])
    absent = sum(1 for r in measured if r["fraction_covered"] == 0)
    bases = sum(r["length"] for r in measured)
    uncovered = sum(r["uncovered_bases"] for r in measured)
    return {
        "regions": total,
        "regions_fully_covered": full,
        "regions_partially_covered": total - full - absent,
        "regions_absent": absent,
        "bases": bases,
        "bases_below_threshold": uncovered,
        "fraction_of_panel_covered": (bases - uncovered) / bases if bases else 0,
        "min_depth": min_depth,
        "all_covered": full == total,
    }


def _bar(fraction):
    """A pure-CSS coverage bar; no scripts, so the file opens anywhere."""
    pct = max(0.0, min(1.0, fraction)) * 100
    tone = "ok" if fraction >= 1 else ("warn" if fraction >= 0.95 else "bad")
    return (f'<div class="bar"><span class="{tone}" '
            f'style="width:{pct:.4f}%"></span></div>')


def write_html(path, sample, bam, bed, measured, summary, params):
    """
    A self-contained report: no CDN, no fonts to fetch, no JavaScript.

    It is opened from the webapp in a new tab and may be attached to an
    email or archived with the run, so it has to render years from now on
    a machine with no network.
    """
    s = summary
    verdict_class = "ok" if s["all_covered"] else (
        "warn" if s["regions_absent"] == 0 else "bad")
    if s["all_covered"]:
        verdict = (f"All {s['regions']} regions reached "
                   f"{s['min_depth']}x across every base.")
        meaning = ("A variant absent from the VCF in these regions is a "
                   "genuine negative as far as depth is concerned.")
    else:
        verdict = (f"{s['regions'] - s['regions_fully_covered']} of "
                   f"{s['regions']} regions did not reach {s['min_depth']}x "
                   f"across every base.")
        meaning = ("In the stretches listed below, a variant would not have "
                   "been called even if it were present, and would not "
                   "appear in the clinical report. Absence of a finding "
                   "there is not evidence of absence.")

    rows = []
    for r in sorted(measured, key=lambda x: (x["fraction_covered"],
                                             -x["length"])):
        gaps = ", ".join(f"{g['start']:,}&ndash;{g['end']:,}"
                         for g in r["gaps"][:12])
        if len(r["gaps"]) > 12:
            gaps += f" &hellip; (+{len(r['gaps']) - 12} more)"
        rows.append(
            "<tr class='{cls}'><td class='mono'>{chrom}:{start:,}-{end:,}</td>"
            "<td>{name}</td><td class='num'>{length:,}</td>"
            "<td class='num'>{mean:.1f}</td><td class='num'>{mind}</td>"
            "<td class='num'>{pct:.2f}%</td><td>{bar}</td>"
            "<td class='mono small'>{gaps}</td></tr>".format(
                cls="" if r["fully_covered"] else (
                    "bad" if r["fraction_covered"] == 0 else "warn"),
                chrom=html.escape(r["chrom"]), start=r["start"], end=r["end"],
                name=html.escape(r["name"] or "&ndash;"), length=r["length"],
                mean=r["mean_depth"], mind=r["min_depth"],
                pct=r["fraction_covered"] * 100, bar=_bar(r["fraction_covered"]),
                gaps=gaps or "&ndash;"))

    doc = f"""<!DOCTYPE html>
<meta charset="utf-8">
<title>Target coverage - {html.escape(sample)}</title>
<style>
  :root {{ --ink:#1c1c1e; --muted:#6b6b70; --line:#e3e3e8; --ok:#1f7a4d;
           --warn:#9a6a00; --bad:#b3261e; --bg:#fbfbfd; }}
  body {{ margin:0; padding:2.5rem 2rem 4rem; background:var(--bg);
          color:var(--ink); font:14px/1.55 -apple-system, "Segoe UI",
          Roboto, Helvetica, Arial, sans-serif; }}
  main {{ max-width:1100px; margin:0 auto; }}
  h1 {{ font-size:1.5rem; margin:0 0 .25rem; }}
  h2 {{ font-size:1.05rem; margin:2.25rem 0 .75rem; }}
  .sub {{ color:var(--muted); margin:0 0 1.75rem; }}
  .verdict {{ border-left:4px solid var(--ok); background:#fff;
              padding:1rem 1.25rem; border-radius:0 8px 8px 0;
              box-shadow:0 1px 2px rgba(0,0,0,.05); }}
  .verdict.warn {{ border-color:var(--warn); }}
  .verdict.bad {{ border-color:var(--bad); }}
  .verdict strong {{ display:block; font-size:1.05rem; margin-bottom:.35rem; }}
  .verdict p {{ margin:.35rem 0 0; color:var(--muted); }}
  dl.kv {{ display:grid; grid-template-columns:max-content 1fr; gap:.3rem 1.25rem;
           margin:0; }}
  dl.kv dt {{ color:var(--muted); }}
  dl.kv dd {{ margin:0; }}
  table {{ border-collapse:collapse; width:100%; background:#fff;
           box-shadow:0 1px 2px rgba(0,0,0,.05); border-radius:8px;
           overflow:hidden; }}
  th, td {{ padding:.5rem .7rem; border-bottom:1px solid var(--line);
            text-align:left; vertical-align:top; }}
  th {{ background:#f4f4f7; font-weight:600; font-size:.82rem;
        text-transform:uppercase; letter-spacing:.03em; color:var(--muted); }}
  td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
  tr.warn td:first-child {{ box-shadow:inset 3px 0 var(--warn); }}
  tr.bad td:first-child {{ box-shadow:inset 3px 0 var(--bad); }}
  .mono {{ font-family:ui-monospace, SFMono-Regular, Menlo, monospace;
           font-size:.85rem; }}
  .small {{ font-size:.8rem; color:var(--muted); }}
  .bar {{ width:110px; height:8px; background:#ececf1; border-radius:4px;
          overflow:hidden; }}
  .bar span {{ display:block; height:100%; }}
  .bar .ok {{ background:var(--ok); }} .bar .warn {{ background:var(--warn); }}
  .bar .bad {{ background:var(--bad); }}
  footer {{ margin-top:2.5rem; color:var(--muted); font-size:.82rem; }}
  code {{ background:#f0f0f4; padding:.1rem .3rem; border-radius:3px; }}
</style>
<main>
  <h1>Target coverage &mdash; {html.escape(sample)}</h1>
  <p class="sub">Whether the panel was sequenced deeply enough for an
    absent variant to mean anything.</p>

  <div class="verdict {verdict_class}">
    <strong>{verdict}</strong>
    <p>{meaning}</p>
  </div>

  <h2>Run</h2>
  <dl class="kv">
    <dt>Sample</dt><dd>{html.escape(sample)}</dd>
    <dt>Alignments</dt><dd class="mono">{html.escape(bam)}</dd>
    <dt>Target regions</dt><dd class="mono">{html.escape(bed)}</dd>
    <dt>Depth threshold</dt><dd>{s['min_depth']}x
      &mdash; the pipeline's own <code>--min-depth</code>, below which a
      call is filtered out</dd>
    <dt>Read filters</dt><dd>MAPQ &ge; {params['min_mapq']},
      base quality &ge; {params['min_baseq']}</dd>
    <dt>Panel</dt><dd>{s['bases']:,} bases in {s['regions']:,} regions;
      {s['bases_below_threshold']:,} below threshold
      ({s['fraction_of_panel_covered'] * 100:.2f}% covered)</dd>
    <dt>Regions</dt><dd>{s['regions_fully_covered']:,} complete,
      {s['regions_partially_covered']:,} partial,
      {s['regions_absent']:,} with no usable coverage at all</dd>
    <dt>Generated</dt><dd>{params['generated']}</dd>
  </dl>

  <h2>Regions</h2>
  <p class="sub">Worst first. Gaps are 1-based inclusive coordinates of the
    stretches below {s['min_depth']}x.</p>
  <table>
    <thead><tr><th>Region</th><th>Name</th><th>Length</th><th>Mean depth</th>
      <th>Min depth</th><th>Covered</th><th></th><th>Gaps below threshold</th>
    </tr></thead>
    <tbody>
      {"".join(rows)}
    </tbody>
  </table>

  <footer>
    <p>Depth measured with <code>{html.escape(' '.join(params['command']))}</code>.</p>
    <p>This is a coverage statement, not a variant report: it says where a
      negative can be trusted, and says nothing about the variants that
      were called.</p>
  </footer>
</main>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


def write_tsv(path, measured):
    """Per-region table, for the reader who wants it in a spreadsheet."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("chrom\tstart\tend\tname\tlength\tmean_depth\tmin_depth\t"
                 "covered_bases\tuncovered_bases\tfraction_covered\t"
                 "fully_covered\tgaps\n")
        for r in measured:
            gaps = ";".join(f"{g['start']}-{g['end']}" for g in r["gaps"])
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\t"
                     f"{r['length']}\t{r['mean_depth']:.2f}\t{r['min_depth']}\t"
                     f"{r['covered_bases']}\t{r['uncovered_bases']}\t"
                     f"{r['fraction_covered']:.6f}\t"
                     f"{int(r['fully_covered'])}\t{gaps}\n")


def run_coverage_check(bam, bed, sample_id, output_dir, min_depth,
                       min_mapq=20, min_baseq=20, dry_run=False):
    """
    Produce the coverage report for one sample.

    Returns a result dict. Every failure is reported in it rather than
    raised: a missing coverage report must never sink a run that has
    already produced its calls -- it is a caveat on the result, not the
    result.
    """
    result = {"sample": sample_id, "bam": bam, "bed": bed, "ok": False,
              "html": None, "json": None, "tsv": None, "summary": None,
              "error": None}

    if shutil.which("samtools") is None:
        result["error"] = "samtools not found on PATH"
        return result
    if not dry_run and not os.path.isfile(bam):
        result["error"] = f"alignments not found: {bam}"
        return result
    if not os.path.isfile(bed):
        result["error"] = f"BED file not found: {bed}"
        return result

    try:
        regions = parse_bed(bed)
    except (OSError, ValueError) as exc:
        result["error"] = str(exc)
        return result

    measured, error = measure(bam, bed, regions, min_depth, min_mapq,
                              min_baseq, dry_run=dry_run)
    if dry_run:
        result["error"] = "dry run: no coverage measured"
        return result
    if error:
        result["error"] = error
        return result

    summary = summarise(measured, min_depth)
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.join(output_dir, f"{sample_id}.coverage")
    params = {
        "min_mapq": min_mapq,
        "min_baseq": min_baseq,
        "command": depth_command(bam, bed, min_mapq, min_baseq),
        "generated": datetime.now(timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"),
    }
    try:
        write_html(base + ".html", sample_id, bam, bed, measured, summary,
                   params)
        write_tsv(base + ".tsv", measured)
        with open(base + ".json", "w", encoding="utf-8") as fh:
            json.dump({"sample": sample_id, "bam": bam, "bed": bed,
                       "params": params, "summary": summary,
                       "regions": measured}, fh, indent=2)
    except OSError as exc:
        result["error"] = f"could not write the report: {exc}"
        return result

    result.update({"ok": True, "html": base + ".html", "json": base + ".json",
                   "tsv": base + ".tsv", "summary": summary})
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Report how much of a BED target the sequencing covered.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--bam", required=True,
                        help="Analysis-ready BAM (post-BQSR).")
    parser.add_argument("--bed", required=True,
                        help="Target regions to check.")
    parser.add_argument("--sample-id", required=True,
                        help="Sample name; used for the output filenames.")
    parser.add_argument("--output-dir", required=True,
                        help="Where the HTML/JSON/TSV reports are written.")
    parser.add_argument("--min-depth", type=int, default=10,
                        help="Depth a base must reach to count as covered "
                             "(default: 10, matching the pipeline's "
                             "--min-depth).")
    parser.add_argument("--min-mapq", type=int, default=20,
                        help="Ignore reads below this mapping quality "
                             "(default: 20).")
    parser.add_argument("--min-baseq", type=int, default=20,
                        help="Ignore bases below this base quality "
                             "(default: 20).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the samtools command without running it.")
    args = parser.parse_args()

    result = run_coverage_check(
        args.bam, args.bed, args.sample_id, args.output_dir,
        args.min_depth, args.min_mapq, args.min_baseq, dry_run=args.dry_run)

    if not result["ok"]:
        print(f"[WARN] Coverage report not produced: {result['error']}")
        return 1
    s = result["summary"]
    print(f"[OK] Coverage report: {result['html']}")
    print(f"[OK] {s['regions_fully_covered']}/{s['regions']} regions fully "
          f"covered at {s['min_depth']}x; "
          f"{s['bases_below_threshold']:,} bases below threshold.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
