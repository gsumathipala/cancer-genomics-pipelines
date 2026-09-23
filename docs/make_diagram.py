#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
Draw the pipeline illustration: every step, with PCGR decomposed.

Why this is a script and not a drawing. The diagram has to agree with the
code, and a diagram that was drawn once by hand stops agreeing the first
time a step is renamed. The step lists below are transcribed from the
banners the engines actually print -- `# STEP n/14:` in
comprehensive_variant_calling.py and banner(n, ...) in fusion_calling.py
-- so when a step changes, this file is the one place to change, and
`--check` will tell you it needs doing.

Two renderers share ONE layout model, so the SVG and the PNG cannot drift
apart. SVG is the master (it scales, and the text stays selectable); the
PNG exists because it renders everywhere without argument.

    python3 docs/make_diagram.py            # write both files
    python3 docs/make_diagram.py --check    # verify against the code
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

W = 2040
H = 2190          # provisional; build() trims it to the content

# -- palette -------------------------------------------------------------
# Deliberately light: this is read on GitHub, printed, and projected. The
# lane colours carry the meaning (which assay a step belongs to), so they
# stay distinguishable in greyscale by differing in lightness, not only
# in hue.
INK        = "#16212e"
MUTED      = "#5a6b7d"
PAGE       = "#ffffff"
FRAME      = "#c9d4de"

DNA_FILL,  DNA_EDGE,  DNA_DARK  = "#e3effa", "#7ba7cf", "#1d4e7c"
RNA_FILL,  RNA_EDGE,  RNA_DARK  = "#e6f5ec", "#7fc09a", "#1d6b42"
PCGR_FILL, PCGR_EDGE, PCGR_DARK = "#f3ebfa", "#b193cf", "#553075"
WARN_FILL, WARN_EDGE, WARN_DARK = "#fdf0e3", "#e0b184", "#8a4f13"
OUT_FILL,  OUT_EDGE,  OUT_DARK  = "#eef1f4", "#a8b6c4", "#33465a"

# =========================================================================
# THE CONTENT -- transcribed from the engines' own step banners
# =========================================================================

DNA_STEPS = [
    ("1",  "QC & cleaning",                "fastp",
     "adapters, quality trim, HTML report"),
    ("2",  "Alignment",                    "BWA-MEM2",
     "FM-index; read groups written here"),
    ("3",  "Duplicate marking",            "GATK MarkDuplicatesSpark",
     "PCR/optical duplicates flagged, not removed"),
    ("4",  "Base quality recalibration",   "GATK BaseRecalibrator + ApplyBQSR",
     "needs dbSNP + known indels"),
    ("5",  "Somatic variant calling",      "GATK Mutect2",
     "tumour-only or matched normal; +PoN, germline resource"),
    ("6",  "Strand-bias modelling",        "GATK LearnReadOrientationModel",
     "FFPE / OxoG artefact prior for step 9"),
    ("7",  "Contamination estimation",     "GATK GetPileupSummaries + CalculateContamination",
     "cross-sample contamination fraction"),
    ("8",  "Microsatellite instability",   "MSIsensor2",
     "MSI-H / MSS from read-length distributions"),
    ("9",  "Filtering",                    "GATK FilterMutectCalls",
     "consumes steps 6 and 7; sets FILTER, keeps rows"),
    ("10", "Known-mutation annotation",    "bcftools annotate  (COSMIC)",
     "optional; licensed data, never redistributed"),
    ("11", "Gene / effect annotation",     "SnpEff",
     "transcript consequence, HGVS"),
    ("12", "Target coverage check",        "samtools depth -a",
     "WHICH REGIONS A NEGATIVE MAY SPEAK FOR"),
    ("13", "Clinical interpretation",      "PCGR      (decomposed at right)",
     "runs in its own conda environment"),
    ("14", "Summary statistics",           "run manifest",
     "counts, versions, every path written"),
]

RNA_STEPS = [
    ("1", "QC & trimming",               "fastp",
     "same cleaning as the DNA arm"),
    ("2", "Splice-aware alignment",      "STAR  --chimOutType WithinBAM Junctions",
     "chimeric detection ON; --sjdbOverhang = read length - 1"),
    ("3", "Coordinate-sorted BAM",       "samtools sort + index",
     "for QC and IGV review"),
    ("4", "Fusion calling",              "Arriba",
     "reads the chimeric alignments STAR wrote"),
    ("5", "RNA library QC",              "rna_qc_report.py",
     "IS A NEGATIVE INTERPRETABLE? depth, unique%, chimeric"),
    ("6", "Fusion report",               "fusion_report.py",
     "ranked HTML + JSON + PCGR fusion TSV"),
]

PCGR_STAGES = [
    ("1",  "Environment & version check",
     "refuses a mismatched reference bundle up front"),
    ("2",  "Input validation",
     "VCF and/or fusion TSV: FusionGene, LeftBreakpoint, RightBreakpoint, SplitReads"),
    ("3",  "VEP annotation",
     "the heavy pass - ~24 GB cache, pinned version, canonical transcript flagged"),
    ("4",  "vcfanno",
     "overlays ClinVar, gnomAD, dbNSFP, CGI, CiVIC, Mutational hotspots"),
    ("5",  "pcgr-summarise",
     "cancer-gene annotation, ONCOGENICITY SCORING, hotspot matching"),
    ("6",  "vcf2maf",
     "MAF conversion for downstream compatibility"),
    ("7",  "CNA / RNA-fusion / RNA-expression sections",
     "each annotated against the biomarker databases, or announced as omitted"),
    ("8",  "Two-hit analysis",
     "tumour suppressors with a somatic variant AND loss of the other allele"),
    ("9",  "OncoKB annotation",
     "only with an API token; otherwise explicitly skipped, NA columns added"),
    ("10", "Report generation",
     "Quarto / R interactive HTML, plus TSV and XLSX"),
]

ONCOGENICITY = [
    (">= 9",      "Oncogenic"),
    ("5 to 8",    "Likely oncogenic"),
    ("-1 to 4",   "Variant of uncertain significance"),
    ("-6 to -2",  "Likely benign"),
    ("<= -7",     "Benign"),
]

# The silent-failure guards. These are the reason the pipeline is shaped
# the way it is, so they belong on the diagram rather than in a footnote.
GUARDS = [
    ("Chimeric output verified",
     "STAR's own chimeric count is compared against the ch:A:1 tags in the "
     "BAM. They can disagree - 116 reported, 0 written - and the fusion "
     "caller reads only the BAM, so it exits cleanly with an empty table "
     "that looks exactly like a true negative."),
    ("Coverage before conclusions",
     "Step 12 states which target regions reached the depth threshold. "
     "Without it, 'no mutation detected' silently covers regions that were "
     "never sequenced deeply enough to speak for."),
    ("TMB denominator measured, not assumed",
     "PCGR defaults to a 34 Mb exome. A 1.9 Mb panel scored against 34 Mb "
     "reads ~18x too low. The footprint is measured from the merged BED "
     "and passed as --effective_target_size_mb."),
    ("An empty fusion file is never written",
     "PCGR rejects one. 'No fusions found' must therefore produce NO file, "
     "so the absence is explicit rather than an empty table."),
    ("Contig naming harmonised",
     "chr1/chrM (UCSC) against 1/MT (Ensembl). Mismatched contigs produce "
     "zero overlaps and no error at all."),
]

OUTPUTS = [
    ("PCGR clinical report",   "HTML + TSV + XLSX, tiered by tumour site"),
    ("Target coverage report", "per-region depth; what a negative covers"),
    ("Fusion report",          "ranked HTML + JSON (RNA arm)"),
    ("RNA library QC",         "whether a negative is interpretable"),
    ("Filtered / annotated VCF", "the call set itself"),
    ("Run manifest (JSON)",    "versions, parameters, every path written"),
]


# =========================================================================
# LAYOUT MODEL -- both renderers consume this and nothing else
# =========================================================================

class Scene:
    def __init__(self):
        self.items = []

    def rect(self, x, y, w, h, fill, edge, r=8, width=1.6, dash=None):
        self.items.append(("rect", dict(x=x, y=y, w=w, h=h, fill=fill,
                                        edge=edge, r=r, width=width,
                                        dash=dash)))

    def text(self, x, y, s, size=15, fill=INK, anchor="start",
             bold=False, mono=False, bound=None):
        # `bound` is the x this text must not cross -- the right edge of
        # the box it belongs to. Set it for anything inside a lane and
        # --check will measure it; leave it None for page-wide banners.
        self.items.append(("text", dict(x=x, y=y, s=s, size=size, fill=fill,
                                        anchor=anchor, bold=bold, mono=mono,
                                        bound=bound)))

    def line(self, x1, y1, x2, y2, stroke=MUTED, width=1.6, dash=None):
        self.items.append(("line", dict(x1=x1, y1=y1, x2=x2, y2=y2,
                                        stroke=stroke, width=width,
                                        dash=dash)))

    def poly(self, points, fill):
        self.items.append(("poly", dict(points=points, fill=fill)))

    def arrow(self, x1, y1, x2, y2, stroke=MUTED, width=1.8, head=8.0):
        """
        A line plus an explicit arrowhead.

        The head is a polygon in the model rather than an SVG marker,
        because a marker has no equivalent in the PNG renderer and the two
        would then disagree about which way an arrow points.
        """
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy) or 1.0
        ux, uy = dx / length, dy / length
        bx, by = x2 - ux * head, y2 - uy * head
        self.line(x1, y1, bx, by, stroke, width)
        px, py = -uy, ux
        self.poly([(x2, y2),
                   (bx + px * head * 0.55, by + py * head * 0.55),
                   (bx - px * head * 0.55, by - py * head * 0.55)], stroke)


def wrap(text, limit):
    """Greedy wrap by character count -- close enough at these sizes."""
    words, lines, cur = text.split(), [], ""
    for word in words:
        trial = (cur + " " + word).strip()
        if len(trial) <= limit:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def build():
    s = Scene()

    # ---- title ---------------------------------------------------------
    s.text(60, 62, "Cancer Genomics Pipelines", size=40, bold=True)
    s.text(60, 98,
           "Complete step decomposition - DNA and RNA arms, with PCGR "
           "broken out stage by stage", size=19, fill=MUTED)
    s.rect(60, 116, 300, 30, WARN_FILL, WARN_EDGE, r=6)
    s.text(74, 137, "RESEARCH USE ONLY", size=15, bold=True, fill=WARN_DARK)
    s.text(380, 137, "Not validated for clinical or diagnostic use.",
           size=15, fill=MUTED)
    s.text(W - 60, 62, "Created by Brainstorm, 2026", size=15,
           fill=MUTED, anchor="end")
    s.text(W - 60, 88, "github.com/gsumathipala/cancer-genomics-pipelines",
           size=13, fill=MUTED, anchor="end", mono=True)
    s.line(60, 162, W - 60, 162, FRAME, 1.4)

    # ---- inputs and the pathway choice ---------------------------------
    y = 186
    s.rect(60, y, 900, 74, OUT_FILL, OUT_EDGE)
    s.text(78, y + 28, "INPUT", size=13, bold=True, fill=OUT_DARK)
    s.text(78, y + 52,
           "Paired FASTQ  -  tumour, optionally a matched normal  |  "
           "reference genome  |  panel profile (one name configures the run)",
           size=15)

    s.rect(1000, y, 980, 74, WARN_FILL, WARN_EDGE)
    s.text(1018, y + 28, "PATHWAY CHOSEN FIRST", size=13, bold=True,
           fill=WARN_DARK)
    s.text(1018, y + 52,
           "DNA  |  RNA  |  Hybrid (both arms, one specimen)   - the web "
           "interface asks before anything else",
           size=15)

    # ---- lane headers --------------------------------------------------
    lane_y = 300
    dna_x, dna_w = 60, 620
    rna_x, rna_w = 706, 436
    pcgr_x, pcgr_w = 1252, 728

    s.rect(dna_x, lane_y, dna_w, 46, DNA_DARK, DNA_DARK, r=8, width=0)
    s.text(dna_x + 18, lane_y + 30, "DNA ARM", size=19, bold=True,
           fill="#ffffff")
    s.text(dna_x + dna_w - 18, lane_y + 30,
           "comprehensive_variant_calling.py   -   env: cancer_pipeline",
           size=12, fill="#cfe2f3", anchor="end", mono=True)

    s.rect(rna_x, lane_y, rna_w, 46, RNA_DARK, RNA_DARK, r=8, width=0)
    s.text(rna_x + 18, lane_y + 30, "RNA ARM", size=19, bold=True,
           fill="#ffffff")
    s.text(rna_x + rna_w - 18, lane_y + 30,
           "fusion_calling.py   -   env: cancer_rna",
           size=12, fill="#d3efe0", anchor="end", mono=True)

    s.rect(pcgr_x, lane_y, pcgr_w, 46, PCGR_DARK, PCGR_DARK, r=8, width=0)
    s.text(pcgr_x + 18, lane_y + 30, "PCGR, DECOMPOSED", size=19,
           bold=True, fill="#ffffff")
    s.text(pcgr_x + pcgr_w - 18, lane_y + 30,
           "pcgr_report.py   -   env: pcgr / pcgrr",
           size=12, fill="#e2d5f0", anchor="end", mono=True)

    # ---- DNA steps -----------------------------------------------------
    step_h, gap = 66, 12
    top = lane_y + 62
    dna_bottom = top
    for i, (num, title, tool, note) in enumerate(DNA_STEPS):
        by = top + i * (step_h + gap)
        dna_bottom = by + step_h
        s.rect(dna_x, by, dna_w, step_h, DNA_FILL, DNA_EDGE)
        s.rect(dna_x, by, 44, step_h, DNA_EDGE, DNA_EDGE, r=8, width=0)
        s.text(dna_x + 22, by + 41, num, size=19, bold=True,
               fill="#ffffff", anchor="middle")
        edge = dna_x + dna_w - 10
        s.text(dna_x + 58, by + 26, title, size=16, bold=True, bound=edge)
        s.text(dna_x + 58, by + 45, tool, size=13, fill=DNA_DARK, mono=True,
               bound=edge)
        s.text(dna_x + 58, by + 61, note, size=12, bound=edge,
               fill=WARN_DARK if note.isupper() else MUTED)
        if i < len(DNA_STEPS) - 1:
            s.arrow(dna_x + 22, by + step_h, dna_x + 22, by + step_h + gap,
                    DNA_EDGE, 2.0)

    # ---- RNA steps -----------------------------------------------------
    rna_bottom = top
    for i, (num, title, tool, note) in enumerate(RNA_STEPS):
        by = top + i * (step_h + gap)
        rna_bottom = by + step_h
        s.rect(rna_x, by, rna_w, step_h, RNA_FILL, RNA_EDGE)
        s.rect(rna_x, by, 44, step_h, RNA_EDGE, RNA_EDGE, r=8, width=0)
        s.text(rna_x + 22, by + 41, num, size=19, bold=True,
               fill="#ffffff", anchor="middle")
        edge = rna_x + rna_w - 10
        s.text(rna_x + 58, by + 26, title, size=16, bold=True, bound=edge)
        s.text(rna_x + 58, by + 45, tool, size=12, fill=RNA_DARK, mono=True,
               bound=edge)
        s.text(rna_x + 58, by + 61, note, size=12, bound=edge,
               fill=WARN_DARK if "?" in note else MUTED)
        if i < len(RNA_STEPS) - 1:
            s.arrow(rna_x + 22, by + step_h, rna_x + 22, by + step_h + gap,
                    RNA_EDGE, 2.0)

    # ---- guards, under the RNA lane ------------------------------------
    gy = rna_bottom + 34
    s.text(rna_x, gy, "SILENT-FAILURE GUARDS", size=15, bold=True,
           fill=WARN_DARK)
    s.text(rna_x, gy + 20,
           "Each catches a failure that yields a plausible, wrong report",
           size=12, fill=MUTED, bound=rna_x + rna_w)
    gy += 36
    for name, why in GUARDS:
        lines = wrap(why, 62)
        box_h = 26 + 15 * len(lines)
        s.rect(rna_x, gy, rna_w, box_h, WARN_FILL, WARN_EDGE, r=6, width=1.3)
        s.text(rna_x + 14, gy + 20, name, size=14, bold=True,
               fill=WARN_DARK, bound=rna_x + rna_w - 10)
        for j, line in enumerate(lines):
            s.text(rna_x + 14, gy + 38 + j * 15, line, size=11.5, fill=INK,
                   bound=rna_x + rna_w - 10)
        gy += box_h + 10

    # ---- PCGR stages ---------------------------------------------------
    py = top
    st_h = 58
    for num, title, note in PCGR_STAGES:
        s.rect(pcgr_x, py, pcgr_w, st_h, PCGR_FILL, PCGR_EDGE)
        s.rect(pcgr_x, py, 40, st_h, PCGR_EDGE, PCGR_EDGE, r=8, width=0)
        s.text(pcgr_x + 20, py + 36, num, size=16, bold=True,
               fill="#ffffff", anchor="middle")
        edge = pcgr_x + pcgr_w - 10
        s.text(pcgr_x + 54, py + 25, title, size=16, bold=True, bound=edge)
        s.text(pcgr_x + 54, py + 45, note, size=12, fill=MUTED, bound=edge)
        py += st_h + 8

    # ---- oncogenicity banding -----------------------------------------
    py += 22
    s.text(pcgr_x, py, "ONCOGENICITY SCORING  (stage 5)", size=15,
           bold=True, fill=PCGR_DARK)
    s.text(pcgr_x, py + 20,
           "ClinGen / CGC / VICC points system. Criteria score toward the "
           "oncogenic or benign pole", size=12, fill=MUTED)
    s.text(pcgr_x, py + 36,
           "(e.g. ONCG_OVS1_A = +8, a null variant in a tumour suppressor; "
           "ONCG_SBS2_A = -4). The total is banded:", size=12, fill=MUTED)
    py += 52
    row_h = 30
    for i, (score, label) in enumerate(ONCOGENICITY):
        ry = py + i * row_h
        s.rect(pcgr_x, ry, 150, row_h, "#ffffff", PCGR_EDGE, r=4, width=1.1)
        s.rect(pcgr_x + 150, ry, pcgr_w - 150, row_h, PCGR_FILL, PCGR_EDGE,
               r=4, width=1.1)
        s.text(pcgr_x + 75, ry + 20, score, size=13, anchor="middle",
               mono=True, bold=True)
        s.text(pcgr_x + 164, ry + 20, label, size=14)
    py += len(ONCOGENICITY) * row_h + 26

    # ---- the two settings that decide whether the report means anything
    s.rect(pcgr_x, py, pcgr_w, 108, WARN_FILL, WARN_EDGE, r=6)
    s.text(pcgr_x + 16, py + 24,
           "ACTIONABILITY TIERS ARE TUMOUR-SITE SPECIFIC", size=14,
           bold=True, fill=WARN_DARK)
    for j, line in enumerate(wrap(
            "AMP/ASCO/CAP tiers are assigned from biomarker evidence FOR "
            "THE SPECIFIED TUMOUR SITE. The same variant can be predictive "
            "in one tissue and merely oncogenic in another - which is why "
            "--tumor_site is a per-sample field on the worksheet, not a "
            "run-wide default.", 88)):
        s.text(pcgr_x + 16, py + 46 + j * 16, line, size=12, fill=INK)
    pcgr_bottom = py + 108

    # ---- both arms converge on PCGR ------------------------------------
    # Drawn rather than implied: the single most common misreading of this
    # design is that the RNA arm has a reporting path of its own. It does
    # not -- PCGR 2.3.2 takes a fusion TSV as a first-class input, so one
    # report covers whichever arms were run.
    # Both feeds enter at PCGR stage 2 (input validation), which is where
    # they genuinely arrive -- a VCF and a fusion TSV are validated by the
    # same stage before anything is annotated.
    channel_x = rna_x + rna_w + 62
    stage2_y = top + (st_h + 8) + st_h / 2

    dna13_y = top + 12 * (step_h + gap) + step_h / 2
    s.line(dna_x + dna_w, dna13_y, channel_x, dna13_y, PCGR_EDGE, 2.2)
    s.text((rna_x + rna_w + channel_x) / 2, dna13_y - 9, "filtered VCF",
           size=11, fill=PCGR_DARK, anchor="middle")

    rna6_y = top + 5 * (step_h + gap) + step_h / 2
    s.line(rna_x + rna_w, rna6_y, channel_x, rna6_y, PCGR_EDGE, 2.2)
    s.text((rna_x + rna_w + channel_x) / 2, rna6_y - 9, "fusion TSV",
           size=11, fill=PCGR_DARK, anchor="middle")

    s.line(channel_x, dna13_y, channel_x, stage2_y, PCGR_EDGE, 2.2)
    s.arrow(channel_x, stage2_y, pcgr_x - 3, stage2_y, PCGR_EDGE, 2.2)

    # ---- outputs -------------------------------------------------------
    oy = max(dna_bottom, gy, pcgr_bottom) + 46
    s.line(60, oy - 24, W - 60, oy - 24, FRAME, 1.4)
    s.text(60, oy, "WHAT A RUN LEAVES BEHIND", size=17, bold=True,
           fill=OUT_DARK)
    s.text(60, oy + 20,
           "Every path below is recorded in the run manifest, and every "
           "one is a hot link on the run's page in the web interface",
           size=13, fill=MUTED)
    oy += 36
    col_w = (W - 120 - 2 * 16) / 3
    for i, (name, what) in enumerate(OUTPUTS):
        cx = 60 + (i % 3) * (col_w + 16)
        cy = oy + (i // 3) * 76
        s.rect(cx, cy, col_w, 64, OUT_FILL, OUT_EDGE, r=6, width=1.3)
        s.text(cx + 14, cy + 26, name, size=15, bold=True, fill=OUT_DARK)
        s.text(cx + 14, cy + 47, what, size=12, fill=MUTED)
    bottom = oy + 2 * 76 - 12

    # ---- footer --------------------------------------------------------
    bottom += 30
    s.line(60, bottom, W - 60, bottom, FRAME, 1.4)
    s.text(60, bottom + 26,
           "Generated by docs/make_diagram.py from the step banners the "
           "engines print. Run it with --check to verify the diagram still "
           "agrees with the code.", size=12, fill=MUTED)
    s.text(60, bottom + 46,
           "Full prose dissection: PIPELINE_ANATOMY.md   |   "
           "Panel profiles: PANELS.md   |   Testing: TESTING.md",
           size=12, fill=MUTED)
    bottom += 70

    # The background goes in FIRST, now that the height is known.
    s.items.insert(0, ("rect", dict(x=0, y=0, w=W, h=bottom, fill=PAGE,
                                    edge=PAGE, r=0, width=0, dash=None)))
    return s, bottom


# =========================================================================
# RENDERER 1: SVG
# =========================================================================

def escape(t):
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def to_svg(scene, path):
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" font-family="DejaVu Sans, Verdana, sans-serif">',
        '<!-- Created by Brainstorm, 2026. Generated by docs/make_diagram.py'
        ' - do not edit by hand. -->',
    ]
    for kind, it in scene.items:
        if kind == "rect":
            dash = f' stroke-dasharray="{it["dash"]}"' if it["dash"] else ""
            out.append(
                f'<rect x="{it["x"]}" y="{it["y"]}" width="{it["w"]}" '
                f'height="{it["h"]}" rx="{it["r"]}" fill="{it["fill"]}" '
                f'stroke="{it["edge"]}" stroke-width="{it["width"]}"{dash}/>')
        elif kind == "line":
            dash = f' stroke-dasharray="{it["dash"]}"' if it["dash"] else ""
            out.append(
                f'<line x1="{it["x1"]}" y1="{it["y1"]}" x2="{it["x2"]}" '
                f'y2="{it["y2"]}" stroke="{it["stroke"]}" '
                f'stroke-width="{it["width"]}"{dash}/>')
        elif kind == "poly":
            pts = " ".join(f"{x},{y}" for x, y in it["points"])
            out.append(f'<polygon points="{pts}" fill="{it["fill"]}"/>')
        else:
            anchor = {"start": "start", "middle": "middle",
                      "end": "end"}[it["anchor"]]
            family = ('font-family="DejaVu Sans Mono, monospace" '
                      if it["mono"] else "")
            weight = 'font-weight="bold" ' if it["bold"] else ""
            out.append(
                f'<text x="{it["x"]}" y="{it["y"]}" font-size="{it["size"]}" '
                f'fill="{it["fill"]}" text-anchor="{anchor}" '
                f'{family}{weight}>{escape(it["s"])}</text>')
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


# =========================================================================
# RENDERER 2: PNG  (Pillow -- so the image renders anywhere, no questions)
# =========================================================================

FONTS = {
    (False, False): "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    (True,  False): "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    (False, True):  "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    (True,  True):  "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf",
}


def to_png(scene, path, scale=2):
    from PIL import Image, ImageDraw, ImageFont

    cache = {}

    def font(size, bold, mono):
        key = (round(size * scale), bold, mono)
        if key not in cache:
            fp = FONTS[(bold, mono)]
            if not os.path.exists(fp):
                fp = FONTS[(bold, False)]
            cache[key] = ImageFont.truetype(fp, key[0])
        return cache[key]

    img = Image.new("RGB", (W * scale, H * scale), PAGE)
    d = ImageDraw.Draw(img)

    for kind, it in scene.items:
        if kind == "rect":
            box = [it["x"] * scale, it["y"] * scale,
                   (it["x"] + it["w"]) * scale, (it["y"] + it["h"]) * scale]
            radius = max(0, int(it["r"] * scale))
            outline = it["edge"] if it["width"] else None
            d.rounded_rectangle(box, radius=radius, fill=it["fill"],
                                outline=outline,
                                width=max(1, int(it["width"] * scale)))
        elif kind == "line":
            d.line([it["x1"] * scale, it["y1"] * scale,
                    it["x2"] * scale, it["y2"] * scale],
                   fill=it["stroke"], width=max(1, int(it["width"] * scale)))
        elif kind == "poly":
            d.polygon([(x * scale, y * scale) for x, y in it["points"]],
                      fill=it["fill"])
        else:
            f = font(it["size"], it["bold"], it["mono"])
            anchor = {"start": "ls", "middle": "ms", "end": "rs"}[it["anchor"]]
            d.text((it["x"] * scale, it["y"] * scale), it["s"],
                   font=f, fill=it["fill"], anchor=anchor)

    img.save(path, "PNG", optimize=True)


# =========================================================================
# --check : does the diagram still agree with the code?
# =========================================================================

def overflows(scene):
    """
    Text drawn wider than the lane it sits in.

    Measured with the same font the PNG uses, so it reflects what is
    actually rendered rather than a character-count guess.
    """
    from PIL import ImageFont
    found, cache = [], {}

    def width(item):
        key = (round(item["size"]), item["bold"], item["mono"])
        if key not in cache:
            fp = FONTS[(item["bold"], item["mono"])]
            fp = fp if os.path.exists(fp) else FONTS[(item["bold"], False)]
            cache[key] = ImageFont.truetype(fp, key[0])
        return cache[key].getlength(item["s"])

    for kind, it in scene.items:
        if kind != "text" or not it.get("bound"):
            continue
        over = it["x"] + width(it) - it["bound"]
        if over > 0:
            found.append(f"{it['s'][:60]!r} overruns its box by "
                         f"{int(over)} px")
    return found


def check():
    """Compare the transcribed step counts against the engines."""
    problems = []

    with open(os.path.join(ROOT, "comprehensive_variant_calling.py"),
              encoding="utf-8") as fh:
        src = fh.read()
    banners = re.findall(r'# STEP (\d+)/(\d+): ([^"\']+)', src)
    if len(banners) != len(DNA_STEPS):
        problems.append(f"DNA: code prints {len(banners)} steps, "
                        f"diagram draws {len(DNA_STEPS)}")
    for (num, _total, title), (dnum, dtitle, _t, _n) in zip(banners,
                                                            DNA_STEPS):
        if num != dnum:
            problems.append(f"DNA step {num} is drawn as {dnum}")

    with open(os.path.join(ROOT, "fusion_calling.py"),
              encoding="utf-8") as fh:
        src = fh.read()
    rna = re.findall(r'banner\((\d+), ', src)
    if len(rna) != len(RNA_STEPS):
        problems.append(f"RNA: code prints {len(rna)} steps, "
                        f"diagram draws {len(RNA_STEPS)}")

    return problems


def main(argv):
    global H
    scene, used = build()
    H = int(used)
    if "--check" in argv:
        problems = check()
        for p in problems:
            print(f"[DRIFT] {p}")
        print("diagram agrees with the code" if not problems
              else f"{len(problems)} mismatch(es)")
        return 1 if problems else 0

    svg = os.path.join(HERE, "pipeline_overview.svg")
    png = os.path.join(HERE, "pipeline_overview.png")
    to_svg(scene, svg)
    to_png(scene, png)
    print(f"wrote {os.path.relpath(svg, ROOT)} and "
          f"{os.path.relpath(png, ROOT)}")
    print(f"canvas {W}x{H}, content reaches y={int(used)}")
    problems = check()
    for p in problems:
        print(f"[DRIFT] {p}")
    for o in overflows(scene):
        print(f"[OVERFLOW] {o}")
        problems.append(o)
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
