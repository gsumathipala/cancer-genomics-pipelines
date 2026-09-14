#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
pdfwriter.py
============
A small, dependency-free PDF writer.

WHY NOT reportlab / WeasyPrint / wkhtmltopdf?
--------------------------------------------
  The analysis scripts in this project deliberately use only the Python
  standard library, and a clinical report generator is a bad place to
  inherit a heavyweight rendering stack (WeasyPrint pulls in cairo/pango;
  wkhtmltopdf is an external binary that is frequently unavailable on
  cluster nodes). Everything this report needs -- headings, key/value
  blocks, simple tables, page breaks and page numbers -- is a few hundred
  lines against the PDF spec using the 14 standard fonts, which require no
  font embedding.

  The trade is that this writer is intentionally limited: left-to-right
  Latin text, no images, no unicode beyond Latin-1. That is enough for a
  variant report and keeps the output auditable.

WHAT A PDF ACTUALLY IS (enough to follow the code)
---------------------------------------------------
  A PDF file is a header, a set of numbered "objects", a cross-reference
  table listing the byte offset of every object, and a trailer pointing at
  the cross-reference table. Readers seek to the trailer first, so those
  byte offsets must be exact -- which is why write() below records the
  offset of each object as it is emitted rather than computing them
  afterwards.

  The object graph produced here:
      1  Catalog          -- document root
      2  Pages            -- the page tree
      3  Font (Helvetica)
      4  Font (Helvetica-Bold)
      5.. Page + Contents pairs, one of each per page
"""

import datetime

# --- Page geometry, in PDF points (1 pt = 1/72 inch) ---------------------
A4 = (595.28, 841.89)
LETTER = (612.0, 792.0)

# Widths of Helvetica characters, in 1/1000 em. Only what we need to wrap
# text sensibly; anything not listed falls back to AVG_WIDTH. Getting this
# roughly right matters because a too-generous estimate overflows the right
# margin, which looks broken in a clinical document.
_HELV_WIDTHS = {
    " ": 278, "!": 278, '"': 355, "#": 556, "$": 556, "%": 889, "&": 667,
    "'": 191, "(": 333, ")": 333, "*": 389, "+": 584, ",": 278, "-": 333,
    ".": 278, "/": 278, ":": 278, ";": 278, "<": 584, "=": 584, ">": 584,
    "?": 556, "@": 1015, "[": 278, "\\": 278, "]": 278, "^": 469, "_": 556,
    "`": 333, "{": 334, "|": 260, "}": 334, "~": 584,
}
for _c in "0123456789":
    _HELV_WIDTHS[_c] = 556
for _c, _w in zip("abcdefghijklmnopqrstuvwxyz",
                  (556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500,
                   222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500,
                   722, 500, 500, 500)):
    _HELV_WIDTHS[_c] = _w
for _c, _w in zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                  (667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667,
                   556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667,
                   944, 667, 667, 611)):
    _HELV_WIDTHS[_c] = _w
AVG_WIDTH = 556


def text_width(text, size, bold=False):
    """
    Approximate the rendered width of `text` in points.

    Bold Helvetica is a little wider than regular; the 1.05 factor is a
    deliberate slight over-estimate so wrapping errs towards breaking a
    line early rather than running past the margin.
    """
    total = sum(_HELV_WIDTHS.get(ch, AVG_WIDTH) for ch in text)
    width = total * size / 1000.0
    return width * 1.05 if bold else width


# WinAnsiEncoding is not Latin-1: it fills 0x80-0x9F with typographic
# punctuation that Latin-1 leaves undefined. Mapping these means an em dash
# or a curly quote pasted from a Word document renders properly instead of
# turning into "?".
_WINANSI_EXTRAS = {
    "€": 0x80, "‚": 0x82, "ƒ": 0x83, "„": 0x84,
    "…": 0x85, "†": 0x86, "‡": 0x87, "ˆ": 0x88,
    "‰": 0x89, "Š": 0x8a, "‹": 0x8b, "Œ": 0x8c,
    "Ž": 0x8e, "‘": 0x91, "’": 0x92, "“": 0x93,
    "”": 0x94, "•": 0x95, "–": 0x96, "—": 0x97,
    "˜": 0x98, "™": 0x99, "š": 0x9a, "›": 0x9b,
    "œ": 0x9c, "ž": 0x9e, "Ÿ": 0x9f,
}


def unsupported_characters(text):
    """
    Return the set of characters in `text` this writer cannot render.

    The standard-14 fonts used here cover WinAnsi only, so scripts such as
    Chinese, Arabic or Cyrillic have no glyph. Callers should surface this
    rather than let it pass: a patient name silently rendered as "???" in a
    clinical report is a data-integrity problem, not a cosmetic one.
    """
    bad = set()
    for ch in str(text):
        if ch in ("\\", "(", ")") or ord(ch) < 32:
            continue
        if ord(ch) < 256 or ch in _WINANSI_EXTRAS:
            continue
        bad.add(ch)
    return bad


def _escape(text):
    """
    Make a Python string safe inside a PDF literal string.

    Backslash and both parentheses are structural in PDF strings and must
    be escaped. Characters with no WinAnsi glyph become "?" -- see
    unsupported_characters(), which lets callers warn about that before
    a report is issued.
    """
    out = []
    for ch in str(text):
        if ch in ("\\", "(", ")"):
            out.append("\\" + ch)
        elif ord(ch) < 32:
            out.append(" ")
        elif ord(ch) < 256:
            out.append(ch)
        elif ch in _WINANSI_EXTRAS:
            out.append(chr(_WINANSI_EXTRAS[ch]))
        else:
            out.append("?")
    return "".join(out)


class PDFDocument:
    """
    Builds a paginated text document and writes it as a PDF.

    Content is added top-to-bottom; the cursor (`self._y`) moves down the
    page and a new page is started automatically when it would cross the
    bottom margin. Nothing is written to disk until save().
    """

    def __init__(self, title="Report", page_size=A4, margin=56.0,
                 footer=None):
        self.width, self.height = page_size
        self.margin = margin
        self.title = title
        self.footer = footer
        self._pages = []        # each entry is a list of content-stream ops
        self._ops = None        # ops for the page currently being written
        self._y = 0.0
        self._new_page()

    # -- page handling ---------------------------------------------------
    def _new_page(self):
        self._ops = []
        self._pages.append(self._ops)
        self._y = self.height - self.margin

    def _need(self, height):
        """Start a new page if `height` more points will not fit."""
        # Leave room for the footer line as well as the bottom margin.
        if self._y - height < self.margin + 24:
            self._new_page()

    @property
    def content_width(self):
        return self.width - 2 * self.margin

    # -- primitives ------------------------------------------------------
    def _draw_text(self, text, x, y, size, bold=False, grey=None):
        font = "/F2" if bold else "/F1"
        colour = f"{grey:.2f} {grey:.2f} {grey:.2f} rg\n" if grey is not None else ""
        reset = "0 0 0 rg\n" if grey is not None else ""
        self._ops.append(
            f"{colour}BT {font} {size} Tf 1 0 0 1 {x:.2f} {y:.2f} Tm "
            f"({_escape(text)}) Tj ET\n{reset}"
        )

    def _wrap(self, text, size, bold, width):
        """Greedy word wrap. Over-long single words are hard-split."""
        words, lines, current = str(text).split(), [], ""
        for word in words:
            trial = f"{current} {word}".strip()
            if text_width(trial, size, bold) <= width or not current:
                # A single word wider than the column still has to go
                # somewhere; break it rather than overflow the margin.
                while text_width(trial, size, bold) > width and len(trial) > 1:
                    cut = max(1, int(len(trial) * width /
                                     max(text_width(trial, size, bold), 1)))
                    lines.append(trial[:cut])
                    trial = trial[cut:]
                current = trial
            else:
                lines.append(current)
                current = word
        if current:
            lines.append(current)
        return lines or [""]

    # -- public content API ----------------------------------------------
    def spacer(self, height=10.0):
        self._need(height)
        self._y -= height

    def rule(self, grey=0.75, above=4.0, below=8.0):
        self._need(above + below + 1)
        self._y -= above
        self._ops.append(
            f"{grey:.2f} {grey:.2f} {grey:.2f} RG 0.6 w "
            f"{self.margin:.2f} {self._y:.2f} m "
            f"{self.width - self.margin:.2f} {self._y:.2f} l S\n"
        )
        self._y -= below

    def heading(self, text, level=1):
        size = {1: 17, 2: 12.5, 3: 10.5}.get(level, 11)
        gap_before = {1: 4, 2: 14, 3: 10}.get(level, 10)
        self._need(size + gap_before + 6)
        self._y -= gap_before
        self._draw_text(text, self.margin, self._y, size, bold=True)
        self._y -= size + 2
        if level == 1:
            self.rule(grey=0.4, above=2, below=6)

    def paragraph(self, text, size=9.5, grey=None, indent=0.0):
        for line in self._wrap(text, size, False, self.content_width - indent):
            self._need(size + 3)
            self._draw_text(line, self.margin + indent, self._y, size,
                            grey=grey)
            self._y -= size + 3

    def key_values(self, pairs, label_width=140.0, size=9.5):
        """
        A two-column block of label/value rows.

        Used for patient and run metadata, where the label column should
        line up and long values must wrap under themselves rather than
        pushing the layout sideways.
        """
        for label, value in pairs:
            value = "-" if value in (None, "") else str(value)
            lines = self._wrap(value, size, False,
                               self.content_width - label_width)
            self._need((size + 3) * len(lines))
            self._draw_text(f"{label}", self.margin, self._y, size,
                            bold=True, grey=0.35)
            for i, line in enumerate(lines):
                if i:
                    self._need(size + 3)
                self._draw_text(line, self.margin + label_width, self._y, size)
                self._y -= size + 3

    def table(self, headers, rows, widths=None, size=8.5, max_rows=None):
        """
        A simple grid with a bold header row and hairline separators.

        `widths` are fractions of the content width; they are normalised so
        a caller cannot accidentally run the table off the page.
        """
        ncols = len(headers)
        if not widths:
            widths = [1.0 / ncols] * ncols
        total = sum(widths) or 1.0
        cols = [self.content_width * w / total for w in widths]

        truncated = 0
        if max_rows is not None and len(rows) > max_rows:
            truncated = len(rows) - max_rows
            rows = rows[:max_rows]

        def emit_header():
            self._need(size + 8)
            x = self.margin
            for head, cw in zip(headers, cols):
                for line in self._wrap(head, size, True, cw - 6)[:1]:
                    self._draw_text(line, x, self._y, size, bold=True)
                x += cw
            self._y -= size + 3
            self._ops.append(
                f"0.75 0.75 0.75 RG 0.5 w {self.margin:.2f} {self._y:.2f} m "
                f"{self.width - self.margin:.2f} {self._y:.2f} l S\n")
            self._y -= 5

        emit_header()
        for row in rows:
            cells = [self._wrap(str(c if c not in (None, "") else "-"),
                                size, False, cw - 6)
                     for c, cw in zip(row, cols)]
            tall = max(len(c) for c in cells)
            if self._y - (size + 2) * tall < self.margin + 24:
                self._new_page()
                emit_header()
            top = self._y
            x = self.margin
            for cell, cw in zip(cells, cols):
                y = top
                for line in cell:
                    self._draw_text(line, x, y, size)
                    y -= size + 2
                x += cw
            self._y = top - (size + 2) * tall - 2

        if truncated:
            self.paragraph(f"... and {truncated} more row(s) not shown.",
                           size=8, grey=0.45)

    # -- output ----------------------------------------------------------
    def _page_objects(self, first_obj):
        """Build the (page, contents) object pairs and their byte content."""
        objects = []
        total = len(self._pages)
        for index, ops in enumerate(self._pages):
            page_num = index + 1
            body = "".join(ops)
            if self.footer:
                body += (
                    f"0.45 0.45 0.45 rg BT /F1 8 Tf 1 0 0 1 "
                    f"{self.margin:.2f} {self.margin - 18:.2f} Tm "
                    f"({_escape(self.footer)}) Tj ET 0 0 0 rg\n"
                )
            label = f"Page {page_num} of {total}"
            body += (
                f"0.45 0.45 0.45 rg BT /F1 8 Tf 1 0 0 1 "
                f"{self.width - self.margin - text_width(label, 8):.2f} "
                f"{self.margin - 18:.2f} Tm ({_escape(label)}) Tj ET 0 0 0 rg\n"
            )
            stream = body.encode("latin-1", "replace")
            content_obj = first_obj + index * 2 + 1
            page_obj = first_obj + index * 2
            objects.append((page_obj,
                            f"<< /Type /Page /Parent 2 0 R "
                            f"/MediaBox [0 0 {self.width:.2f} {self.height:.2f}] "
                            f"/Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> "
                            f"/Contents {content_obj} 0 R >>".encode("latin-1")))
            objects.append((content_obj,
                            b"<< /Length " + str(len(stream)).encode() +
                            b" >>\nstream\n" + stream + b"\nendstream"))
        return objects

    def save(self, path):
        """Write the document. Returns `path`."""
        first_page_obj = 5
        page_objects = self._page_objects(first_page_obj)
        kids = " ".join(f"{first_page_obj + i * 2} 0 R"
                        for i in range(len(self._pages)))

        stamp = datetime.datetime.now().strftime("D:%Y%m%d%H%M%S")
        objects = {
            1: b"<< /Type /Catalog /Pages 2 0 R >>",
            2: (f"<< /Type /Pages /Count {len(self._pages)} "
                f"/Kids [{kids}] >>").encode("latin-1"),
            3: b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica "
               b"/Encoding /WinAnsiEncoding >>",
            4: b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold "
               b"/Encoding /WinAnsiEncoding >>",
        }
        for num, body in page_objects:
            objects[num] = body

        info_num = max(objects) + 1
        objects[info_num] = (
            f"<< /Title ({_escape(self.title)}) "
            f"/Producer (cancer-genomics-pipelines webapp) "
            f"/CreationDate ({stamp}) >>").encode("latin-1")

        out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = {}
        for num in sorted(objects):
            offsets[num] = len(out)
            out += f"{num} 0 obj\n".encode("latin-1")
            out += objects[num]
            out += b"\nendobj\n"

        # Cross-reference table: readers rely on these byte offsets being
        # exact, so they are captured above as each object is emitted.
        xref_at = len(out)
        count = max(objects) + 1
        out += f"xref\n0 {count}\n".encode("latin-1")
        out += b"0000000000 65535 f \n"
        for num in range(1, count):
            out += f"{offsets[num]:010d} 00000 n \n".encode("latin-1")
        out += (f"trailer\n<< /Size {count} /Root 1 0 R "
                f"/Info {info_num} 0 R >>\nstartxref\n{xref_at}\n"
                f"%%EOF\n").encode("latin-1")

        with open(path, "wb") as fh:
            fh.write(out)
        return path
