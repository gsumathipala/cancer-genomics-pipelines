# Created by Brainstorm, 2026.
"""
Shared fixtures: synthetic inputs small enough to build in a temp dir.

Everything here is fabricated. No test reads real sequencing data, a real
reference or a real run directory -- partly for speed, and partly because
a test suite that needed patient data would be a test suite nobody could
run on a laptop.
"""

import os
import sys
import tempfile
import unittest

# The modules under test sit at the repository root and in webapp/. Tests
# are run from anywhere, so the paths are resolved from this file rather
# than from the working directory.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBAPP = os.path.join(ROOT, "webapp")
for path in (ROOT, WEBAPP):
    if path not in sys.path:
        sys.path.insert(0, path)


def have_flask():
    """Flask backs the web interface only; its absence skips those tests."""
    try:
        import flask                       # noqa: F401
        return True
    except ImportError:
        return False


class TempCase(unittest.TestCase):
    """A TestCase with a scratch directory that cleans itself up."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

    def path(self, *parts):
        """An absolute path inside the scratch directory."""
        full = os.path.join(self.tmp, *parts)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        return full

    def write(self, name, text):
        """Write a file into the scratch directory and return its path."""
        full = self.path(name)
        with open(full, "w", encoding="utf-8") as handle:
            handle.write(text)
        return full

    def touch(self, name):
        return self.write(name, "")


# --- synthetic inputs ---------------------------------------------------

def bed(rows):
    """A BED body from (chrom, start, end) triples."""
    return "".join(f"{c}\t{s}\t{e}\tr{i}\n"
                   for i, (c, s, e) in enumerate(rows))


def fai(rows):
    """A .fai body from (name, length) pairs. Only column 1 is ever read."""
    return "".join(f"{name}\t{length}\t0\t60\t61\n" for name, length in rows)


# Arriba's column order, trimmed to what fusion_report.py actually reads.
ARRIBA_HEADER = ("gene1\tgene2\tbreakpoint1\tbreakpoint2\tsite1\tsite2\t"
                 "type\tsplit_reads1\tsplit_reads2\tdiscordant_mates\t"
                 "confidence\treading_frame\n")


def arriba_row(gene1, gene2, bp1="chr1:100", bp2="chr2:200",
               split1=6, split2=5, discordant=3, confidence="medium",
               kind="translocation", frame="in-frame"):
    return (f"{gene1}\t{gene2}\t{bp1}\t{bp2}\tsplice-site\tsplice-site\t"
            f"{kind}\t{split1}\t{split2}\t{discordant}\t{confidence}\t"
            f"{frame}\n")


def star_log(input_reads=1000000, unique_pct=85.0, too_short_pct=5.0,
             chimeric=120, splices=50000):
    """A STAR Log.final.out carrying the fields the QC report parses."""
    return (
        "                                 Started job on |\tJan 01\n"
        f"                          Number of input reads |\t{input_reads}\n"
        f"                      Uniquely mapped reads % |\t{unique_pct}%\n"
        f"       % of reads mapped to multiple loci |\t5.00%\n"
        f"            % of reads unmapped: too short |\t{too_short_pct}%\n"
        f"                     Number of splices: Total |\t{splices}\n"
        f"           Number of splices: Annotated (sjdb) |\t{splices}\n"
        f"                     Number of chimeric reads |\t{chimeric}\n"
        "                     Average input read length |\t200\n"
    )
