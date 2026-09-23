# Created by Brainstorm, 2026.
"""
The fusion report: ranking by clinical salience, and the PCGR handoff.

Arriba ranks by evidence quality, which puts a high-confidence read-through
between two housekeeping genes above a medium-confidence EML4--ALK. The
reordering here is the whole point of the module, so it is what these tests
pin down.
"""

import unittest

from helpers import ARRIBA_HEADER, TempCase, arriba_row, read_json, read_text
import fusion_report as fr


class TestClassification(unittest.TestCase):

    def _row(self, **kwargs):
        base = {"gene1": "A", "gene2": "B", "breakpoint1": "chr1:1",
                "breakpoint2": "chr2:2", "split_reads1": "5",
                "split_reads2": "4", "discordant_mates": "2",
                "confidence": "medium"}
        base.update(kwargs)
        return fr.classify(base)

    def test_actionable_gene_is_recognised(self):
        self.assertEqual(self._row(gene1="EML4", gene2="ALK")
                         ["_actionable_genes"], ["ALK"])

    def test_same_gene_event_is_flagged_with_its_note(self):
        row = self._row(gene1="MET", gene2="MET")
        self.assertTrue(row["_same_gene_event"])
        self.assertIn("exon 14", row["_splice_note"])

    def test_supporting_reads_include_discordant_mates(self):
        # The report's own figure is TOTAL evidence for the event. It is
        # deliberately not the same number PCGR is given -- see
        # TestPcgrHandoff.test_pcgr_split_reads_exclude_discordant_mates.
        row = self._row(split_reads1="6", split_reads2="5",
                        discordant_mates="99")
        self.assertEqual(row["_supporting_reads"], 110)

    def test_intergenic_gene_names_are_parsed(self):
        row = self._row(gene1="FOO(123),BAR", gene2="ALK")
        self.assertEqual(row["_pair"], "FOO--ALK")

    def test_confidence_bar_marks_rather_than_drops(self):
        row = self._row(confidence="low")
        self.assertFalse(row["_meets_confidence"])
        self.assertEqual(row["_confidence"], "low")


class TestRanking(TempCase):

    def _report(self, rows, **kwargs):
        path = self.write("f.tsv", ARRIBA_HEADER + "".join(rows))
        return fr.run_fusion_report(path, "S", self.path("out/x"), **kwargs)

    def test_actionable_outranks_a_higher_confidence_artefact(self):
        result = self._report([
            arriba_row("GAPDH", "ACTB", confidence="high",
                       split1=40, split2=38),
            arriba_row("EML4", "ALK", confidence="medium"),
        ])
        data = read_json(result["json"])
        self.assertEqual(data["fusions"][0]["_pair"], "EML4--ALK")

    def test_same_gene_splice_event_outranks_a_plain_actionable_fusion(self):
        result = self._report([
            arriba_row("EML4", "ALK", confidence="medium"),
            arriba_row("MET", "MET", confidence="medium"),
        ])
        data = read_json(result["json"])
        self.assertEqual(data["fusions"][0]["_pair"], "MET--MET")

    def test_summary_counts_are_consistent(self):
        result = self._report([
            arriba_row("EML4", "ALK", confidence="medium"),
            arriba_row("FOO", "BAR", confidence="low"),
        ])
        summary = result["summary"]
        self.assertEqual(summary["total_called"], 2)
        self.assertEqual(summary["at_or_above_confidence"], 1)
        self.assertEqual(summary["involving_actionable_gene"], 1)


class TestPcgrHandoff(TempCase):

    def _convert(self, rows, min_confidence="medium"):
        path = self.write("f.tsv", ARRIBA_HEADER + "".join(rows))
        result = fr.run_fusion_report(path, "S", self.path("out/y"),
                                      min_confidence=min_confidence)
        return result

    def test_schema_matches_what_pcgr_requires(self):
        result = self._convert([arriba_row("EML4", "ALK")])
        header = read_text(result["pcgr_tsv"]).splitlines()[0].split("\t")
        self.assertEqual(header, ["FusionGene", "LeftBreakpoint",
                                  "RightBreakpoint", "SplitReads"])

    def test_gene_pair_uses_the_double_dash_separator(self):
        result = self._convert([arriba_row("EML4", "ALK")])
        body = read_text(result["pcgr_tsv"]).splitlines()[1]
        self.assertTrue(body.split("\t")[0] == "EML4--ALK")

    def test_calls_below_the_bar_are_excluded_from_pcgr(self):
        # Unlike the HTML report, which marks them: a weak call handed to
        # PCGR becomes an entry in a clinical interpretation with no way
        # to show it was weak.
        result = self._convert([arriba_row("EML4", "ALK", confidence="low")])
        self.assertIsNone(result["pcgr_tsv"])

    def test_pcgr_split_reads_exclude_discordant_mates(self):
        # PCGR's column is named SplitReads and its
        # --fusion_min_split_reads threshold is written against split
        # reads; handing it the report's total would push every fusion
        # past that threshold.
        result = self._convert([arriba_row("EML4", "ALK", split1=6,
                                           split2=5, discordant=99)])
        row = read_text(result["pcgr_tsv"]).splitlines()[1].split("\t")
        self.assertEqual(row[3], "11")

    def test_no_qualifying_fusions_writes_no_file(self):
        # PCGR rejects an EMPTY fusion file outright, so "nothing found"
        # must produce no file rather than a header-only one.
        result = self._convert([])
        self.assertIsNone(result["pcgr_tsv"])


class TestDegenerateInput(TempCase):

    def test_header_only_file_is_handled(self):
        path = self.write("e.tsv", ARRIBA_HEADER)
        result = fr.run_fusion_report(path, "S", self.path("out/a"))
        self.assertTrue(result["ok"])
        self.assertEqual(result["summary"]["total_called"], 0)

    def test_missing_file_reports_rather_than_raises(self):
        result = fr.run_fusion_report(self.path("nope.tsv"), "S",
                                      self.path("out/b"))
        self.assertFalse(result["ok"])
        self.assertIn("could not read", result["error"])

    def test_truncated_row_does_not_crash(self):
        path = self.write("t.tsv", ARRIBA_HEADER + "EML4\n")
        result = fr.run_fusion_report(path, "S", self.path("out/c"))
        self.assertTrue(result["ok"])


if __name__ == "__main__":
    unittest.main()
