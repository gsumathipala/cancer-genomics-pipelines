# Created by Brainstorm, 2026.
"""
RNA library adequacy: the report that decides whether a negative counts.

A degraded FFPE library produces a clean, well-formed, EMPTY fusion table
that is indistinguishable from a true negative. These tests pin the checks
that tell the two apart.
"""

import json
import unittest

from helpers import TempCase, star_log
import rna_qc_report as rq
import align_rna


class TestStarLogParsing(TempCase):

    def test_percentages_and_counts_are_typed_correctly(self):
        path = self.write("Log.final.out", star_log(unique_pct=85.5,
                                                    chimeric=116))
        stats = align_rna.parse_star_log(path)
        self.assertEqual(stats["Uniquely mapped reads %"], 85.5)
        self.assertEqual(stats["Number of chimeric reads"], 116)
        self.assertIsInstance(stats["Number of input reads"], int)

    def test_a_missing_log_returns_empty_not_an_exception(self):
        self.assertEqual(align_rna.parse_star_log(self.path("nope")), {})


class TestChecks(unittest.TestCase):

    def _states(self, measured, thresholds=None):
        checks = rq.evaluate(measured, thresholds or {})
        return {c["name"]: c["state"] for c in checks}

    def test_chimeric_check_can_actually_fail(self):
        # It was written `chimeric is not None` INSIDE `if chimeric is not
        # None`, so it was always true and the one failure it exists to
        # catch -- chimeric detection switched off -- could never raise.
        self.assertEqual(self._states({"chimeric_reads": 0})
                         ["Chimeric detection active"], "concern")
        self.assertEqual(self._states({"chimeric_reads": 120})
                         ["Chimeric detection active"], "pass")

    def test_library_size_floor_is_applied(self):
        states = self._states({"input_reads_millions": 4.0},
                              {"min_reads_millions": 20.0})
        self.assertEqual(states["Library size"], "concern")

    def test_degraded_library_signature_is_flagged(self):
        self.assertEqual(self._states({"too_short_pct": 44.0})
                         ["Reads lost as 'too short'"], "concern")

    def test_mitochondrial_and_rrna_fractions(self):
        states = self._states({"mito_fraction": 0.41, "rrna_fraction": 0.5})
        self.assertEqual(states["Mitochondrial fraction"], "concern")
        self.assertEqual(states["rRNA fraction"], "concern")


class TestVerdict(unittest.TestCase):

    def test_nothing_measured_is_unknown_not_a_pass(self):
        verdict = rq.verdict([])
        self.assertEqual(verdict["state"], "unknown")
        self.assertIn("uninterpreted", verdict["text"])

    def test_any_concern_blocks_a_reportable_negative(self):
        verdict = rq.verdict([{"name": "Library size", "state": "concern",
                               "detail": ""}])
        self.assertEqual(verdict["state"], "concern")
        self.assertIn("NOT", verdict["text"])

    def test_all_clear_says_a_negative_is_a_negative(self):
        verdict = rq.verdict([{"name": "Library size", "state": "pass",
                               "detail": ""}])
        self.assertEqual(verdict["state"], "ok")


class TestReport(TempCase):

    def test_report_runs_with_nothing_measurable(self):
        result = rq.run_rna_qc("S", None, None, None, self.path("out/x"))
        self.assertTrue(result["ok"])
        self.assertEqual(result["metrics"]["verdict"]["state"], "unknown")

    def test_unmeasurable_things_are_named_not_omitted(self):
        result = rq.run_rna_qc("S", None, None, None, self.path("out/y"))
        joined = " ".join(result["metrics"]["not_measured"])
        # Their absence must not read as a pass.
        self.assertIn("duplication rate", joined)
        self.assertIn("DV200", joined)

    def test_rrna_intervals_are_derived_from_the_gtf(self):
        gtf = self.write("g.gtf",
                         'chr1\tT\tgene\t1\t500\t.\t+\t.\tgene_id "r"; '
                         'gene_type "rRNA";\n'
                         'chr1\tT\tgene\t900\t999\t.\t+\t.\tgene_id "p"; '
                         'gene_type "protein_coding";\n')
        path, count = rq.rrna_intervals(gtf, self.path("out/r.bed"))
        self.assertEqual(count, 1)
        # GTF is 1-based inclusive; BED is 0-based half-open.
        self.assertEqual(open(path).read().strip(), "chr1\t0\t500")

    def test_a_gtf_with_no_rrna_yields_no_interval_file(self):
        gtf = self.write("g2.gtf",
                         'chr1\tT\tgene\t1\t50\t.\t+\t.\tgene_type "x";\n')
        path, count = rq.rrna_intervals(gtf, self.path("out/s.bed"))
        self.assertIsNone(path)
        self.assertEqual(count, 0)


if __name__ == "__main__":
    unittest.main()
