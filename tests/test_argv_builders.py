# Created by Brainstorm, 2026.
"""
Command construction: every flag the web interface emits must exist.

Command EXECUTION is not tested -- see tests/__init__.py for why. What is
tested is that a flag the app emits is a flag the target script accepts,
which is the failure that otherwise surfaces as a run aborting minutes in.
"""

import os
import re
import subprocess
import sys
import unittest

from helpers import TempCase, have_flask

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def script_flags(script):
    """Every long flag a script's --help advertises."""
    out = subprocess.run([sys.executable, os.path.join(ROOT, script),
                          "--help"], capture_output=True, text=True).stdout
    return set(re.findall(r"(--[a-z0-9][a-z0-9-]+)", out))


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestEmittedFlagsExist(TempCase):

    def _assert_known(self, script, argv):
        known = script_flags(script)
        used = {token for token in argv if token.startswith("--")}
        self.assertEqual(used - known, set(),
                         f"{script} does not accept these")

    def test_dna_pipeline_argv(self):
        from jobs import build_pipeline_argv
        argv = build_pipeline_argv("py", "comprehensive_variant_calling.py", {
            "output_dir": "/o", "reference": "hg38", "reference_dir": "/rd",
            "input_dir": "/i", "tumour_sample": "T", "tumour_r1": "/1",
            "tumour_r2": "/2", "normal_sample": "N", "normal_r1": "/3",
            "normal_r2": "/4", "panel": "p", "panel_file": "a",
            "cosmic": "/c", "dbsnp": "/d", "germline_resource": "/g",
            "panel_of_normals": "/p", "contamination_resource": "/cc",
            "msi_models": "/m", "intervals": "/iv", "known_indels": "/k",
            "interval_padding": "100", "coverage_bed": "/cb",
            "coverage_min_depth": "20", "min_depth": "50",
            "min_allele_fraction": "0.05", "mutect2_extra_args": "-z",
            "threads": "8", "two_pass": "on", "resume": "on",
            "dry_run": "on", "skip_steps": "qc"})
        self._assert_known("comprehensive_variant_calling.py", argv)

    def test_rna_pipeline_argv(self):
        from jobs import build_rna_argv
        argv = build_rna_argv("py", "fusion_calling.py", {
            "output_dir": "/o", "reference": "/r.fa", "gtf": "/g.gtf",
            "panel": "p", "panel_file": "a b", "star_index": "/i",
            "reference_dir": "/rd", "read_length": "150",
            "arriba_resources": "/ar", "manifest": "/m.json",
            "fusion_min_confidence": "medium", "min_fusion_reads": "3",
            "strandedness": "unstranded", "rna_min_reads_millions": "20",
            "rna_min_unique_mapped_pct": "60", "min_read_length": "35",
            "star_extra_args": "-x", "arriba_extra_args": "-y",
            "threads": "4", "dry_run": "on", "skip_steps": "qc"})
        self._assert_known("fusion_calling.py", argv)

    def test_pcgr_argv_for_both_shapes(self):
        from jobs import build_pcgr_argv, build_rna_pcgr_argv
        vcf = self.touch("calls.vcf.gz")
        dna, _ = build_pcgr_argv(
            "py", ROOT, {"tumour": "T", "filtered_vcf": vcf},
            {"pcgr_refdata_dir": "/rd", "vep_dir": "/v",
             "pcgr_assay": "TARGETED", "pcgr_tumour_site": "9",
             "pcgr_estimate_tmb": "on", "pcgr_estimate_msi": "on",
             "pcgr_estimate_signatures": "on"}, "/o")
        self._assert_known("pcgr_report.py", dna)

        fusions = self.touch("f.tsv")
        rna, _ = build_rna_pcgr_argv(
            "py", ROOT,
            {"samples": {"S": {"fusion_report": {"pcgr_tsv": fusions}}}},
            {"pcgr_refdata_dir": "/rd", "pcgr_tumour_site": "9",
             "min_fusion_reads": "3"}, "/o")
        self._assert_known("pcgr_report.py", rna)


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestPcgrInputShapes(TempCase):

    def test_fusions_alone_need_no_vep_cache(self):
        # PCGR requires a VEP cache only when given a VCF; a fusion-only
        # report has nothing to annotate.
        from jobs import build_rna_pcgr_argv
        fusions = self.touch("f.tsv")
        argv, _ = build_rna_pcgr_argv(
            "py", ROOT,
            {"samples": {"S": {"fusion_report": {"pcgr_tsv": fusions}}}},
            {"pcgr_refdata_dir": "/rd"}, "/o")
        self.assertIn("--input-rna-fusion", argv)
        self.assertNotIn("--vep-dir", argv)

    def test_no_fusions_means_no_pcgr_run(self):
        # PCGR rejects an empty fusion file; "nothing found" is a result.
        from jobs import build_rna_pcgr_argv
        argv, note = build_rna_pcgr_argv(
            "py", ROOT, {"samples": {"S": {}}},
            {"pcgr_refdata_dir": "/rd"}, "/o")
        self.assertIsNone(argv)
        self.assertIn("nothing for PCGR", note)


class TestPcgrCommand(unittest.TestCase):

    def test_command_builds_for_vcf_only_fusion_only_and_both(self):
        from pcgr_report import build_pcgr_command
        fusion_only = build_pcgr_command(None, "/o", "S", "/rd",
                                         rna_fusion="/f.tsv")
        self.assertIn("--input_rna_fusion", fusion_only)
        self.assertNotIn("--input_vcf", fusion_only)

        both = build_pcgr_command("/in.vcf", "/o", "S", "/rd",
                                  vep_dir="/vep", rna_fusion="/f.tsv")
        self.assertIn("--input_vcf", both)
        self.assertIn("--input_rna_fusion", both)

    def test_tmb_denominator_is_only_emitted_when_known(self):
        from pcgr_report import build_pcgr_command
        without = build_pcgr_command("/in.vcf", "/o", "S", "/rd",
                                     vep_dir="/v")
        self.assertNotIn("--effective_target_size_mb", without)
        with_size = build_pcgr_command("/in.vcf", "/o", "S", "/rd",
                                       vep_dir="/v",
                                       effective_target_size_mb=1.94)
        self.assertIn("1.94", with_size)


if __name__ == "__main__":
    unittest.main()
