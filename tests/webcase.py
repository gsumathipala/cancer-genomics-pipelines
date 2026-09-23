# Created by Brainstorm, 2026.
"""
A TestCase that boots the web interface against a scratch directory.

The app reads module-level config and owns a JobManager, so each test
class gets its own runs directory and allowed roots. Jobs are submitted
with dry_run set, so nothing executes a real tool.
"""

import os
import unittest

from helpers import TempCase, have_flask


@unittest.skipUnless(have_flask(), "Flask is not installed")
class WebCase(TempCase):

    def setUp(self):
        super().setUp()
        import sys
        sys.argv = ["app.py"]
        import app as app_module
        self.app_module = app_module
        self.runs = os.path.join(self.tmp, "runs")
        app_module.app.config["ALLOWED_ROOTS"] = [
            os.path.realpath(self.tmp), os.path.expanduser("~/data")]
        app_module.app.config["RUNS_DIR"] = self.runs
        app_module.manager = app_module.jobs.JobManager(self.runs)
        self.client = app_module.app.test_client()

    # -- fixtures a run needs ------------------------------------------
    def fastqs(self, *samples):
        made = {}
        for sample in samples:
            made[sample] = (self.touch(f"fq/{sample}_R1_001.fastq.gz"),
                            self.touch(f"fq/{sample}_R2_001.fastq.gz"))
        return made

    def reference(self):
        ref = self.write("ref.fa", ">chr1\nACGT\n")
        self.write("ref.fa.fai", "chr1\t4\t6\t4\t5\n")
        return ref

    def panel_bed(self):
        return self.write("panel.bed", "chr1\t1\t3\tR\n")

    def gtf(self):
        return self.write("g.gtf",
                          'chr1\tT\tgene\t1\t3\t.\t+\t.\tgene_id "x"; '
                          'gene_type "rRNA";\n')

    def star_index(self):
        for name in ("SA", "SAindex", "Genome", "genomeParameters.txt"):
            self.touch(f"idx/{name}")
        return os.path.join(self.tmp, "idx")

    def job_count(self):
        return len(self.app_module.manager.all())
