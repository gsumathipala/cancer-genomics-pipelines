# Created by Brainstorm, 2026.
"""
A TestCase that boots the web interface against a scratch directory.

The app reads module-level config and owns a JobManager, so each test
class gets its own runs directory and allowed roots. Jobs are submitted
with dry_run set, so nothing executes a real tool.
"""

import os
import time
import unittest

from helpers import TempCase, have_flask


def closing_client_class():
    """
    A test client that reads each response fully and then closes it.

    A route that serves a report calls send_file() with a PATH, so it is
    Werkzeug that opens the file and attaches it to the response. Under a
    real server that handle is released when the response has been sent --
    the WSGI server closes the iterable. The TEST client never sends
    anything, so unless the response is closed the handle survives until
    the garbage collector happens to notice, which is what raised
    ResourceWarning during the route tests.

    Forcing buffered=True makes Werkzeug consume the response and close it
    inside the request call, so a handle cannot outlive the test that
    opened it. This costs nothing here -- the fixtures are a few hundred
    bytes -- and it keeps a real leak visible: with this in place, a
    ResourceWarning from a route test means the APP held a handle open,
    not the harness.
    """
    from flask.testing import FlaskClient

    class ClosingClient(FlaskClient):
        def open(self, *args, **kwargs):
            kwargs.setdefault("buffered", True)
            return super().open(*args, **kwargs)

    return ClosingClient


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
        app_module.app.test_client_class = closing_client_class()
        self.client = app_module.app.test_client()

    def tearDown(self):
        # Jobs run on a daemon thread and keep writing into the run
        # directory. TempCase removes that directory, so a test that ends
        # while a dry run is still finishing deletes files from under the
        # worker -- an intermittent teardown error that has nothing to do
        # with the test that happens to be running when it lands. Wait for
        # the queue to empty first.
        self.drain()
        super().tearDown()

    def drain(self, timeout=10.0):
        """Block until no job is queued, running or reporting."""
        manager = getattr(self.app_module, "manager", None)
        if manager is None:
            return
        deadline = time.monotonic() + timeout
        while manager.busy() and time.monotonic() < deadline:
            time.sleep(0.02)
        worker = getattr(manager, "_worker", None)
        if worker is not None and worker.is_alive():
            worker.join(timeout=max(0.0, deadline - time.monotonic()))

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
