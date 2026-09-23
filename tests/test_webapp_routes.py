# Created by Brainstorm, 2026.
"""
Web interface: every page renders, and paths stay inside the allowed roots.

The path confinement tests matter most. The form asks for server-side
paths, which means the browser can post any string it likes.
"""

import os
import re
import unittest

from webcase import WebCase


PAGES = ("/", "/new", "/new/dna", "/new/rna", "/new/hybrid",
         "/new/worksheet", "/panels", "/databases")


class TestPages(WebCase):

    def test_every_page_renders(self):
        for url in PAGES:
            self.assertEqual(self.client.get(url).status_code, 200, url)

    def test_the_first_entry_point_is_the_chooser_not_a_form(self):
        # "Start a run" on the home page once skipped the choice entirely
        # and dropped the operator into the DNA form, silently deciding
        # the assay for them.
        body = self.client.get("/").data.decode()
        match = re.search(r'href="([^"]+)"[^>]*>Start a run<', body)
        if match:
            self.assertEqual(match.group(1), "/new")

    def test_every_url_for_target_resolves(self):
        endpoints = {r.endpoint for r in self.app_module.app.url_map.iter_rules()}
        templates = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "webapp", "templates")
        for name in sorted(os.listdir(templates)):
            with open(os.path.join(templates, name)) as handle:
                src = handle.read()
            for target in set(re.findall(r"url_for\(\s*'([a-z_]+)'", src)):
                self.assertIn(target, endpoints, f"{name} -> {target}")

    def test_every_form_posts_to_a_post_capable_endpoint(self):
        methods = {r.endpoint: r.methods
                   for r in self.app_module.app.url_map.iter_rules()}
        templates = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "webapp", "templates")
        for name in sorted(os.listdir(templates)):
            with open(os.path.join(templates, name)) as handle:
                src = handle.read()
            for target in re.findall(
                    r'<form[^>]*method="post"[^>]*action="\{\{\s*url_for\(\s*\'([a-z_]+)\'',
                    src, re.I):
                self.assertIn("POST", methods.get(target, set()),
                              f"{name} -> {target}")


class TestPathConfinement(WebCase):

    def test_a_path_outside_the_allowed_roots_is_refused(self):
        with self.assertRaises(ValueError):
            self.app_module.safe_path("/etc/passwd", must_exist=True)

    def test_traversal_is_resolved_before_the_check(self):
        with self.assertRaises(ValueError):
            self.app_module.safe_path(os.path.join(self.tmp, "..", "..",
                                                   "etc", "passwd"))

    def test_a_path_inside_a_root_is_accepted(self):
        inside = self.touch("ok.txt")
        self.assertEqual(self.app_module.safe_path(inside, must_exist=True),
                         os.path.realpath(inside))


class TestArtefactAddressing(WebCase):

    def _run_with_outputs(self):
        out = os.path.join(self.tmp, "out")
        for rel in ("coverage/S_B.coverage.html", "pcgr/r.html",
                    "fastp_reports/S_B.fastp.html"):
            self.write(f"out/{rel}", "<html>x</html>")
        job = self.app_module.manager.submit(
            "x.py", ["py", "x.py"], {"patient_id": "P"},
            meta={"output_dir": out}, assay="dna")
        # The job writes its log and state file on a worker thread, and
        # those are artefacts too. Counting artefacts while it is still
        # running counts a moving target -- the cause of an intermittent
        # "5 != 4" here that had nothing to do with what the test checks.
        self.drain()
        return job, out

    def test_ids_survive_new_files_appearing_mid_run(self):
        # The reports are meant to be readable WHILE the run is going, so
        # the set grows under an already-rendered page. Addressed by
        # position, a link then served a different file than its label.
        import artefacts
        job, out = self._run_with_outputs()
        before = {i["id"]: i["label"]
                  for i in artefacts.collect(job.snapshot(), job.run_dir, out)}
        self.write("out/coverage/S_A.coverage.html", "<html>new</html>")
        after = {i["id"]: i["label"]
                 for i in artefacts.collect(job.snapshot(), job.run_dir, out)}
        for key, label in before.items():
            if key in after:
                self.assertEqual(after[key], label, "an id changed meaning")
        self.assertEqual(len(after), len(before) + 1)

    def test_resolve_refuses_a_path_outside_the_run(self):
        import artefacts
        evil = [{"path": "/etc/passwd", "action": "download", "id": "abc"}]
        path, item = artefacts.resolve(evil, "abc", self.tmp, self.tmp)
        self.assertIsNone(path)

    def test_every_link_on_the_job_page_serves(self):
        job, _ = self._run_with_outputs()
        body = self.client.get(f"/job/{job.id}").data.decode()
        links = re.findall(r'href="(/job/[^"]*/file/[a-f0-9]+)"', body)
        self.assertGreater(len(links), 0)
        for url in links:
            self.assertEqual(self.client.get(url).status_code, 200, url)

    def test_an_unknown_artefact_id_is_a_404(self):
        job, _ = self._run_with_outputs()
        self.assertEqual(
            self.client.get(f"/job/{job.id}/file/deadbeef0000").status_code,
            404)


if __name__ == "__main__":
    unittest.main()
