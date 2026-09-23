# Created by Brainstorm, 2026.
"""
The multi-sample worksheet.

Two properties are load-bearing and both are tested by COUNTING JOBS:
nothing is queued unless every row validates, and each row carries its own
patient. A worksheet that queued four of twelve rows and then rejected the
fifth would leave a half-started batch with no clear way back.
"""

import os
import time
import unittest

from webcase import WebCase


class WorksheetCase(WebCase):

    def setUp(self):
        super().setUp()
        self.fq = self.fastqs("P001", "P002", "P003")
        self.shared = {
            "output_dir": os.path.join(self.runs, "RUN"),
            "reference": self.reference(),
            "panel_bed": self.panel_bed(),
            "gtf": self.gtf(),
            "star_index": self.star_index(),
            "pcgr_tumour_site": "15",
            "referring_clinician": "Dr Shared",
            "specimen_type": "FFPE block",
            "threads": "2", "dry_run": "on", "form_rendered": "1",
            "minimal_run": "on",
        }

    def row(self, index, sample, assay="dna", **extra):
        data = {
            f"row-{index}-include": "on",
            f"row-{index}-sample": sample,
            f"row-{index}-assay": assay,
            f"row-{index}-patient_id": f"MRN-{index}",
            f"row-{index}-r1": self.fq[sample][0],
            f"row-{index}-r2": self.fq[sample][1],
        }
        if assay == "hybrid":
            data[f"row-{index}-rna_r1"] = self.fq[sample][0]
            data[f"row-{index}-rna_r2"] = self.fq[sample][1]
        for key, value in extra.items():
            data[f"row-{index}-{key}"] = value
        return data

    def post(self, *rows, **shared_overrides):
        payload = dict(self.shared)
        payload.update(shared_overrides)
        for row in rows:
            payload.update(row)
        return self.client.post("/submit/worksheet", data=payload)

    def settle(self, seconds=15):
        """Wait for the serial queue to drain the dry runs."""
        deadline = time.time() + seconds
        while time.time() < deadline:
            if not self.app_module.manager.busy():
                return
            time.sleep(0.4)


class TestRowParsing(WorksheetCase):

    def test_blank_and_unticked_rows_are_ignored(self):
        parse = self.app_module.parse_worksheet_rows
        self.assertEqual(parse({"row-0-sample": "", "row-0-include": "on"}), [])
        self.assertEqual(parse({"row-0-sample": "A"}), [])   # not ticked
        self.assertEqual(len(parse({"row-0-sample": "A",
                                    "row-0-include": "on"})), 1)

    def test_input_dir_is_derived_from_the_row_files(self):
        # The handlers need one even when the files are named explicitly.
        row = {"sample": "A", "assay": "dna",
               "r1": "/data/run/A_R1.fastq.gz", "r2": "/data/run/A_R2.fastq.gz"}
        values = self.app_module.worksheet_row_to_form(row, {"output_dir": "/o"})
        self.assertEqual(values["input_dir"], "/data/run")

    def test_auto_discover_is_never_inherited_by_a_row(self):
        # It would override the row's explicit files and reintroduce the
        # ambiguity the worksheet exists to remove.
        values = self.app_module.worksheet_row_to_form(
            {"sample": "A", "assay": "dna", "r1": "/d/A_R1.fastq.gz"},
            {"output_dir": "/o", "auto_discover": "on", "batch_mode": "on"})
        self.assertNotIn("auto_discover", values)
        self.assertNotIn("batch_mode", values)

    def test_each_row_gets_its_own_output_directory(self):
        a = self.app_module.worksheet_row_to_form(
            {"sample": "A", "assay": "dna"}, {"output_dir": "/o"})
        b = self.app_module.worksheet_row_to_form(
            {"sample": "B", "assay": "dna"}, {"output_dir": "/o"})
        self.assertNotEqual(a["output_dir"], b["output_dir"])


class TestNothingQueuedOnError(WorksheetCase):

    def _refuses(self, *rows, **shared):
        before = self.job_count()
        response = self.post(*rows, **shared)
        self.assertEqual(self.job_count(), before,
                         "rows were queued despite an invalid worksheet")
        return response

    def test_a_bad_path_in_one_row_stops_the_whole_worksheet(self):
        bad = self.row(1, "P002")
        bad["row-1-r1"] = "/nonexistent/x.fastq.gz"
        self._refuses(self.row(0, "P001"), bad)

    def test_duplicate_sample_names_are_refused(self):
        # Rows write into directories named after the sample, so duplicates
        # would overwrite each other's results.
        response = self._refuses(self.row(0, "P001"),
                                 {**self.row(1, "P002"),
                                  "row-1-sample": "P001"})
        self.assertIn("more than once", response.data.decode())

    def test_a_row_with_no_patient_identifier_is_refused(self):
        row = self.row(0, "P001")
        row["row-0-patient_id"] = ""
        self._refuses(row, patient_id="", specimen_id="")

    def test_an_invalid_tumour_site_names_its_row(self):
        response = self._refuses(self.row(0, "P001", pcgr_tumour_site="99"))
        self.assertIn("Row 1", response.data.decode())

    def test_an_empty_worksheet_is_refused(self):
        self._refuses()


class TestPerRowIdentity(WorksheetCase):

    def test_each_row_becomes_a_run_with_its_own_patient(self):
        response = self.post(
            self.row(0, "P001", patient_id="MRN-A", specimen_id="SP-A"),
            self.row(1, "P002", patient_id="MRN-B", specimen_id="SP-B"))
        self.assertEqual(response.status_code, 302)
        self.settle()
        patients = {j.meta.get("tumour_sample"): j.patient
                    for j in self.app_module.manager.all()}
        self.assertEqual(patients["P001"]["patient_id"], "MRN-A")
        self.assertEqual(patients["P002"]["patient_id"], "MRN-B")

    def test_blank_row_fields_inherit_the_header(self):
        self.post(self.row(0, "P001", specimen_type="Bone marrow"),
                  self.row(1, "P002"))
        self.settle()
        by_sample = {j.meta.get("tumour_sample"): j.patient
                     for j in self.app_module.manager.all()}
        self.assertEqual(by_sample["P001"]["specimen_type"], "Bone marrow")
        self.assertEqual(by_sample["P002"]["specimen_type"], "FFPE block")
        self.assertEqual(by_sample["P002"]["referring_clinician"], "Dr Shared")

    def test_tumour_site_is_per_row_with_a_header_fallback(self):
        # PCGR tiers actionability AGAINST the site: a shared site would
        # tier every sample as though it came from the same organ.
        self.post(self.row(0, "P001", pcgr_tumour_site="6"),
                  self.row(1, "P002"))
        self.settle()
        sites = {j.meta.get("tumour_sample"): j.meta.get("pcgr_tumour_site")
                 for j in self.app_module.manager.all()}
        self.assertEqual(str(sites["P001"]), "6")
        self.assertEqual(str(sites["P002"]), "15")

    def test_a_hybrid_row_queues_two_linked_runs(self):
        self.post(self.row(0, "P001", assay="hybrid"))
        self.settle(20)
        jobs = self.app_module.manager.all()
        assays = sorted(j.snapshot()["assay"] for j in jobs)
        self.assertEqual(assays, ["dna", "rna"])
        self.assertTrue(all(j.meta.get("paired_run_id") for j in jobs))

    def test_mixed_assays_in_one_worksheet(self):
        self.post(self.row(0, "P001", assay="dna"),
                  self.row(1, "P002", assay="rna"))
        self.settle()
        assays = sorted(j.snapshot()["assay"]
                        for j in self.app_module.manager.all())
        self.assertEqual(assays, ["dna", "rna"])


class TestAnalysisToggles(WorksheetCase):

    def test_the_page_ships_tmb_and_signatures_ticked(self):
        # form_rendered means "take the boxes literally". The worksheet
        # posted that marker while having no boxes, so every toggle read as
        # unticked and every run silently produced no TMB.
        body = self.client.get("/new/worksheet").data.decode()
        for name in ("pcgr_estimate_tmb", "pcgr_estimate_signatures"):
            block = body.split(f'name="{name}"')[1][:40]
            self.assertIn("checked", block, name)

    def test_toggles_reach_the_queued_run(self):
        self.post(self.row(0, "P001"),
                  pcgr_estimate_tmb="on", pcgr_estimate_signatures="on")
        self.settle()
        job = self.app_module.manager.all()[0]
        self.assertEqual(job.meta.get("pcgr_estimate_tmb"), "on")


class TestScan(WorksheetCase):

    def test_scanning_a_folder_fills_one_row_per_sample(self):
        response = self.client.post("/worksheet/scan", data={
            "scan_dir": os.path.join(self.tmp, "fq"), "scan_assay": "dna"})
        body = response.data.decode()
        self.assertEqual(response.status_code, 200)
        for sample in ("P001", "P002", "P003"):
            self.assertIn(sample, body)

    def test_scanning_outside_the_allowed_roots_is_refused(self):
        response = self.client.post("/worksheet/scan",
                                    data={"scan_dir": "/etc"})
        self.assertIn("outside the allowed", response.data.decode())


if __name__ == "__main__":
    unittest.main()
