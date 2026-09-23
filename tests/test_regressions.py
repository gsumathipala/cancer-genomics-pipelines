# Created by Brainstorm, 2026.
"""
One named test per bug that actually reached this code.

The value of a regression test is knowing WHICH failure it prevents, so
every case below names its bug and, where the bug was a silent one, says
what the user would have seen instead of an error. A silent failure is
the class this pipeline treats as most serious: a wrong answer that looks
like a right one costs more than a crash, because nobody goes looking.
"""

import os
import sys
import unittest

from helpers import TempCase, have_flask, star_log
import panel_profiles


class TestFalsyValuesSurvive(unittest.TestCase):
    """
    BUG: `if value:` treated a legitimate 0 as "not set".

    A minimum allele fraction of 0 (report everything) and a padding of 0
    (use the BED exactly) are real, deliberate settings. Tested for truth
    rather than for None, they were silently replaced by the default --
    the run then filtered variants the user had asked to keep.
    """

    def test_zero_is_applied_not_dropped(self):
        class Args:
            min_allele_fraction = None
            interval_padding = None
        profile = {"id": "p", "name": "P",
                   "settings": {"min_allele_fraction": 0.0,
                                "interval_padding": 0}}
        args = Args()
        applied, _ = panel_profiles.apply_profile(args, profile, set())
        self.assertEqual(args.min_allele_fraction, 0.0)
        self.assertEqual(args.interval_padding, 0)
        self.assertEqual(len(applied), 2)


class TestRnaProfilesCarryNoTargetBed(unittest.TestCase):
    """
    BUG: DNA-only settings leaked into the RNA engine.

    A capture BED describes where DNA was enriched. Applied to an RNA run
    it means nothing, but the flag existed, so it was passed and quietly
    narrowed nothing while looking like it narrowed something.
    """

    def test_rna_settings_exclude_dna_only_keys(self):
        for key in ("panel_bed", "coverage_bed", "interval_padding",
                    "min_allele_fraction"):
            self.assertNotIn(key, panel_profiles.RNA_SETTINGS)

    def test_rna_settings_include_the_library_floors(self):
        for key in ("rna_min_reads_millions", "rna_min_unique_mapped_pct"):
            self.assertIn(key, panel_profiles.RNA_SETTINGS)


class TestChimericOutputMismatch(TempCase):
    """
    BUG: STAR reported 116 chimeric reads and wrote 0 into the BAM.

    The fusion caller reads only the BAM, so it found nothing and exited
    0. The run produced a clean, empty, well-formed fusion table for a
    specimen that may have carried a fusion -- identical in every
    downstream artefact to a true negative.
    """

    def setUp(self):
        super().setUp()
        import align_rna
        self.align_rna = align_rna
        self._real = align_rna.bam_has_chimeric_alignments

    def tearDown(self):
        self.align_rna.bam_has_chimeric_alignments = self._real
        super().tearDown()

    def test_disagreement_is_reported(self):
        self.align_rna.bam_has_chimeric_alignments = lambda bam: False
        notes = self.align_rna.verify_chimeric_output(
            {"bam": "/x/Aligned.bam", "prefix": "/x/"},
            {"Number of chimeric reads": 116})
        self.assertEqual(len(notes), 1)
        self.assertIn("116", notes[0])
        # It must say the specimen is not the explanation, or the reader
        # will file it as a true negative anyway.
        self.assertIn("not a property of the", notes[0])

    def test_agreement_is_silent(self):
        self.align_rna.bam_has_chimeric_alignments = lambda bam: True
        self.assertEqual(self.align_rna.verify_chimeric_output(
            {"bam": "/x/a.bam", "prefix": "/x/"},
            {"Number of chimeric reads": 116}), [])

    def test_no_claim_means_nothing_to_check(self):
        # Zero chimeric reads with zero in the BAM is consistent, not a bug.
        self.align_rna.bam_has_chimeric_alignments = lambda bam: False
        self.assertEqual(self.align_rna.verify_chimeric_output(
            {"bam": "/x/a.bam", "prefix": "/x/"},
            {"Number of chimeric reads": 0}), [])

    def test_junction_file_is_named_when_it_holds_the_evidence(self):
        # The junctions surviving proves detection worked and only the
        # in-BAM channel failed -- that distinction is the whole fix.
        self.write("Chimeric.out.junction", "chr1\t1\t+\n")
        self.align_rna.bam_has_chimeric_alignments = lambda bam: False
        notes = self.align_rna.verify_chimeric_output(
            {"bam": os.path.join(self.tmp, "a.bam"),
             "prefix": os.path.join(self.tmp, "")},
            {"Number of chimeric reads": 116})
        self.assertIn("confirms detection worked", notes[0])


class TestStarIndexDiscovery(TempCase):
    """
    BUG: the default index path was constructed, not discovered.

    The installer always appends the read length (star_hg38_150). The
    fallback invented "star_hg38", found nothing there, and told the user
    to build an index that already existed one directory along.
    """

    def test_existing_index_is_found_without_a_read_length(self):
        import align_rna
        os.makedirs(os.path.join(self.tmp, "star_hg38_150"))
        # Every sentinel is the sentinel discovery looks for -- see
        # installed_star_indexes(); an index without it is incomplete.
        for sentinel in align_rna.INDEX_SENTINELS:
            self.write(os.path.join("star_hg38_150", sentinel), "x")
        found = align_rna.default_star_index(self.tmp, None)
        self.assertEqual(os.path.basename(found), "star_hg38_150")

    def test_a_half_built_index_is_not_offered(self):
        # Found while writing the test above: discovery accepted a
        # directory on one sentinel file while the alignment step demands
        # four, so an interrupted build could be auto-selected and then
        # rejected -- which reads as a contradiction rather than as
        # "that index is unfinished". Both now ask the same question.
        import align_rna
        os.makedirs(os.path.join(self.tmp, "star_hg38_100"))
        self.write(os.path.join("star_hg38_100", "SAindex"), "x")
        self.assertEqual(align_rna.installed_star_indexes(self.tmp), [])

    def test_an_explicit_read_length_still_names_its_own_index(self):
        import align_rna
        self.assertTrue(align_rna.default_star_index(self.tmp, 100)
                        .endswith("star_hg38_100"))


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestManifestResolution(TempCase):
    """
    BUGS in the report step's two lookups, both of which failed by
    reporting something true and misleading rather than by crashing.
    """

    def test_vcf_preference_yields_to_existence(self):
        # BUG: the COSMIC-annotated VCF is preferred, but preference is not
        # availability. Naming it and testing only it meant a run whose
        # COSMIC step was skipped reported "no VCF on disk" while its
        # filtered calls sat right there.
        from jobs import manifest_vcf
        filtered = self.touch("filtered.vcf.gz")
        self.assertEqual(
            manifest_vcf({"cosmic_vcf": "/gone/cosmic.vcf.gz",
                          "filtered_vcf": filtered}), filtered)

    def test_preference_is_honoured_when_both_exist(self):
        from jobs import manifest_vcf
        cosmic = self.touch("cosmic.vcf.gz")
        filtered = self.touch("filtered.vcf.gz")
        self.assertEqual(manifest_vcf({"cosmic_vcf": cosmic,
                                       "filtered_vcf": filtered}), cosmic)

    def test_manifest_lookup_follows_the_assay(self):
        # BUG: the glob was DNA-only, so an RNA run's report step reported
        # "no run manifest found" -- true of the pattern it searched for,
        # and entirely misleading about the run.
        from jobs import latest_manifest
        self.write("rna_manifest_1.json", '{"assay": "rna"}')
        self.assertIsNone(latest_manifest(self.tmp, assay="dna"))
        self.assertEqual(latest_manifest(self.tmp, assay="rna"),
                         {"assay": "rna"})


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestArtefactIdentityIsStable(TempCase):
    """
    BUG: report links addressed artefacts by list position.

    The list grows while a run is in flight. A link captured at step 3
    pointed at a different file by step 6 -- the pathologist clicked
    "coverage report" and was shown someone else's artefact. Wrong data
    presented confidently is the worst failure this system can produce,
    so identity is now content-addressed.
    """

    def test_id_depends_on_the_path_not_on_position(self):
        from artefacts import artefact_id
        first = artefact_id(os.path.join(self.tmp, "a", "report.html"),
                            self.tmp, self.tmp)
        self.write("b.txt", "x")
        again = artefact_id(os.path.join(self.tmp, "a", "report.html"),
                            self.tmp, self.tmp)
        self.assertEqual(first, again)

    def test_different_files_get_different_ids(self):
        from artefacts import artefact_id
        self.assertNotEqual(
            artefact_id(os.path.join(self.tmp, "a.html"), self.tmp, self.tmp),
            artefact_id(os.path.join(self.tmp, "b.html"), self.tmp, self.tmp))


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestJobReleasesItsPipe(TempCase):
    """
    BUG: a finished job held its subprocess pipe open.

    The run loop read the pipe to EOF but never closed it, and CPython
    releases the read end only when the Popen object is collected. The
    manager deliberately keeps every Job for the life of the server, so
    finished runs stay listable -- which meant each completed run also
    kept a file descriptor, and a server working through a long worksheet
    accumulated them until it ran out.

    Not a wrong answer, so it belongs in a different class of bug from the
    rest of this module: it degrades a long-lived server rather than
    misreporting a specimen. It is pinned here because it is invisible
    until the limit is hit, and then presents as something else entirely.
    """

    def test_the_pipe_is_closed_once_the_run_is_over(self):
        import jobs
        run_dir = os.path.join(self.tmp, "run")
        os.makedirs(run_dir)
        job = jobs.Job("t1", run_dir, sys.executable,
                       [sys.executable, "-c", "print('hello')"],
                       {"patient_id": "P"}, {})
        job.run()
        self.assertEqual(job.returncode, 0)
        self.assertTrue(job._proc.stdout.closed,
                        "the job finished still holding its pipe open")

    def test_output_still_reaches_the_log(self):
        # Closing the pipe must not cost the log its contents -- the log is
        # the only record of what a run actually did.
        import jobs
        run_dir = os.path.join(self.tmp, "run2")
        os.makedirs(run_dir)
        job = jobs.Job("t2", run_dir, sys.executable,
                       [sys.executable, "-c", "print('marker-line')"],
                       {"patient_id": "P"}, {})
        job.run()
        with open(job.log_path, encoding="utf-8") as fh:
            self.assertIn("marker-line", fh.read())


@unittest.skipUnless(have_flask(), "Flask is not installed")
class TestEveryHandlerAcceptsValidateOnly(unittest.TestCase):
    """
    BUG: submit_hybrid() had no validate_only parameter.

    The worksheet validates each row before queueing any of them. The
    missing parameter raised TypeError during QUEUEING -- after earlier
    rows had already started -- so a worksheet failed halfway through,
    which is the worst possible moment for a batch to fail.
    """

    def test_all_three_handlers_expose_the_parameter(self):
        import inspect
        import app
        for name in ("submit", "submit_rna", "submit_hybrid"):
            handler = getattr(app, name)
            params = inspect.signature(handler).parameters
            self.assertIn("validate_only", params,
                          f"{name}() cannot be dry-run by the worksheet")


if __name__ == "__main__":
    unittest.main()
