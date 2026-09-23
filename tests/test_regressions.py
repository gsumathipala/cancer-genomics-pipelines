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

from helpers import TempCase, have_flask
from webcase import WebCase
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
    GUARD: STAR's chimeric count and the BAM must agree.

    The fusion caller reads only the BAM. If STAR counted chimeric reads
    but none reached it, the run would produce a clean, empty, well-formed
    fusion table -- identical in every downstream artefact to a true
    negative. These cases pin the guard's LOGIC with the detector mocked;
    TestChimericDetectorOnARealBam pins the detector itself, which is the
    half that was wrong.
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


def _find_samtools():
    """samtools on PATH, or in any conda env under the usual roots."""
    import glob
    import shutil
    found = shutil.which("samtools")
    if found:
        return os.path.dirname(found)
    for root in ("~/miniconda3", "~/anaconda3", "~/miniforge3"):
        hits = sorted(glob.glob(os.path.expanduser(
            f"{root}/envs/*/bin/samtools")))
        if hits:
            return os.path.dirname(hits[0])
    return None


@unittest.skipUnless(_find_samtools(), "samtools not installed anywhere")
class TestChimericDetectorOnARealBam(TempCase):
    """
    BUG: the chimeric-output guard fired on EVERY run with chimeric reads.

    It looked for STAR's ch:A:1 tag, which STAR writes only when
    --outSAMattributes asks for it -- and this pipeline never did. So it
    reported "no chimeric alignments in the BAM" whenever there were any.
    Found by a real run on simulated reads: the error fired while Arriba,
    reading that same BAM, called EML4::ALK and BCR::ABL1 at their exact
    breakpoints. The previous tests mocked the detector, which is how it
    got through; these build real BAMs and run real samtools.
    """

    HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr2\tLN:50000000\n"
    # A chimeric read the way STAR WithinBAM writes it: a primary record
    # and a SUPPLEMENTARY one (flag 0x800), linked by SA tags -- and, as
    # in every BAM this pipeline makes, no ch tag.
    CHIMERIC = (
        "r1\t65\tchr2\t100\t255\t50M50S\t=\t500\t0\t{seq}\t{q}"
        "\tSA:Z:chr2,40000,+,50S50M,255,0;\n"
        "r1\t2113\tchr2\t40000\t255\t50H50M\t=\t500\t0\t{half}\t{hq}"
        "\tSA:Z:chr2,100,+,50M50S,255,0;\n")
    ORDINARY = "r2\t0\tchr2\t2000\t255\t100M\t*\t0\t0\t{seq}\t{q}\n"

    def setUp(self):
        super().setUp()
        self._path = os.environ.get("PATH", "")
        os.environ["PATH"] = _find_samtools() + os.pathsep + self._path

    def tearDown(self):
        os.environ["PATH"] = self._path
        super().tearDown()

    def _bam(self, body):
        import subprocess
        fill = dict(seq="A" * 100, q="I" * 100, half="A" * 50, hq="I" * 50)
        sam = self.write("in.sam", self.HEADER + body.format(**fill))
        bam = os.path.join(self.tmp, "out.bam")
        subprocess.run(["samtools", "view", "-b", "-o", bam, sam],
                       check=True, capture_output=True)
        return bam

    def test_supplementary_chimeric_records_are_found(self):
        import align_rna
        bam = self._bam(self.ORDINARY + self.CHIMERIC)
        self.assertTrue(align_rna.bam_has_chimeric_alignments(bam))

    def test_a_bam_without_them_is_reported_as_such(self):
        import align_rna
        bam = self._bam(self.ORDINARY)
        self.assertFalse(align_rna.bam_has_chimeric_alignments(bam))

    def test_the_guard_is_silent_when_the_bam_agrees(self):
        # The exact false alarm: STAR counted chimeras, they ARE in the
        # BAM, and the guard must say nothing.
        import align_rna
        bam = self._bam(self.ORDINARY + self.CHIMERIC)
        self.assertEqual(align_rna.verify_chimeric_output(
            {"bam": bam, "prefix": os.path.join(self.tmp, "")},
            {"Number of chimeric reads": 225}), [])


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


def _function_source(module_file, name):
    """The exact source of a top-level function, for copy comparison."""
    import ast
    path = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), module_file)
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(text, node)
    raise AssertionError(f"{name} not found in {module_file}")


class TestTumourAndNormalAreNeverGuessed(unittest.TestCase):
    """
    BUG: with --auto-discover the tumour was the alphabetically FIRST sample.

    The engine read "assign the first pair as tumour, second as normal" and
    ignored --tumour-sample altogether. Ordinary names put the normal first
    -- PT01-N before PT01-T, "normal" before "tumour" -- so the NORMAL was
    analysed as the tumour. Mutect2 then found every real somatic mutation
    in the "normal" and filtered it as germline, and the run finished with
    no error and a short, plausible, wrong list. Found by running the real
    pipeline on simulated reads with planted BRAF/EGFR/KRAS/TP53/PIK3CA
    mutations, where it printed "Tumour : NORMAL_S2".

    Two smaller bugs sat beside it: --normal-sample could not be given at
    all without --normal-r1/r2, so in the discovery modes the normal could
    not be named; and a tumour name that was NOT found only warned, then
    fell back to the first sample -- a typo became a swapped run.
    """

    def setUp(self):
        import align_reads
        import comprehensive_variant_calling
        self.copies = (comprehensive_variant_calling.resolve_roles,
                       align_reads.resolve_roles)

    def test_both_copies_are_identical(self):
        # The scripts are self-contained by design, so the rule is written
        # twice. Two copies of a safety rule are only safe while identical.
        self.assertEqual(
            _function_source("comprehensive_variant_calling.py",
                             "resolve_roles"),
            _function_source("align_reads.py", "resolve_roles"))

    def _each(self, *args, **kwargs):
        return [fn(*args, **kwargs) for fn in self.copies]

    def _refused(self, *args, **kwargs):
        for fn in self.copies:
            with self.assertRaises(ValueError):
                fn(*args, **kwargs)

    def test_the_normal_is_not_taken_from_alphabetical_order(self):
        # The exact case that swapped: the normal sorts first.
        for result in self._each({"PT01-N", "PT01-T"}, "PT01-T", "PT01-N"):
            self.assertEqual(result, ("PT01-T", "PT01-N"))

    def test_several_samples_and_no_tumour_named_is_refused(self):
        self._refused({"PT01-N", "PT01-T"})

    def test_the_normal_is_never_inferred_from_what_is_left(self):
        # Two samples in a folder may be two different people.
        self._refused({"PT01-N", "PT01-T"}, "PT01-T")

    def test_a_name_that_is_not_there_is_refused_not_replaced(self):
        self._refused({"PT01-N", "PT01-T"}, "PT01-X", "PT01-N")
        self._refused({"PT01-N", "PT01-T"}, "PT01-T", "PT01-Q")

    def test_one_sample_and_no_names_is_tumour_only(self):
        for result in self._each({"ONLY"}):
            self.assertEqual(result, ("ONLY", None))

    def test_a_declared_tumour_only_run_ignores_the_neighbours(self):
        for result in self._each({"A", "B", "C"}, "B", tumour_only=True):
            self.assertEqual(result, ("B", None))

    def test_contradictions_are_refused(self):
        self._refused({"A", "B"}, "A", "A")                  # same sample
        self._refused({"A", "B"}, None, "B")                 # normal, no tumour
        self._refused({"A", "B"}, "A", "B", tumour_only=True)

    def test_the_engine_accepts_a_named_normal_when_discovering(self):
        # Used to be rejected unless --normal-r1/r2 were given too, which
        # made naming the normal impossible in the discovery modes.
        import comprehensive_variant_calling as cvc
        parser = cvc.build_parser()
        with __import__("tempfile").TemporaryDirectory() as d:
            args = parser.parse_args(
                ["-i", d, "-o", d, "-r", "hg38", "--auto-discover",
                 "--tumour-sample", "PT01-T", "--normal-sample", "PT01-N"])
            cvc.validate_args(args, parser)          # must not exit


class TestEveryLaneIsUsed(TempCase):
    """
    BUG: a lane-split sample was analysed on lane 1 only.

    Discovery reports such a sample with lane 1 in "r1" and every lane in
    "lanes". Three places used "r1" alone: the web interface's batch mode
    and its worksheet scan (copied into an explicit --tumour-r1), and the
    RNA engine's auto-discovery, which had no lane merging at all. A
    four-lane library ran on a quarter of its reads: a quarter of the
    depth, low-fraction variants and fusions lost, and -- on the RNA side
    -- a library QC that then blamed the specimen for being shallow.
    """

    def _lanes(self, sample="PT01_S1", lanes=4):
        import gzip
        for lane in range(1, lanes + 1):
            for read in (1, 2):
                path = os.path.join(
                    self.tmp, f"{sample}_L00{lane}_R{read}_001.fastq.gz")
                with gzip.open(path, "wt") as fh:
                    fh.write(f"@r{lane}\nACGT\n+\nIIII\n")
        return os.path.join(self.tmp, f"{sample}_L001_R1_001.fastq.gz")

    def test_a_lane_is_recognised_as_part_of_its_sample(self):
        import comprehensive_variant_calling as cvc
        lane1 = self._lanes()
        self.assertEqual(cvc.lane_split_owner(lane1), ("PT01_S1", 4))

    def test_a_standalone_file_is_not_mistaken_for_a_lane(self):
        import comprehensive_variant_calling as cvc
        alone = self.touch("SOLO_S1_R1_001.fastq.gz")
        self.touch("SOLO_S1_R2_001.fastq.gz")
        self.assertIsNone(cvc.lane_split_owner(alone))

    def test_the_engine_refuses_one_lane_given_explicitly(self):
        import comprehensive_variant_calling as cvc
        lane1 = self._lanes()
        parser = cvc.build_parser()
        args = parser.parse_args(
            ["-i", self.tmp, "-o", self.tmp, "-r", "hg38",
             "--tumour-sample", "PT01_S1", "--tumour-r1", lane1,
             "--tumour-r2", lane1.replace("_R1_", "_R2_")])
        with self.assertRaises(SystemExit):
            with __import__("contextlib").redirect_stderr(
                    __import__("io").StringIO()):
                cvc.validate_args(args, parser)

    def test_the_rna_engine_merges_every_lane(self):
        import gzip
        import fusion_calling
        import comprehensive_variant_calling as cvc
        self._lanes()
        pairs, _ = cvc.find_paired_fastqs(self.tmp, "*_R1_*.fastq*",
                                          "*_R2_*.fastq*")
        merged, written = fusion_calling.merge_lanes(
            pairs[0], os.path.join(self.tmp, "merge"))
        with gzip.open(merged["r1"], "rt") as fh:
            reads = [line for line in fh if line.startswith("@")]
        # Every lane's read is present, not just lane 1's.
        self.assertEqual(sorted(reads), ["@r1\n", "@r2\n", "@r3\n", "@r4\n"])
        self.assertEqual(len(written), 2)


class TestTheWebInterfaceNamesTheSamples(WebCase):
    """The same two bugs, as the web interface met them."""

    def _dna_form(self, **extra):
        form = {
            "output_dir": os.path.join(self.runs, "OUT"),
            "reference": self.reference(),
            "panel_bed": self.panel_bed(),
            "input_dir": os.path.join(self.tmp, "fq"),
            "patient_id": "MRN-1", "threads": "2", "dry_run": "on",
            "form_rendered": "1", "minimal_run": "on",
        }
        form.update(extra)
        return form

    def _queued_argv(self):
        jobs = self.app_module.manager.all()
        return jobs[-1].argv if jobs else None

    def test_the_pair_checkbox_needs_the_tumour_named(self):
        self.fastqs("PT01-N", "PT01-T")
        response = self.client.post("/submit", data=self._dna_form(
            auto_discover="on", second_is_matched_normal="on"))
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"TUMOUR", response.data)
        self.assertEqual(self.job_count(), 0)

    def test_the_pair_checkbox_sends_both_names(self):
        # The normal sorts first; the tumour must still be the tumour.
        self.fastqs("PT01-N", "PT01-T")
        response = self.client.post("/submit", data=self._dna_form(
            auto_discover="on", second_is_matched_normal="on",
            tumour_sample="PT01-T"))
        self.assertEqual(response.status_code, 302, response.data[:400])
        argv = self._queued_argv()
        self.assertEqual(argv[argv.index("--tumour-sample") + 1], "PT01-T")
        self.assertEqual(argv[argv.index("--normal-sample") + 1], "PT01-N")

    def test_batch_runs_each_sample_by_name(self):
        # By name with auto-discover kept on, so the engine merges lanes;
        # never pair["r1"], which is lane 1 of a lane-split sample.
        self.fastqs("S1", "S2")
        response = self.client.post("/submit", data=self._dna_form(
            auto_discover="on", batch_mode="on"))
        self.assertEqual(response.status_code, 302, response.data[:400])
        for job in self.app_module.manager.all():
            self.assertIn("--auto-discover", job.argv)
            self.assertIn("--tumour-only", job.argv)
            self.assertNotIn("--tumour-r1", job.argv)

    def test_a_worksheet_row_on_one_lane_runs_the_whole_sample(self):
        import gzip
        for lane in (1, 2):
            for read in (1, 2):
                path = os.path.join(
                    self.tmp, f"LS_S1_L00{lane}_R{read}_001.fastq.gz")
                with gzip.open(path, "wt") as fh:
                    fh.write("@r\nACGT\n+\nIIII\n")
        values = self.app_module.worksheet_row_to_form(
            {"sample": "my label", "assay": "dna",
             "r1": os.path.join(self.tmp, "LS_S1_L001_R1_001.fastq.gz"),
             "r2": os.path.join(self.tmp, "LS_S1_L001_R2_001.fastq.gz")},
            {"output_dir": "/o"})
        self.assertEqual(values.get("auto_discover"), "1")
        self.assertEqual(values["tumour_sample"], "LS_S1")
        self.assertEqual(values.get("tumour_only"), "1")
        self.assertNotIn("tumour_r1", values)

    def test_one_lane_typed_into_the_form_is_refused(self):
        import gzip
        for lane in (1, 2):
            for read in (1, 2):
                path = os.path.join(
                    self.tmp, "fq", f"LS_S1_L00{lane}_R{read}_001.fastq.gz")
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with gzip.open(path, "wt") as fh:
                    fh.write("@r\nACGT\n+\nIIII\n")
        lane1 = os.path.join(self.tmp, "fq", "LS_S1_L001_R1_001.fastq.gz")
        response = self.client.post("/submit", data=self._dna_form(
            tumour_sample="LS_S1", tumour_r1=lane1,
            tumour_r2=lane1.replace("_R1_", "_R2_")))
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"lanes of", response.data)
        self.assertEqual(self.job_count(), 0)


class TestThePanelBedIsNeverDropped(WebCase):
    """
    BUG: the panel BED was used only when a vendor profile was ALSO chosen.

    The form's panel BED fills both the calling intervals and the coverage
    BED -- but that copy sat after a "no profile selected" early return. A
    laboratory running its own panel filled in the BED and had it silently
    dropped: Mutect2 called genome-wide, step 12 made no coverage statement
    ("absent variants cannot be told apart from unsequenced regions"), and
    the TMB denominator was never measured. Found by submitting a real run
    through the web form, with no profile, and reading step 12's log.
    """

    def test_a_custom_panel_bed_reaches_the_pipeline(self):
        fq = self.fastqs("T1")
        bed = self.panel_bed()
        response = self.client.post("/submit", data={
            "output_dir": os.path.join(self.runs, "OUT"),
            "reference": self.reference(), "panel_bed": bed,
            "input_dir": os.path.dirname(fq["T1"][0]),
            "tumour_sample": "T1", "tumour_r1": fq["T1"][0],
            "tumour_r2": fq["T1"][1], "patient_id": "MRN-1",
            "threads": "2", "dry_run": "on", "form_rendered": "1",
            "minimal_run": "on"})               # note: no "panel" at all
        self.assertEqual(response.status_code, 302, response.data[:300])
        argv = self.app_module.manager.all()[-1].argv
        self.assertEqual(argv[argv.index("--intervals") + 1],
                         os.path.realpath(bed))
        self.assertEqual(argv[argv.index("--coverage-bed") + 1],
                         os.path.realpath(bed))


class TestAnRnaRunIsOneSpecimen(WebCase):
    """
    BUG: an RNA run could hold several patients' specimens.

    Auto-discovery (or a manifest) that found several libraries ran them
    ALL in one job under the one patient typed on the form, and sent only
    the alphabetically first to PCGR -- the other specimens' fusions filed
    against the wrong person. The DNA form had long refused this; the RNA
    form never checked.
    """

    def _rna_form(self, **extra):
        form = {
            "output_dir": os.path.join(self.runs, "RNA"),
            "reference": self.reference(), "gtf": self.gtf(),
            "star_index": self.star_index(),
            "input_dir": os.path.join(self.tmp, "fq"),
            "patient_id": "MRN-1", "threads": "2", "dry_run": "on",
            "form_rendered": "1", "auto_discover": "on"}
        form.update(extra)
        return form

    def test_several_unnamed_libraries_are_refused(self):
        self.fastqs("LIB_A", "LIB_B")
        response = self.client.post("/submit/rna", data=self._rna_form())
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"libraries found", response.data)
        self.assertEqual(self.job_count(), 0)

    def test_a_named_library_runs_alone(self):
        self.fastqs("LIB_A", "LIB_B")
        response = self.client.post("/submit/rna",
                                    data=self._rna_form(sample="LIB_B"))
        self.assertEqual(response.status_code, 302, response.data[:300])
        argv = self.app_module.manager.all()[-1].argv
        self.assertEqual(argv[argv.index("--sample") + 1], "LIB_B")

    def test_the_report_is_for_the_named_sample_not_the_first(self):
        from jobs import build_rna_pcgr_argv
        tsv_a, tsv_b = self.touch("a.tsv"), self.touch("b.tsv")
        manifest = {"samples": {
            "LIB_A": {"fusion_report": {"pcgr_tsv": tsv_a}},
            "LIB_B": {"fusion_report": {"pcgr_tsv": tsv_b}}}}
        argv, _ = build_rna_pcgr_argv(
            "py", ".", manifest,
            {"pcgr_refdata_dir": "/rd", "sample": "LIB_B"}, "/o")
        self.assertEqual(argv[argv.index("--sample-id") + 1], "LIB_B")
        # And with several and none named, no report is attributed at all.
        argv, note = build_rna_pcgr_argv(
            "py", ".", manifest, {"pcgr_refdata_dir": "/rd"}, "/o")
        self.assertIsNone(argv)
        self.assertIn("none named", note)


class TestQueuedRunsSurviveARestart(TempCase):
    """
    BUG: a run still QUEUED when the server restarted vanished.

    The run list is rebuilt from each run's state file, and a queued job
    had none -- it was first written when the job started. After a restart
    (which install.md prescribes after every --update) every queued
    worksheet row but the running one was simply gone, with no error and
    no trace in the list. load_existing() even had a branch to mark queued
    jobs "interrupted"; it could never fire.
    """

    def test_a_queued_job_is_listed_as_interrupted_after_a_restart(self):
        import threading
        import jobs
        runs = os.path.join(self.tmp, "runs")
        manager = jobs.JobManager(runs)
        # A worker that is "busy" forever, so the job stays queued -- the
        # state it would be in when the server went down.
        blocker = threading.Event()
        manager._worker = threading.Thread(target=blocker.wait, daemon=True)
        manager._worker.start()
        try:
            job = manager.submit("x.py", ["py", "x.py"],
                                 {"patient_id": "P"}, meta={}, assay="dna")
            self.assertEqual(job.status, "queued")

            restarted = jobs.JobManager(runs)
            restarted.load_existing()
            found = restarted.get(job.id)
            self.assertIsNotNone(found, "the queued run vanished")
            self.assertEqual(found.status, "interrupted")
        finally:
            blocker.set()


class TestAnUnjudgedLibraryIsNotCertified(unittest.TestCase):
    """
    BUG: with no depth floor, RNA QC certified any library as "negative-capable".

    The depth and mapping floors are assay-specific, so they come from a
    panel profile or the operator. With neither, those two checks were
    silently LEFT OUT and the verdict read "every measured check cleared its
    bar -- an empty fusion table is a negative result". Found on a real
    web run of a 0.6-million-read library, certified exactly so.
    """

    MEASURED = {"input_reads_millions": 0.6, "uniquely_mapped_pct": 95.1,
                "too_short_pct": 0.1, "splice_junctions_total": 400000,
                "chimeric_reads": 225}

    def test_no_floor_means_not_judged_and_no_negative(self):
        import rna_qc_report as r
        checks = r.evaluate(self.MEASURED, {"min_reads_millions": None,
                                            "min_unique_mapped_pct": None})
        states = {c["name"]: c["state"] for c in checks}
        self.assertEqual(states["Library size"], "unassessed")
        self.assertEqual(states["Uniquely mapped"], "unassessed")
        v = r.verdict(checks)
        self.assertEqual(v["state"], "partial")
        self.assertIn("cannot yet be reported as a negative", v["text"])

    def test_with_a_floor_a_shallow_library_is_a_concern(self):
        import rna_qc_report as r
        checks = r.evaluate(self.MEASURED, {"min_reads_millions": 20,
                                            "min_unique_mapped_pct": 60})
        self.assertEqual(r.verdict(checks)["state"], "concern")


class TestTheOrchestratorGivesEveryStageTheSameNames(unittest.TestCase):
    """
    BUG: the orchestrator's stages could disagree about which sample is which.

    Stage 2 and stage 3 each resolve tumour and normal for themselves, and
    the orchestrator's own documented example named neither -- so each
    guessed, alphabetically, getting PT01-N/PT01-T backwards. Names are now
    given once, at the top, and forwarded; names given the old way (stage 3
    only) are copied to stage 2; and names given both ways are refused.
    """

    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def _dry(self, *extra):
        import subprocess
        import sys
        return subprocess.run(
            [sys.executable, os.path.join(self.ROOT, "pipeline_orchestrator.py"),
             "--dry-run", "--stage1-args", "-i /x -o /q"] + list(extra),
            capture_output=True, text=True)

    def _stage_lines(self, out):
        return [l for l in out.splitlines()
                if "align_reads.py" in l or
                "comprehensive_variant_calling.py" in l]

    def test_top_level_names_reach_both_stages(self):
        run = self._dry("--tumour-sample", "PT01-T", "--normal-sample",
                        "PT01-N", "--stage2-args", "-o /a --reference hg38",
                        "--stage3-args", "-o /r --reference hg38")
        lines = self._stage_lines(run.stdout)
        self.assertEqual(len(lines), 2, run.stdout[-800:])
        for line in lines:
            self.assertIn("--tumour-sample PT01-T", line)
            self.assertIn("--normal-sample PT01-N", line)

    def test_names_given_to_stage_3_only_are_copied_to_stage_2(self):
        run = self._dry("--stage2-args", "-o /a --reference hg38",
                        "--stage3-args", "-o /r --reference hg38 "
                        "--tumour-sample PT01-T --normal-sample PT01-N")
        for line in self._stage_lines(run.stdout):
            self.assertIn("--tumour-sample PT01-T", line)
            self.assertIn("--normal-sample PT01-N", line)

    def test_names_given_both_ways_are_refused(self):
        run = self._dry("--tumour-sample", "PT01-T",
                        "--stage3-args", "-o /r --tumour-sample PT01-X")
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("Give them once", run.stderr)

    def test_an_rna_run_has_no_normal(self):
        run = self._dry("--assay", "rna", "--normal-sample", "PT01-N",
                        "--stage3-args", "-o /r")
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("no matched normal", run.stderr)


class TestTheInstallerComparesIndexAndAnnotation(TempCase):
    """
    BUG: the installer reported an index/annotation mismatch as all green.

    The default GENCODE release moved from v44 to v50. The v50 GTF was
    downloaded; the STAR-index step saw the index's sentinel and skipped;
    --check reported "GENCODE v50 annotation ok" and "STAR index ok" -- a
    pair it never compared, for an index built from v44.
    """

    def _index(self, built_from):
        import json
        index = os.path.join(self.tmp, "star_hg38_150")
        os.makedirs(index)
        if built_from:
            self.write(os.path.join(index, "pipeline_index.json"),
                       json.dumps({"gtf": f"/d/{built_from}"}))
        return index

    def test_a_different_release_is_named(self):
        import install_pipeline
        index = self._index("gencode.v44.primary_assembly.annotation.gtf")
        note = install_pipeline.annotation_mismatch(
            index, "/d/gencode.v50.primary_assembly.annotation.gtf")
        self.assertIn("gencode.v44", note)
        self.assertIn("--only star-index --force", note)

    def test_a_matching_release_is_silent(self):
        import install_pipeline
        index = self._index("gencode.v50.primary_assembly.annotation.gtf")
        self.assertIsNone(install_pipeline.annotation_mismatch(
            index, "/x/gencode.v50.primary_assembly.annotation.gtf"))

    def test_no_record_is_not_invented_into_a_mismatch(self):
        import install_pipeline
        self.assertIsNone(install_pipeline.annotation_mismatch(
            self._index(None), "/d/gencode.v50.x.gtf"))


class TestStep13FindsPcgrsOwnEnvironment(TempCase):
    """
    BUG: step 13 of the command-line pipeline could never produce a report.

    PCGR is always installed in its own conda environment and the pipeline
    runs in another, so `pcgr` is never on PATH there -- and step 13 checked
    PATH, printed "skipping", and --pcgr-refdata-dir quietly did nothing on
    every standard install. It now finds the environment and runs PCGR
    activated in it (CONDA_PREFIX set, which is what PCGR needs to find its
    VEP plugins). Found by a real CLI run whose step 13 skipped.
    """

    def setUp(self):
        super().setUp()
        self._env = dict(os.environ)
        root = os.path.join(self.tmp, "conda")
        self.pcgr_env = os.path.join(root, "envs", "pcgr")
        self.write(os.path.join(self.pcgr_env, "bin", "pcgr"), "#!/bin/sh\n")
        os.environ["CONDA_ROOT"] = root
        # No pcgr on PATH -- the situation step 13 is always in.
        os.environ["PATH"] = os.path.join(self.tmp, "empty")
        # The search looks beside the running interpreter and under ~
        # FIRST (as the web app's does), which on a real machine finds the
        # real env. Point both into scratch so only CONDA_ROOT can match.
        import sys
        self._exe = sys.executable
        sys.executable = os.path.join(self.tmp, "py", "bin", "python3")
        os.environ["HOME"] = self.tmp

    def tearDown(self):
        import sys
        sys.executable = self._exe
        os.environ.clear()
        os.environ.update(self._env)
        super().tearDown()

    def test_the_environment_is_found(self):
        import pcgr_report
        self.assertEqual(pcgr_report.locate_pcgr_env(), self.pcgr_env)

    def test_pcgr_runs_activated_in_it(self):
        import contextlib
        import io as _io
        import pcgr_report
        out = _io.StringIO()
        with contextlib.redirect_stdout(out):
            result = pcgr_report.run_pcgr(
                "/in.vcf", self.tmp, "S1", "/refdata", vep_dir="/vep",
                dry_run=True)
        printed = out.getvalue()
        self.assertTrue(result)                # not the old "skipped" None
        self.assertIn(os.path.join(self.pcgr_env, "bin", "pcgr"), printed)
        self.assertEqual(pcgr_report.pcgr_environment(self.pcgr_env)
                         ["CONDA_PREFIX"], self.pcgr_env)


class TestTheRnaFormOffersAMatchingAnnotation(WebCase):
    """
    BUG: the RNA form offered the newest GTF, not the index's own.

    Found on the installed machine: the STAR index was built from GENCODE
    v44, a v50 GTF had been downloaded beside it, and the form pre-filled
    v50 -- so every RNA run paired STAR's junctions from one release with
    Arriba's annotation from another. Two sort bugs sat in the same
    function: versions and read lengths were compared as TEXT, so
    "gencode.v9" beat "gencode.v44" and "star_hg38_75" beat
    "star_hg38_150"; and an index was offered on one sentinel file.
    """

    def _data(self):
        import json
        data = os.path.join(self.tmp, "data")
        gencode = os.path.join(data, "references", "gencode")
        os.makedirs(gencode)
        for release in ("9", "44", "50"):
            self.write(os.path.join(
                gencode, f"gencode.v{release}.primary_assembly."
                         f"annotation.gtf"), "x")
        for length in ("75", "150"):
            for name in ("SA", "SAindex", "Genome", "genomeParameters.txt"):
                self.write(os.path.join(data, "references",
                                        f"star_hg38_{length}", name), "x")
        # A half-built index: one sentinel, and the longest read length.
        self.write(os.path.join(data, "references", "star_hg38_250",
                                "SAindex"), "x")
        return data, gencode, json

    def test_the_index_s_own_annotation_is_offered(self):
        data, gencode, json = self._data()
        built = os.path.join(gencode,
                             "gencode.v44.primary_assembly.annotation.gtf")
        self.write(os.path.join(data, "references", "star_hg38_150",
                                "pipeline_index.json"),
                   json.dumps({"gtf": built, "read_length": 150}))
        found = self.app_module.default_rna_resources(data)
        self.assertEqual(found["gtf"], built)          # not the newer v50

    def test_without_a_record_versions_compare_as_numbers(self):
        data, _gencode, _ = self._data()
        found = self.app_module.default_rna_resources(data)
        self.assertIn(".v50.", found["gtf"])           # not v9

    def test_read_lengths_compare_as_numbers_and_half_built_is_skipped(self):
        data, _gencode, _ = self._data()
        found = self.app_module.default_rna_resources(data)
        self.assertTrue(found["star_index"].endswith("star_hg38_150"))
        self.assertEqual(found["read_length"], "150")


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
