# Created by Brainstorm, 2026.
"""
Panel profiles: the registry, the precedence rule, and the BED arithmetic.

The precedence rule -- explicit flag > profile > default -- is the whole
contract of the profile system. If it inverts, a run is configured by
something the operator cannot see, which is the failure the system exists
to prevent.
"""

import argparse
import unittest

from helpers import TempCase, bed, fai, read_text
import panel_profiles as pp


class TestRegistry(unittest.TestCase):

    def test_builtins_all_validate(self):
        registry = pp.available_panels()
        self.assertGreater(len(registry), 15)
        for profile in registry.values():
            self.assertIn(profile["chemistry"], pp.CHEMISTRIES)
            self.assertTrue(profile["id"])
            self.assertTrue(profile["name"])

    def test_assay_type_follows_chemistry(self):
        # Derived, never declared -- the two can then never disagree.
        self.assertEqual(pp.assay_type({"chemistry": "rna-capture"}), "rna")
        self.assertEqual(pp.assay_type({"chemistry": "amplicon"}), "dna")
        self.assertEqual(pp.assay_type({}), "dna")

    def test_resolve_by_id_and_alias(self):
        registry = pp.available_panels()
        by_id = pp.resolve_panel("illumina-tso500", registry)
        by_alias = pp.resolve_panel("tso500", registry)
        self.assertEqual(by_id["id"], by_alias["id"])

    def test_unknown_panel_raises_with_a_suggestion(self):
        registry = pp.available_panels()
        with self.assertRaises(pp.PanelError) as caught:
            pp.resolve_panel("illumina-tso5000", registry)
        self.assertIn("Did you mean", str(caught.exception))

    def test_every_rna_profile_is_routed_to_the_rna_branch(self):
        for profile in pp.available_panels().values():
            if profile["chemistry"].startswith("rna-"):
                self.assertEqual(pp.assay_type(profile), "rna", profile["id"])

    def test_hybrid_partners_are_mutual_and_opposite(self):
        registry = pp.available_panels()
        for profile in registry.values():
            partner_id = pp.partner_id(profile)
            if not partner_id:
                continue
            partner = registry[partner_id]
            self.assertEqual(pp.partner_id(partner), profile["id"])
            self.assertNotEqual(pp.assay_type(profile),
                                pp.assay_type(partner))


class TestProfileValidation(TempCase):

    def _load(self, body):
        return pp.load_panel_file(self.write("p.json", body))

    def test_unknown_setting_key_is_an_error(self):
        # A typo that silently did nothing is the failure this prevents.
        with self.assertRaises(pp.PanelError) as caught:
            self._load('{"id":"x","name":"X","manufacturer":"Y",'
                       '"chemistry":"amplicon","settings":{"min_dept":10}}')
        self.assertIn("min_depth", str(caught.exception))

    def test_unknown_chemistry_is_an_error(self):
        with self.assertRaises(pp.PanelError):
            self._load('{"id":"x","name":"X","manufacturer":"Y",'
                       '"chemistry":"capture"}')

    def test_user_profile_overriding_a_builtin_inherits_its_aliases(self):
        path = self.write("tso.json",
                          '{"id":"illumina-tso500","name":"Local",'
                          '"manufacturer":"Illumina",'
                          '"chemistry":"hybrid-capture"}')
        registry = pp.available_panels(extra_files=[path])
        # '--panel tso500' must keep working after a site overrides the id.
        self.assertEqual(pp.resolve_panel("tso500", registry)["name"], "Local")


class TestPrecedence(unittest.TestCase):

    def _args(self, **kwargs):
        parser = argparse.ArgumentParser()
        parser.add_argument("--min-depth", type=int, default=0)
        parser.add_argument("--interval-padding", type=int, default=0)
        parser.add_argument("--skip-steps", nargs="*", default=[])
        namespace = parser.parse_args([])
        for key, value in kwargs.items():
            setattr(namespace, key, value)
        return parser, namespace

    def test_explicit_flag_beats_the_profile(self):
        parser, args = self._args()
        args.min_depth = 99
        profile = {"settings": {"min_depth": 50}}
        applied, overridden = pp.apply_profile(args, profile, {"min_depth"})
        self.assertEqual(args.min_depth, 99)
        self.assertEqual(overridden, [("min_depth", 50)])
        self.assertEqual(applied, [])

    def test_profile_fills_what_the_command_line_omitted(self):
        parser, args = self._args()
        profile = {"settings": {"min_depth": 50}}
        applied, _ = pp.apply_profile(args, profile, set())
        self.assertEqual(args.min_depth, 50)
        self.assertEqual(applied, [("min_depth", 50)])

    def test_skip_steps_is_merged_never_replaced(self):
        # A profile skipping dedup plus a command line skipping qc must
        # produce a run that skips BOTH.
        parser, args = self._args()
        args.skip_steps = ["qc"]
        profile = {"settings": {"skip_steps": ["dedup"]}}
        pp.apply_profile(args, profile, set())
        self.assertEqual(sorted(args.skip_steps), ["dedup", "qc"])

    def test_rename_maps_a_setting_onto_a_different_dest(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--min-depth", type=int, default=0)
        args = parser.parse_args([])
        profile = {"settings": {"coverage_min_depth": 100}}
        applied, _ = pp.apply_profile(
            args, profile, set(), only={"coverage_min_depth"},
            rename={"coverage_min_depth": "min_depth"})
        self.assertEqual(args.min_depth, 100)
        # Reported under the profile's own name, which is what you'd edit.
        self.assertEqual(applied, [("coverage_min_depth", 100)])

    def test_explicit_dests_reads_both_flag_spellings(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--min-depth", type=int, default=0)
        parser.add_argument("-o", "--output-dir", default=None)
        named = pp.explicit_dests(parser, ["--min-depth", "5", "-o", "/x"])
        self.assertEqual(named, {"min_depth", "output_dir"})
        self.assertEqual(pp.explicit_dests(parser, ["--min-depth=5"]),
                         {"min_depth"})


class TestBedArithmetic(TempCase):

    def test_overlapping_regions_are_merged_before_summing(self):
        # Vendor BEDs overlap their probe tiles. Summing unmerged inflates
        # the footprint and therefore deflates TMB, silently.
        path = self.write("a.bed", bed([("chr1", 100, 200),
                                        ("chr1", 150, 260)]))
        measured = pp.bed_footprint(path)
        self.assertEqual(measured[0], 160)      # 100..260, not 100+110

    def test_touching_regions_merge(self):
        path = self.write("b.bed", bed([("chr1", 0, 100), ("chr1", 100, 200)]))
        self.assertEqual(pp.bed_footprint(path)[0], 200)

    def test_junk_rows_are_ignored_not_fatal(self):
        path = self.write("c.bed",
                          "track name=x\n#comment\nchr1\tNaN\t200\n"
                          "chr1\t500\t400\nchr1\t10\t20\n")
        self.assertEqual(pp.bed_footprint(path)[0], 10)

    def test_a_file_that_is_not_a_bed_returns_none(self):
        self.assertIsNone(pp.footprint_mb(self.write("d.tsv", "a\tb\n")))

    def test_target_size_note_flags_a_large_disagreement(self):
        note = pp.target_size_note(0.5, 1.94)
        self.assertIn("differs", note)
        self.assertIsNone(pp.target_size_note(None, 1.0))
        self.assertIn("consistent", pp.target_size_note(1.9, 1.94))

    def test_tmb_advice_fires_only_below_the_floor(self):
        self.assertIsNone(pp.tmb_advice(35.0, True))
        self.assertIsNone(pp.tmb_advice(0.5, False))   # not requested
        self.assertIn("do not report", pp.tmb_advice(0.5, True))


class TestContigNaming(TempCase):

    def _reference(self, names):
        ref = self.write("ref.fa", ">x\n")
        self.write("ref.fa.fai", fai([(n, 1000) for n in names]))
        return ref

    def test_ensembl_bed_against_ucsc_reference_is_translated(self):
        ref = self._reference(["chr1", "chr7"])
        source = self.write("v.bed", bed([("1", 100, 200), ("7", 10, 20)]))
        out, note = pp.harmonise_bed_contigs(source, ref, self.path("out/x"))
        self.assertNotEqual(out, source)
        self.assertIn("renamed", note)
        with open(out) as handle:
            self.assertTrue(
                all(line.startswith("chr") for line in handle if line.strip()))

    def test_the_original_file_is_never_edited(self):
        ref = self._reference(["chr1"])
        source = self.write("v.bed", bed([("1", 100, 200)]))
        before = read_text(source)
        pp.harmonise_bed_contigs(source, ref, self.path("out/y"))
        self.assertEqual(read_text(source), before)

    def test_a_matching_bed_is_left_alone(self):
        ref = self._reference(["chr1"])
        source = self.write("m.bed", bed([("chr1", 100, 200)]))
        out, note = pp.harmonise_bed_contigs(source, ref, self.path("out/z"))
        self.assertEqual(out, source)
        self.assertIsNone(note)

    def test_contigs_absent_from_the_reference_are_dropped(self):
        ref = self._reference(["chr1"])
        source = self.write("v.bed", bed([("1", 100, 200), ("MT", 1, 50)]))
        out, note = pp.harmonise_bed_contigs(source, ref, self.path("out/w"))
        self.assertIn("dropped", note)
        self.assertNotIn("chrM", read_text(out))

    def test_a_bed_for_the_wrong_genome_says_so(self):
        ref = self._reference(["chr1"])
        source = self.write("v.bed", bed([("scaffold_9", 1, 50)]))
        out, note = pp.harmonise_bed_contigs(source, ref, self.path("out/v"))
        self.assertEqual(out, source)
        self.assertIn("wrong BED", note)


class TestSavedProfiles(unittest.TestCase):

    def test_zero_valued_settings_survive_being_saved(self):
        # `value in (None, "", [], False)` also matches 0, which silently
        # dropped --interval-padding 0 -- the value an amplicon profile
        # most needs to state.
        parser = argparse.ArgumentParser()
        parser.add_argument("--interval-padding", type=int, default=0)
        args = parser.parse_args([])
        args.interval_padding = 0
        profile = pp.profile_from_args(args, "x", chemistry="amplicon")
        self.assertEqual(profile["settings"].get("interval_padding"), 0)

    def test_paths_are_never_saved_into_a_profile(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--intervals", default=None)
        args = parser.parse_args([])
        args.intervals = "/data/panel.bed"
        profile = pp.profile_from_args(args, "x")
        self.assertNotIn("intervals", profile["settings"])

    def test_run_local_skip_steps_are_stripped(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--skip-steps", nargs="*", default=[])
        args = parser.parse_args([])
        args.skip_steps = ["qc", "dedup"]
        profile = pp.profile_from_args(args, "x", chemistry="amplicon")
        # 'qc' describes this run's inputs, not the kit.
        self.assertEqual(profile["settings"]["skip_steps"], ["dedup"])


if __name__ == "__main__":
    unittest.main()
