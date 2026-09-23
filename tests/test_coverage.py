# Created by Brainstorm, 2026.
"""
The coverage report: the module that lets a negative mean something.
"""

import unittest

from helpers import TempCase, bed
import coverage_report as cr


class TestDepthCommand(unittest.TestCase):

    def test_reports_every_position_including_uncovered(self):
        # Without -a, zero-depth positions are simply absent from the
        # output -- the exact silence this report exists to break.
        command = cr.depth_command("s.bam", "p.bed", 20, 20)
        self.assertIn("-a", command)
        self.assertEqual(command[:2], ["samtools", "depth"])

    def test_quality_floors_are_passed_through(self):
        command = cr.depth_command("s.bam", "p.bed", 30, 25)
        self.assertEqual(command[command.index("-Q") + 1], "30")
        self.assertEqual(command[command.index("-q") + 1], "25")


class TestBedParsing(TempCase):

    def test_regions_are_read_with_their_names(self):
        path = self.write("p.bed", bed([("chr1", 100, 200),
                                        ("chr2", 5, 25)]))
        regions = cr.parse_bed(path)
        self.assertEqual(len(regions), 2)

    def test_header_and_comment_lines_are_skipped(self):
        path = self.write("h.bed",
                          "track name=x\n#c\nchr1\t1\t10\tR\n")
        self.assertEqual(len(cr.parse_bed(path)), 1)


if __name__ == "__main__":
    unittest.main()
