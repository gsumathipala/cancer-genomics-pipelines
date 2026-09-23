# Created by Brainstorm, 2026.
"""
Test suite for the Cancer Genomics Pipelines.

WHY STDLIB unittest AND NOT pytest
----------------------------------
  requirements.txt promises that the analysis scripts "import nothing
  outside the Python standard library and will run against a bare Python
  3.8+ with this file never installed". A test suite that needed pytest
  would quietly withdraw that promise: a contributor could no longer verify
  the code without first installing something.

  So the suite uses unittest, which is always there:

      python3 -m unittest discover -s tests -v
      python3 run_tests.py

  The one exception is the web-interface tests, which need Flask because
  the web interface needs Flask. Those skip themselves cleanly when it is
  absent rather than failing.

WHAT THIS SUITE DOES NOT DO
---------------------------
  No test runs STAR, GATK, Arriba or PCGR. Those need the tools, about 30
  GB of reference data, and minutes to hours per run -- a suite that
  needed those is a suite nobody runs, and an unrun suite catches nothing.

  Command CONSTRUCTION is tested exhaustively; command EXECUTION is not.
  Be clear about what that means: the suite would NOT have caught the STAR
  bug where chimeric reads never reached the BAM, because only a real run
  could. It says nothing about biological correctness either -- see
  RNA_SCOPE.md and CNV_SCOPE.md for what validation still requires.
"""
