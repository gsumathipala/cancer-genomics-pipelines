#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
Run the test suite.

    python3 run_tests.py              # everything
    python3 run_tests.py regressions  # one module, by its short name
    python3 run_tests.py -v           # per-test names

Nothing needs installing. The suite uses the standard library's unittest
for the same reason the analysis scripts import nothing else: a user who
has just cloned this should be able to check that it works before they
build an environment. Tests that need Flask skip themselves when Flask is
absent, and say so rather than failing.

No sequencing tool is invoked. The suite tests the decisions -- which
flags get built, which file gets picked, what a report says about a
result -- and never the aligners and callers themselves, which are
third-party software with their own test suites and which would turn a
two-second check into an overnight one.
"""

import os
import sys
import unittest
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "tests")


def main(argv):
    verbosity = 2 if ("-v" in argv or "--verbose" in argv) else 1
    names = [a for a in argv if not a.startswith("-")]

    # tests/ is put first on the path so the modules can import both the
    # bundle's scripts and their own helpers without being a package.
    sys.path.insert(0, TESTS)
    sys.path.insert(0, HERE)

    # Leaked file handles are collected, and reported, at whatever moment
    # the garbage collector notices -- which is usually inside an unrelated
    # test. Turning them into errors would therefore fail the wrong test;
    # they are counted instead and reported once, at the end, against the
    # run as a whole. A leak does not produce a wrong answer, but it does
    # degrade a server that is deliberately holding every finished job, so
    # it should not be something the suite prints and nobody reads.
    leaks = []
    warnings.simplefilter("always", ResourceWarning)
    previous_hook = warnings.showwarning

    def record(message, category, filename, lineno, *args, **kwargs):
        if issubclass(category, ResourceWarning):
            leaks.append(f"{filename}:{lineno}: {message}")
        else:
            previous_hook(message, category, filename, lineno,
                          *args, **kwargs)

    warnings.showwarning = record

    loader = unittest.TestLoader()
    if names:
        suite = unittest.TestSuite()
        for name in names:
            module = name if name.startswith("test_") else f"test_{name}"
            suite.addTests(loader.loadTestsFromName(module))
    else:
        suite = loader.discover(TESTS, pattern="test_*.py", top_level_dir=TESTS)

    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Skips are reported explicitly. A suite that quietly skipped half of
    # itself and printed OK is the same failure mode this pipeline guards
    # against everywhere else: a clean result that means less than it looks.
    if result.skipped:
        print(f"\n{len(result.skipped)} test(s) skipped:")
        for case, reason in result.skipped:
            print(f"  - {case.id().split('.')[-1]}: {reason}")

    warnings.showwarning = previous_hook
    if leaks:
        print(f"\n{len(leaks)} leaked file handle(s) -- a resource the code "
              f"opened and did not close:")
        for leak in sorted(set(leaks))[:20]:
            print(f"  - {leak}")

    return 0 if (result.wasSuccessful() and not leaks) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
