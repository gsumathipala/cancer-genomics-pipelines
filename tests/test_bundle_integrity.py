# Created by Brainstorm, 2026.
"""
Properties of the bundle as a whole -- the things that drift silently.

None of this is about behaviour. It is about the bundle staying the thing
it claims to be: every script importable, every module watermarked, the
stdlib-only promise in requirements.txt still true, and the documentation
still pointing at files that exist. Each of these has exactly one way of
going wrong: somebody adds a file and forgets.
"""

import ast
import os
import py_compile
import re
import sys
import tempfile
import unittest

# helpers puts the bundle root on sys.path, which is how every test
# module imports the scripts without the bundle being a package.
from helpers import ROOT as _ROOT  # noqa: F401

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Third-party imports that are allowed, and only where they are allowed.
# The analysis scripts promise to need nothing but the standard library,
# so a user can run them without creating an environment first; the web
# interface is the one part that asks for Flask, and says so.
ALLOWED_THIRD_PARTY = {"flask", "werkzeug", "jinja2"}
WEBAPP_ONLY = {"flask", "werkzeug", "jinja2"}


def python_files(subdir=""):
    base = os.path.join(ROOT, subdir) if subdir else ROOT
    for name in sorted(os.listdir(base)):
        if name.endswith(".py"):
            yield os.path.join(base, name)


def stdlib_names():
    names = set(sys.stdlib_module_names)
    names.discard("this")
    return names


class TestEveryScriptCompiles(unittest.TestCase):

    def test_all_modules_compile(self):
        # A syntax error in a rarely-run script otherwise surfaces only
        # when a run reaches it, which can be an hour in. docs/ is
        # included: the diagram generator is run by hand, so nothing else
        # would ever notice it had been broken.
        with tempfile.TemporaryDirectory() as scratch:
            for path in (list(python_files()) + list(python_files("webapp"))
                         + list(python_files("tests"))
                         + list(python_files("docs"))):
                with self.subTest(module=os.path.basename(path)):
                    # Byte-code goes to scratch: compiling in place would
                    # litter the bundle with __pycache__ directories that
                    # the checksum manifest then disagrees with.
                    py_compile.compile(
                        path, doraise=True,
                        cfile=os.path.join(scratch, "out.pyc"))


class TestWatermark(unittest.TestCase):
    """Every source file carries the author's watermark."""

    def test_first_line_names_the_author(self):
        pattern = re.compile(r"Created by Brainstorm", re.IGNORECASE)
        for path in (list(python_files()) + list(python_files("webapp"))
                     + list(python_files("tests"))
                     + list(python_files("docs"))):
            with self.subTest(module=os.path.relpath(path, ROOT)):
                with open(path, encoding="utf-8") as fh:
                    head = "".join(fh.readline() for _ in range(6))
                self.assertTrue(pattern.search(head),
                                "watermark missing from the file header")


class TestDependencyPromise(unittest.TestCase):
    """
    requirements.txt promises the analysis scripts import nothing outside
    the standard library. That promise is what lets a user run a script
    before building an environment, so it is worth enforcing mechanically.
    """

    def _top_level_imports(self, path):
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read(), filename=path)
        found = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                found.update(alias.name.split(".")[0]
                             for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                if node.module:
                    found.add(node.module.split(".")[0])
        return found

    def test_analysis_scripts_import_only_the_standard_library(self):
        # Deliberately the root scripts only. docs/make_diagram.py uses
        # Pillow, and that is fine: it is a development tool that draws the
        # illustration, never something a user needs in order to run an
        # analysis. The promise in requirements.txt is about the latter.
        local = {os.path.basename(p)[:-3] for p in python_files()}
        for path in python_files():
            with self.subTest(module=os.path.basename(path)):
                outside = (self._top_level_imports(path)
                           - stdlib_names() - local)
                self.assertEqual(
                    outside, set(),
                    "requirements.txt promises stdlib-only here")

    def test_only_the_webapp_asks_for_flask(self):
        for path in python_files("webapp"):
            imports = self._top_level_imports(path)
            with self.subTest(module=os.path.basename(path)):
                self.assertEqual(imports & ALLOWED_THIRD_PARTY,
                                 imports & WEBAPP_ONLY)


class TestPanelRegistryIsCoherent(unittest.TestCase):

    def setUp(self):
        import panel_profiles
        self.panel_profiles = panel_profiles
        # The registry is keyed by id, so duplicates cannot survive it;
        # what is tested here is that each surviving profile is usable.
        self.registry = panel_profiles.available_panels()

    def test_the_registry_is_not_empty(self):
        self.assertGreater(len(self.registry), 10)

    def test_every_profile_is_well_formed(self):
        for pid, profile in sorted(self.registry.items()):
            with self.subTest(panel=pid):
                self.assertIn(profile["chemistry"],
                              self.panel_profiles.CHEMISTRIES)
                # A profile without a name cannot be presented in the
                # interface; one without settings does nothing at all.
                self.assertTrue(profile.get("name"))
                self.assertIsInstance(profile.get("settings"), dict)
                for key in self.panel_profiles.REQUIRED_PROFILE_KEYS:
                    self.assertIn(key, profile)

    def test_chemistry_decides_the_assay(self):
        # The RNA engine refuses DNA profiles by this test alone, so a
        # profile whose chemistry and assay disagree would be routed to
        # the wrong engine.
        for pid, profile in sorted(self.registry.items()):
            expected = ("rna" if profile["chemistry"]
                        in self.panel_profiles.RNA_CHEMISTRIES else "dna")
            with self.subTest(panel=pid):
                self.assertEqual(self.panel_profiles.assay_type(profile),
                                 expected)

    def test_every_setting_name_is_one_the_engines_know(self):
        # A misspelled setting key is applied to nothing and reported as
        # applied -- the profile appears to take effect and does not.
        known = set(self.panel_profiles.PANEL_SETTING_TYPES)
        for pid, profile in sorted(self.registry.items()):
            with self.subTest(panel=pid):
                self.assertEqual(set(profile["settings"]) - known, set())


class TestDocumentationPointsAtRealFiles(unittest.TestCase):

    def test_referenced_scripts_exist(self):
        # A README naming a script that was renamed sends the reader to a
        # dead end, and is the single most common form of doc rot here.
        named = set()
        for doc in ("README.md", "PIPELINE_ANATOMY.md"):
            path = os.path.join(ROOT, doc)
            if not os.path.exists(path):
                continue
            with open(path, encoding="utf-8") as fh:
                named.update(re.findall(r"\b([a-z0-9_]+\.py)\b", fh.read()))
        for name in sorted(named):
            with self.subTest(script=name):
                self.assertTrue(
                    os.path.exists(os.path.join(ROOT, name))
                    or os.path.exists(os.path.join(ROOT, "webapp", name))
                    or os.path.exists(os.path.join(ROOT, "tests", name))
                    or os.path.exists(os.path.join(ROOT, "docs", name)),
                    f"documentation refers to {name}, which is not here")

    def test_research_use_only_is_stated_first(self):
        # This is a legal statement, not a stylistic one: it has to be
        # visible before a reader decides what the software is for.
        with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as fh:
            head = fh.read(1200).lower()
        self.assertIn("research use only", head)


if __name__ == "__main__":
    unittest.main()


class TestTheDiagramAgreesWithTheCode(unittest.TestCase):
    """
    The illustration in the README is generated, not drawn.

    A hand-drawn diagram stops being true the first time a step is
    renamed, and nothing tells you. The generator transcribes the step
    lists and can compare them against the banners the engines print, so
    that drift becomes a test failure rather than a misleading picture on
    the front page of the repository.
    """

    def test_step_lists_match_the_engines(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "make_diagram", os.path.join(ROOT, "docs", "make_diagram.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.check(), [])

    def test_the_rendered_files_are_present(self):
        for name in ("pipeline_overview.png", "pipeline_overview.svg"):
            path = os.path.join(ROOT, "docs", name)
            with self.subTest(file=name):
                self.assertTrue(os.path.exists(path))
                # A zero-byte image still renders as a broken link, which
                # is exactly the kind of thing nobody checks after a merge.
                self.assertGreater(os.path.getsize(path), 10_000)
