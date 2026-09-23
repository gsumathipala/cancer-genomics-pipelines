<!-- Created by Brainstorm, 2026. -->

# Testing

**Research use only.** Nothing here validates this software for clinical
use. A passing suite says the code does what its authors intended; it says
nothing about whether that intention is right for your assay.

```bash
python3 run_tests.py            # everything, about 20 seconds
python3 run_tests.py -v         # name every test as it runs
python3 run_tests.py regressions   # one module
```

Nothing needs installing first. The suite uses the standard library's
`unittest`, for the same reason the analysis scripts import nothing else:
somebody who has just cloned this should be able to check that it works
before they spend an hour building a conda environment. Tests that need
Flask skip themselves, announce the skip, and do not fail.

---

## What is tested, and what deliberately is not

No aligner, caller or annotator is ever invoked. BWA, STAR, GATK, Arriba
and PCGR are third-party software with their own test suites; running them
here would turn a twenty-second check into an overnight one and would test
somebody else's code.

What is tested is **everything this bundle decides**:

| Tier | Module | What it pins |
|---|---|---|
| Configuration | `test_panel_profiles.py` | Which settings a panel applies, what an explicit flag overrides, how a BED footprint becomes a TMB denominator |
| Reporting | `test_fusion_report.py`, `test_rna_qc.py` | How a fusion is ranked, what counts as a supporting read, when a library is called inadequate |
| Interface | `test_webapp_routes.py`, `test_worksheet.py` | Every route answers, a worksheet row becomes a valid run, a bad row is rejected before anything is queued |
| Commands | `test_argv_builders.py`, `test_coverage.py` | Every flag the interface emits is a flag the target script accepts |
| Regressions | `test_regressions.py` | One named test per bug that actually reached this code |
| Bundle | `test_bundle_integrity.py` | Everything compiles, everything is watermarked, the stdlib-only promise holds, the documentation names files that exist |

The command tier deserves a word. It does not run commands; it compares
the flags the web interface emits against the flags each script's `--help`
advertises. A mismatch there is the failure that otherwise appears twenty
minutes into a run as `unrecognized arguments`, with the alignment already
paid for.

---

## Why the regression tier is the important one

Every test in `test_regressions.py` names the bug it prevents, and where
that bug was a silent one, says what the user would have seen instead of an
error. That matters more here than it does in most software, because of
the kind of wrong answer this pipeline can produce.

A crash is cheap. Somebody sees it, and nobody reports a result. The
expensive failure is the one that produces a clean, well-formed, entirely
plausible report that is wrong — and every bug in that module was one of
those:

- STAR counted 116 chimeric reads and wrote none of them into the BAM. The
  fusion caller reads only the BAM, found nothing, and exited zero. The
  output was an empty fusion table, indistinguishable from a true negative
  for a specimen that may well have carried a fusion.
- The report step looked for the COSMIC-annotated VCF by name and gave up
  when it was absent, reporting "no VCF on disk" while the filtered calls
  sat beside it.
- An RNA run's report step searched for a DNA manifest, did not find one,
  and skipped itself with a message that was true of the pattern it
  searched for and misleading about the run.
- Report links addressed artefacts by position in a list that grows while
  a run is in flight, so a link captured early could serve a different
  file later.
- A minimum allele fraction of `0` — report everything — was read as "not
  set" and replaced by the default, which then filtered out the variants
  the user had explicitly asked to keep.

None of these would have failed a run. All of them would have produced a
report. That is why they are each pinned by name.

One bug in that module is of a different kind and is marked as such: a
finished job held its subprocess pipe open. The manager deliberately keeps
every `Job` for the life of the server so finished runs stay listable, so
each completed run also kept a file descriptor — a server working through
a long worksheet accumulated them until it ran out. It reports no wrong
answer; it degrades a long-lived server, and it is invisible until the
limit is hit, at which point it presents as something else entirely.

## Leaked file handles fail the run

`run_tests.py` counts `ResourceWarning`s and fails the run on any that
appear, listing them at the end. They are not turned into errors, because
a leaked handle is reported whenever the garbage collector happens to
notice it — usually inside an unrelated test, which would then fail for
something it did not do.

This is why the suite reads files through `read_text()` and `read_json()`
in `tests/helpers.py` rather than `open(path).read()`: the tempting
one-liner leaks a handle, and a suite that enforces a rule has to keep
it. The check found nine such handles in the tests themselves the first
time it ran.

---

## Adding a test

Put it in the tier it belongs to, and say in a comment what it is for. A
test whose name and body do not explain which failure it prevents becomes,
within a year, something nobody dares delete and nobody understands.

`tests/helpers.py` has the fixtures: `TempCase` for a scratch directory
that cleans itself up, plus builders for BED files, FASTA indexes, Arriba
output and STAR logs. `tests/webcase.py` boots the Flask app against a
scratch run directory, so route tests never touch real runs.

**If a test fails on arrival, that is a finding, not a nuisance.** Two of
the cases here began as failures and were traced to the code rather than
the test:

- Two different definitions of "supporting reads" in `fusion_report.py`,
  one including discordant mates and one excluding them. Both were
  correct — PCGR's `SplitReads` column is defined against split reads
  alone — but neither was labelled, so a reader comparing the HTML table
  against the PCGR input saw two numbers for the same fusion. Both are now
  labelled, and both definitions are pinned by a test.
- STAR index discovery accepted a directory on one sentinel file while the
  alignment step demanded four, so a half-built index could be selected
  automatically and then rejected a moment later. Both now ask the same
  question.

Do not adjust a test to make it pass until you know which side is wrong.
