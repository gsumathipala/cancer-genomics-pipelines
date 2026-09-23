#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
Grade a run against the known-answer truth set. Exit 0 only on a clean pass.

    check_truth.py dna TRUTH_TSV FILTERED_VCF [--tumour-only]
    check_truth.py rna TRUTH_TSV ARRIBA_FUSIONS_TSV

DNA. Every somatic truth variant must be called at its exact position and
allele, FILTER=PASS, with a tumour allele fraction near 0.30 (the fraction
make_truth.sh builds in). The germline SNP must be ABSENT from a matched-
normal run -- the normal exists to subtract it -- and is expected, at ~0.50,
in a tumour-only run, where nothing can tell a novel germline variant from
a somatic one. Any other PASS call is reported as a false positive.

RNA. Every truth fusion must be called with both breakpoints exact; any
other call is reported.

The tumour's column is found from Mutect2's own ##tumor_sample header, not
assumed from position -- getting that backwards is precisely the kind of
failure this is meant to catch.

Standard library only.
"""
import csv
import gzip
import sys

AF_TOLERANCE = 0.08       # around the constructed 0.30 / 0.50


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_truth(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def vcf_calls(path):
    """{(chrom, pos, ref, alt): (filter, tumour_af)}"""
    calls, tumour, samples = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##tumor_sample="):
                tumour = line.strip().split("=", 1)[1]
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                keys = f[8].split(":")
                col = samples.index(tumour) if tumour in samples else 0
                values = dict(zip(keys, f[9 + col].split(":")))
                af = float(values.get("AF", "nan").split(",")[0])
                for alt in f[4].split(","):
                    calls[(f[0], f[1], f[3], alt)] = (f[6], af)
    return calls


def check_dna(truth_path, vcf_path, tumour_only):
    truth = read_truth(truth_path)
    calls = vcf_calls(vcf_path)
    failures, expected = 0, set()
    for t in truth:
        key = (t["chrom"], t["pos"], t["ref"], t["alt"])
        expected.add(key)
        got = calls.get(key)
        want_present = t["class"] == "somatic" or tumour_only
        target = 0.30 if t["class"] == "somatic" else 0.50
        if want_present:
            ok = bool(got) and got[0] == "PASS" and \
                abs(got[1] - target) <= AF_TOLERANCE
            state = (f"PASS af={got[1]:.3f}" if got else "NOT CALLED")
        else:
            ok = not got or got[0] != "PASS"
            state = "correctly absent" if ok else f"CALLED af={got[1]:.3f}"
        failures += not ok
        print(f"{'ok  ' if ok else 'FAIL'}  {t['label']:26} {state}")
    extras = [k for k, v in calls.items()
              if v[0] == "PASS" and k not in expected]
    for k in extras:
        print(f"FAIL  false positive          {':'.join(k)}")
    return failures + len(extras)


def check_rna(truth_path, arriba_path):
    truth = read_truth(truth_path)
    with open(arriba_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    called = {(r["breakpoint1"], r["breakpoint2"]):
              f"{r['#gene1']}::{r['gene2']}" for r in rows}
    failures = 0
    for t in truth:
        name = called.get((t["breakpoint1"], t["breakpoint2"]))
        ok = name == t["fusion"]
        failures += not ok
        print(f"{'ok  ' if ok else 'FAIL'}  {t['fusion']:12} "
              f"{t['breakpoint1']} -> {t['breakpoint2']}"
              f"{'' if ok else '  NOT CALLED at these breakpoints'}")
    wanted = {(t["breakpoint1"], t["breakpoint2"]) for t in truth}
    for bp, name in called.items():
        if bp not in wanted:
            failures += 1
            print(f"FAIL  extra call   {name} {bp[0]} -> {bp[1]}")
    return failures


def main(argv):
    if len(argv) < 3 or argv[0] not in ("dna", "rna"):
        print(__doc__)
        return 2
    if argv[0] == "dna":
        failures = check_dna(argv[1], argv[2], "--tumour-only" in argv)
    else:
        failures = check_rna(argv[1], argv[2])
    print(f"\n{'PASS' if not failures else f'FAIL ({failures})'}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
