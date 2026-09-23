# Created by Brainstorm, 2026.
"""
Simulated DNA specimens with a KNOWN answer, built from real hg38.

Writes three haplotypes around five hotspots and a panel BED, plus
truth.tsv, the list of calls a correct run must make:

  somatic   BRAF V600E, KRAS G12D, TP53 R273H, EGFR L858R, PIK3CA H1047R,
            and a 3-bp deletion in TP53 (intronic -- it tests indel calling,
            not interpretation)
  germline  one heterozygous SNP near KRAS, present in tumour AND normal

Every coordinate is checked against the reference base before anything is
written, so a wrong coordinate fails here rather than as a missed call.
make_truth.sh turns the haplotypes into reads at exact allele fractions.

Usage: build_dna_truth.py OUT_DIR   (needs samtools on PATH)
"""
import os, subprocess, sys

REF = os.path.expanduser(
    "~/data/references/hg38/Homo_sapiens_assembly38.fasta")
OUT = sys.argv[1]
PAD = 1000          # sequence simulated around each hotspot
TARGET = 150        # panel BED half-width around each hotspot

SOMATIC = [  # chrom, pos (1-based), ref, alt, label
    ("chr7", 140753336, "A", "T", "BRAF p.V600E"),
    ("chr12", 25245350, "C", "T", "KRAS p.G12D"),
    ("chr17", 7673802, "C", "T", "TP53 p.R273H"),
    ("chr7", 55191822, "T", "G", "EGFR p.L858R"),
    ("chr3", 179234297, "A", "G", "PIK3CA p.H1047R"),
]
# A somatic deletion 120 bp downstream of the TP53 hotspot (inside the
# TP53 target), to prove indels survive every stage too.
DEL_POS = 7673802 + 120
GERMLINE_OFFSET = -90   # het SNP 90 bp upstream of the KRAS hotspot


def fetch(region):
    out = subprocess.run(["samtools", "faidx", REF, region],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.splitlines()[1:]).upper()


def main():
    os.makedirs(OUT, exist_ok=True)
    regions = []
    for chrom, pos, ref, alt, label in SOMATIC:
        start = pos - PAD
        seq = list(fetch(f"{chrom}:{start}-{pos + PAD}"))
        assert seq[PAD] == ref, (label, seq[PAD])
        regions.append(dict(chrom=chrom, start=start, pos=pos, ref=ref,
                            alt=alt, label=label, seq=seq))

    truth = []
    ref_hap, germ_hap, tum_hap = {}, {}, {}
    for r in regions:
        name = f"{r['chrom']}_{r['start']}"
        base = r["seq"][:]
        germ = base[:]
        tum = base[:]
        # somatic SNV
        tum[PAD] = r["alt"]
        truth.append((r["chrom"], r["pos"], r["ref"], r["alt"], "somatic",
                      r["label"]))
        if r["label"].startswith("KRAS"):
            i = PAD + GERMLINE_OFFSET
            g_ref = base[i]
            g_alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[g_ref]
            germ[i] = g_alt
            tum[i] = g_alt
            truth.append((r["chrom"], r["start"] + i, g_ref, g_alt,
                          "germline", "het SNP near KRAS"))
        if r["label"].startswith("TP53"):
            i = PAD + 120
            anchor = base[i - 1]
            deleted = "".join(base[i:i + 3])
            del tum[i:i + 3]
            truth.append((r["chrom"], r["start"] + i - 1, anchor + deleted,
                          anchor, "somatic", "3-bp deletion in TP53"))
        ref_hap[name] = "".join(base)
        germ_hap[name] = "".join(germ)
        tum_hap[name] = "".join(tum)

    def write_fa(path, seqs):
        with open(path, "w") as fh:
            for n, s in seqs.items():
                fh.write(f">{n}\n")
                for k in range(0, len(s), 60):
                    fh.write(s[k:k + 60] + "\n")

    write_fa(f"{OUT}/ref_hap.fa", ref_hap)
    write_fa(f"{OUT}/germ_hap.fa", germ_hap)
    write_fa(f"{OUT}/tum_hap.fa", tum_hap)

    with open(f"{OUT}/panel.bed", "w") as fh:
        for r in sorted(regions, key=lambda r: (r["chrom"], r["pos"])):
            fh.write(f"{r['chrom']}\t{r['pos'] - TARGET - 1}\t"
                     f"{r['pos'] + TARGET + 150}\t{r['label'].split()[0]}\n")
    with open(f"{OUT}/truth.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tclass\tlabel\n")
        for t in truth:
            fh.write("\t".join(map(str, t)) + "\n")
    for t in truth:
        print(*t, sep="\t")


main()
