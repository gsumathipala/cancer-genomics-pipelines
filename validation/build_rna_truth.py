# Created by Brainstorm, 2026.
"""
Simulated RNA library with two KNOWN fusions, from real hg38 + GENCODE v44.

  EML4::ALK  variant 1 : EML4 exons 1-13 joined to ALK exons 20-end
  BCR::ABL1  e14a2     : BCR exons 1-14 joined to ABL1 exons 2-end

Background: ~300 canonical protein-coding transcripts, plus the WILD-TYPE
transcripts of all four partners, so the library looks like a library --
Arriba refuses a library with no normal reads, which is why a fusion-only
simulation cannot test it.

GENCODE v44 is used because the STAR index this was first run against was
built from it: transcripts from the same release as the index. Override
with GTF=... in the environment when yours differs.

Usage: build_rna_truth.py OUT_DIR   (needs samtools on PATH)
"""
import collections, os, random, subprocess, sys

REF = os.path.expanduser("~/data/references/hg38/Homo_sapiens_assembly38.fasta")
GTF = os.environ.get("GTF") or os.path.expanduser(
    "~/data/references/gencode/gencode.v44.primary_assembly.annotation.gtf")
OUT = sys.argv[1]
PARTNERS = {"EML4", "ALK", "BCR", "ABL1"}
random.seed(7)
PRIMARY = {f"chr{i}" for i in list(range(1, 23)) + ["X"]}

tx_exons = collections.defaultdict(list)   # tid -> [(chrom,start,end,strand,exon_no)]
tx_gene = {}
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if f[2] != "exon" or 'Ensembl_canonical' not in f[8] or \
                'gene_type "protein_coding"' not in f[8]:
            continue
        a = f[8]
        gene = a.split('gene_name "')[1].split('"')[0]
        tid = a.split('transcript_id "')[1].split('"')[0]
        exon_no = int(a.split('exon_number ')[1].split(';')[0])
        tx_exons[tid].append((f[0], int(f[3]), int(f[4]), f[6], exon_no))
        tx_gene[tid] = gene

by_gene = {g: t for t, g in tx_gene.items()}
assert PARTNERS <= set(by_gene), PARTNERS - set(by_gene)
others = sorted(g for g in by_gene if g not in PARTNERS
                and tx_exons[by_gene[g]][0][0] in PRIMARY)
background = random.sample(others, 300)

# Fetch every exon needed in ONE samtools call.
wanted = [by_gene[g] for g in list(PARTNERS) + background]
regions = []
for tid in wanted:
    for chrom, s, e, strand, n in tx_exons[tid]:
        regions.append(f"{chrom}:{s}-{e}")
with open(f"{OUT}/regions.txt", "w") as fh:
    fh.write("\n".join(regions) + "\n")
fa = subprocess.run(["samtools", "faidx", REF, "-r", f"{OUT}/regions.txt"],
                    capture_output=True, text=True, check=True).stdout
seqs, name = {}, None
for line in fa.splitlines():
    if line.startswith(">"):
        name = line[1:]; seqs[name] = []
    else:
        seqs[name].append(line)
seqs = {k: "".join(v).upper() for k, v in seqs.items()}
RC = str.maketrans("ACGTN", "TGCAN")

def exon_seq(chrom, s, e, strand):
    x = seqs[f"{chrom}:{s}-{e}"]
    return x.translate(RC)[::-1] if strand == "-" else x

def transcript(tid, first=None, last=None):
    """Spliced mRNA (5'->3'), optionally only exon numbers first..last."""
    ex = sorted(tx_exons[tid], key=lambda r: r[4])
    return "".join(exon_seq(c, s, e, st) for c, s, e, st, n in ex
                   if (first is None or n >= first) and (last is None or n <= last))

def breakpoint(tid, exon_no, side):
    """Genomic coordinate of the fusion junction on the exon's boundary."""
    for c, s, e, st, n in tx_exons[tid]:
        if n == exon_no:
            # 5' partner: junction at the exon's 3' end; 3' partner: 5' end
            if side == "5p":
                return c, (e if st == "+" else s)
            return c, (s if st == "+" else e)

fusions = [("EML4", 13, "ALK", 20), ("BCR", 14, "ABL1", 2)]
with open(f"{OUT}/transcripts_bg.fa", "w") as bg, \
        open(f"{OUT}/transcripts_fusion.fa", "w") as fu, \
        open(f"{OUT}/truth.tsv", "w") as truth:
    for g in background + sorted(PARTNERS):
        bg.write(f">{g}\n{transcript(by_gene[g])}\n")
    truth.write("fusion\tbreakpoint1\tbreakpoint2\n")
    for g5, e5, g3, e3 in fusions:
        t5, t3 = by_gene[g5], by_gene[g3]
        seq = transcript(t5, last=e5) + transcript(t3, first=e3)
        fu.write(f">{g5}--{g3}\n{seq}\n")
        c1, p1 = breakpoint(t5, e5, "5p")
        c2, p2 = breakpoint(t3, e3, "3p")
        truth.write(f"{g5}::{g3}\t{c1}:{p1}\t{c2}:{p2}\n")
        print(f"{g5}(ex{e5})::{g3}(ex{e3})  {c1}:{p1} -> {c2}:{p2}  "
              f"transcript {len(seq)} nt  [{t5} / {t3}]")
print(f"background: {len(background)} genes + 4 wild-type partners")
