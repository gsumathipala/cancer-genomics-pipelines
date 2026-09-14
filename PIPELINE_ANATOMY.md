<!-- Created by Brainstorm, 2026. -->
# Anatomy of a cancer genomics pipeline

*How this one is put together, and what it should teach you about building
your own.*

This is not the user guide — that is [README.md](README.md). This is the
document I would want if I had been handed a pile of bioinformatics tools
and told to turn them into something a hospital could run. It explains
**what each tool is for, how they chain together, and — mostly — the
things that go wrong quietly.**

Read it if you are a scientist who can code a bit, or a programmer who has
landed in a genomics lab. You do not need to know the tools already.

---

## 1. Why this kind of software is unusual

Most software tells you when it fails. A web request 500s, a compiler
errors, a test goes red. Bioinformatics pipelines mostly do not, and that
single fact should shape every design decision you make.

Three properties conspire:

**Runs are long.** Aligning one exome takes hours. A mistake in step 2 is
discovered at step 9, tomorrow. So errors must be caught *before* the
expensive part, not wherever they happen to surface.

**There is no ground truth.** When your web app returns the wrong user, the
user notices. When your pipeline returns the wrong variant list, nobody
notices — there is nothing to compare it against. The output of a broken
run and a correct run look *identical*: a well-formed VCF with plausible
numbers in it.

**Every tool is permissive.** Hand a variant caller the wrong reference and
it will not refuse; it will call variants. Hand it an empty target region
and it will call nothing and exit zero. The tools are libraries, not
guardians. The guarding is your job.

> **The first principle.** In this domain, the dangerous bug is not the one
> that crashes. It is the one that produces a confident, well-formatted,
> wrong answer — and there is a *lot* of surface area for those.

---

## 2. The organising idea: make silence impossible

Almost every design decision in this codebase comes from one question:

> *If this goes wrong, how would anyone find out?*

If the honest answer is "they wouldn't", that is the thing to fix — not by
making it impossible (often you can't), but by making it **loud**. Three
tactics recur throughout, and they are the transferable part:

1. **Check early, where it is cheap.** A PCGR reference bundle without a
   VEP cache is refused at argument-parsing time, not five hours later when
   the report step starts.
2. **Say what you did *not* do.** Every report lists the things it could not
   measure, by name. A missing section reads as "nothing found"; a named
   absence reads as "never measured". Those are wildly different clinically.
3. **Record provenance.** Every run writes a manifest: every setting, every
   reference file, every tool version. Not for you — for whoever reads the
   result in a year and needs to know whether to trust it.

---

## 3. The map

Two branches share a front end, a genome, and a configuration system.

```
                         ┌──────────────────────┐
     FASTQ files  ──────►│  fastq_qc_clean.py   │  STAGE 1 (shared)
                         │  trim, filter, QC    │
                         └──────────┬───────────┘
                                    │ run manifest (JSON)
                    ┌───────────────┴────────────────┐
                    ▼                                ▼
        ┌───────────────────────┐        ┌────────────────────────┐
        │   align_reads.py      │        │   fusion_calling.py    │
        │   bwa-mem2            │        │   ├─ align_rna.py      │
        │   STAGE 2 (DNA only)  │        │   │  STAR (splice-     │
        └──────────┬────────────┘        │   │  aware, chimeric)  │
                   │ BAM                 │   ├─ Arriba (fusions)  │
                   ▼                     │   ├─ rna_qc_report.py  │
   ┌───────────────────────────────┐     │   └─ fusion_report.py  │
   │ comprehensive_variant_calling │     └───────────┬────────────┘
   │  dedup → BQSR → Mutect2 →     │                 │
   │  filter → COSMIC → SnpEff →   │                 │
   │  coverage → PCGR              │                 │
   └───────────────┬───────────────┘                 │
                   │ VCF                             │ fusions.tsv
                   └──────────────┬──────────────────┘
                                  ▼
                          ┌───────────────┐
                          │     PCGR      │  clinical interpretation
                          │  tiers, TMB,  │  (takes either, or BOTH
                          │  fusions      │   for one specimen)
                          └───────────────┘
```

`pipeline_orchestrator.py` runs the chain; `--assay dna|rna` picks the
branch. `webapp/` is a Flask front end over the same scripts.

**Why two engines and not one with a flag?** Because after alignment,
*nothing* is shared. Duplicate marking, base recalibration, and somatic
calling are either meaningless or actively wrong on RNA. A single engine
with `if assay == "rna"` scattered through fourteen steps would be a file
nobody could reason about. Separate engines, shared front end.

---

## 4. The cast of tools

You do not need to memorise these. You need to know *what question each
one answers*.

| Tool | Question it answers | Branch |
|---|---|---|
| **fastp** | Are the reads clean? Trim adapters, low-quality ends, poly-G | both |
| **FastQC** | What do the reads look like? (report only) | both |
| **bwa-mem2** | Where in the genome does each read come from? | DNA |
| **STAR** | Same — but allowing a read to span an exon–exon junction | RNA |
| **samtools** | The Swiss army knife: sort, index, depth, stats | both |
| **GATK MarkDuplicates** | Which reads are PCR copies of one molecule? | DNA |
| **GATK BQSR** | Are the instrument's quality scores actually right? | DNA |
| **GATK Mutect2** | Which positions differ from the reference, somatically? | DNA |
| **GATK FilterMutectCalls** | Which of those calls survive scrutiny? | DNA |
| **bcftools** | VCF manipulation, annotation, statistics | DNA |
| **SnpEff** | What gene and what consequence? | DNA |
| **Ensembl VEP** | Same, more thoroughly (run by PCGR) | DNA |
| **MSIsensor2** | Is this tumour microsatellite unstable? | DNA |
| **Arriba** | Which reads suggest two genes have fused? | RNA |
| **PCGR** | Which findings are clinically actionable? | both |

Plus four reports written here rather than borrowed, because nothing
existing answered the question: `coverage_report.py`,
`rna_qc_report.py`, `fusion_report.py`, and the PDF writer.

---

## 5. Walkthrough: the DNA branch

Fourteen steps. Here is what each is *for*, and what it costs you when it
is wrong.

**1. QC and trimming (fastp).** Adapter sequence is the sequencer reading
past the end of a short fragment. Left in, it aligns somewhere and looks
like sequence. Poly-G is a NextSeq-specific artefact: its two-colour
chemistry reads "no signal" as G, so a dying cluster produces a run of Gs.

**2. Alignment (bwa-mem2).** Each read is placed in the genome. Everything
downstream is a claim about what happened at a position, so a
mis-alignment is a wrong answer with full confidence behind it.

**3. Duplicate marking.** PCR amplifies each original molecule many times.
Twenty reads from one molecule are *one* observation, not twenty, and
counting them as twenty makes an artefact look like a confident variant.
Marked, not removed — see the "tag, don't delete" pearl below.

**4. BQSR.** The instrument's quality scores are systematically biased.
BQSR measures the bias against known variant sites and corrects it.

**5. Somatic calling (Mutect2).** Proposes candidate variants. It is
deliberately permissive; filtering is a separate step.

**6. Read-orientation model.** FFPE fixation deaminates cytosine, producing
C>T/G>A artefacts on one strand. This models that so the filter can
separate chemistry from biology.

**7. Contamination estimate.** How much of this "tumour" is another
sample? Without measuring it, filtering assumes zero.

**8. MSI.** Microsatellite instability — a therapy-selection marker.

**9. Filtering.** Turns candidates into calls. Then a depth floor, then the
foldback tagger (below).

**10–11. Annotation.** COSMIC identifiers, then gene and consequence.

**12. Coverage.** *Which parts of the panel were actually sequenced deeply
enough to have found a variant.* The most under-appreciated step here.

**13. PCGR.** Clinical interpretation: tiers, TMB, signatures.

**14. Summary.**

---

## 6. Walkthrough: the RNA branch

Six steps, and the differences from DNA are the lesson.

**1. QC.** Same tool, lower length floor — FFPE RNA is fragmented, and a
50 bp floor inherited from DNA throws away much of a usable library.

**2. Alignment (STAR).** A read crossing an exon–exon junction has **no
contiguous match** anywhere in the genome. bwa-mem2 cannot represent that;
it soft-clips the read back to one exon or discards it. Silently. And
those are exactly the reads that could evidence a fusion.

**3. Sort.** For QC and IGV only — the fusion caller wants the *unsorted*
BAM, because it reads mate pairs together and sorting separates them.

**4. Fusion calling (Arriba).** Reads whose two halves land on two
different genes.

**5. Library adequacy.** Does this library support a negative result at
all?

**6. Fusion report.** Ranking by clinical importance, not caller
confidence.

**No duplicate marking, and no flag to enable it.** On RNA, twenty
identical reads from a highly expressed gene are twenty genuine
observations of an abundant transcript. Removing them destroys the signal.
There is no correct value for that setting, so the setting does not exist.

---

## 7. The pearls

Each of these cost somebody a day. They are grouped by the kind of mistake
they represent, because the *categories* generalise even when the specific
tools do not.

### I. Reference data is where runs die

**Contig naming will get you at least once.** Two conventions exist for
naming chromosomes: UCSC (`chr1`, `chrM`) and Ensembl (`1`, `MT`). They
are not interchangeable, and mixing them does not always error.

- A COSMIC VCF in Ensembl naming against a UCSC genome annotates
  **nothing at all** — and the run reports success.
- A panel BED in the wrong naming makes GATK abort (loud, fine) but makes
  the *coverage report* find no reads in any region — which reads as "none
  of the panel was covered". A confident, completely wrong clinical
  statement.

> **Pearl:** detect the convention, translate, and say so. Never edit the
> user's file — write a corrected copy and name it in the log.

**Verify downloads by size, not by existence.** A truncated 20 GB
reference passes `os.path.exists()` forever. Every download here checks
the server's `Content-Length` before renaming the file into place.

**Version-lock what is baked in.** The STAR index is built from a genome
*and* an annotation *and* a read length. All three are baked in. An index
built for 100 bp reads works fine on 150 bp reads and quietly loses
junction sensitivity. So the build records what it was made from, and
later runs compare.

**Dependency versions can fail in the middle.** GATK's Spark tools do not
run on JDK 24+, and the failure lands at duplicate marking — *after*
alignment has burned two hours. `gatk --version` passes either way. So the
installer parses `java -version` and enforces the supported range up
front.

### II. The chemistry decides the algorithm

This is the one that most surprises programmers: **the correct code
depends on how the DNA was prepared in the lab.**

**Amplicon panels must not have duplicates marked.** Amplicon reads begin
and end at primer coordinates *by construction*, so every read from one
amplicon looks like a duplicate of every other. MarkDuplicates flags
nearly the whole library; the caller then works from a handful of
surviving reads. A catastrophic, silent loss of depth whose only symptom
is an implausibly thin VCF.

**BQSR needs data.** It fits an error model from covered sites. On a 20 kb
hotspot panel there are nowhere near enough, and the model it produces is
noise applied to real quality scores.

**Interval padding depends on chemistry.** Capture probes pull in fragments
that extend past the target, so a variant at the first base of an exon
needs ~100 bp of flanking sequence to be callable. On an amplicon panel
the target *is* the amplicon, and padding walks straight into primer
sequence — where the bases come from the oligo, not the patient.

> **Pearl:** if a setting's correct value depends on the wet lab, do not
> give it a default and hope. Name the assay, and derive the settings from
> that. This is what `panel_profiles.py` exists for.

### III. Every denominator is a lie until you check it

**The TMB story is the perfect worked example.** Tumour mutational burden
is mutations per megabase — a count divided by the size of the region you
sequenced. PCGR assumes **34 Mb** for a "targeted" assay unless told
otherwise, because that is exome-sized.

Run a 4.5 Mb panel and your TMB is ~7.5× too low. Run a 2 Mb panel and it
is ~17× too low. **PCGR does not warn you.** It prints a number.

The fix here: measure the footprint from the BED actually being used —
merging overlapping probe tiles first, since summing them unmerged
inflates the denominator and deflates TMB in the same silent direction.

**And know when a number is not interpretable at all.** Below ~1 Mb of
target, a single extra artefact moves TMB by several mutations per Mb. The
pipeline switches TMB off below that threshold and *says why*, rather than
printing a meaningless figure.

> **Pearl:** find every place your system divides by something. Ask where
> that denominator came from. In science code it is usually a default that
> was right for the author's data.

### IV. Absence of evidence is not evidence of absence

**This is the most important idea in the whole codebase.**

A variant caller reports variants. It does not report "I looked here and
found nothing" versus "I never looked here". Both appear in the output as
*nothing at all*.

So a region the sequencing never reached and a region that is genuinely
wild type are **identical in the VCF, identical in the report, and
identical to the clinician reading it.** One means "no mutation"; the
other means "no information". Those are opposite clinical conclusions.

Hence `coverage_report.py`: measure every target region against the depth
a call would need, and name the ones that fall short.

The RNA branch has the same hole and worse. A degraded FFPE RNA extraction
produces a clean, well-formed, **empty** fusion table — indistinguishable
from a specimen that carries no fusion. Hence `rna_qc_report.py`, which
ends not with a table of numbers but with a sentence:

> *3 of 6 checks did not clear their bar. An empty fusion table from this
> library should NOT be reported as a negative without explaining these.*

> **Pearl:** in any system that reports findings, ask what an empty result
> means. If "we found nothing" and "we couldn't look" render identically,
> you have a serious problem and it will not show up in testing.

And a corollary that shapes the whole interface: **library adequacy is
printed before the findings, everywhere** — in the HTML reports, in the
PDF, in the ordering of links on the run page. A reader who sees an empty
table first has already concluded "negative" before reaching the caveat.

### V. Artefacts that look like biology

**The foldback insertions.** During library prep, a fragment end can fold
back on itself and be extended, producing a sequence followed by its own
reverse complement. Aligned, that looks like an *insertion*.

These are not rare and not obviously wrong. On the panel this was built
for:

- 95.6% of insertions longer than 10 bp were reverse-complement copies of
  the adjacent reference (0.2% in the forward direction)
- 1,882 of 2,047 insertions were this artefact
- they accounted for **206 of 270 PASS frameshift/stop calls**, including
  BRCA1, BRCA2, PALB2, PTEN and ATM frameshifts that PCGR tiered as
  oncogenic

No generic filter catches them. They are not low-depth, not strand-biased,
not in tandem repeats, and they reach allele fractions above 0.3.

> **Pearl:** when a filter is doing nothing, ask whether the artefact you
> are worried about actually has the property you are filtering on.

**Off-target reads.** Every capture protocol leaves stray reads scattered
across the genome. Without a target restriction, the caller walks the whole
genome and calls on them. Measured on that same 4.5 Mb panel:

| | raw candidates | PASS calls |
|---|---|---|
| without `--intervals` | 635,414 | 127,081 |
| with `--intervals` | 20,761 | ~1,400 |

Two orders of magnitude of noise, and every one of those 127,081 looks
like a result.

**Tag, don't delete.** Notice that the foldback filter *tags* rather than
removes. Genuine foldback inversions occur in cancer genomes, so the
pattern is *evidence* of artefact, not proof. Every tested insertion gets
an `INFO/FBMATCH` value recording the match fraction, so a reviewer sees
the evidence instead of trusting a threshold they cannot inspect.

> **Pearl:** deleting data because it is *probably* wrong destroys the
> reviewer's ability to disagree with you. Annotate instead, and let the
> threshold be visible.

### VI. The interface is part of the instrument

**Two patients in one directory.** The auto-discovery mode originally took
the first FASTQ pair as the tumour and the second as its matched normal.
Point it at a directory holding two *different patients* and it produces
somatic calls on the genetic difference between two unrelated people:
their differing germline variants reported as somatic, their shared real
mutations subtracted. The output looks entirely ordinary and PCGR will
happily tier it.

The web interface therefore never submits an ambiguous auto-discover run —
it resolves the samples itself and makes you say what they are.

**Validate everything before starting anything.** The hybrid (DNA+RNA)
submission validates *both* halves before queueing either. Before that
fix, a bad RNA half left a DNA run already executing behind an error page
reading "Nothing was started".

**Never address things by position.** Run outputs were linked by their
index in a list of files. But those reports are meant to be readable
*while the run is still going*, so the list grows under an already-rendered
page — and a link labelled "Read quality — SAMPLE_B" started serving
SAMPLE_B's *coverage* report. HTTP 200, wrong file, no error anywhere.

> **Pearl:** if a collection can change between rendering a link and
> following it, the link must address content, not position. This class of
> bug is invisible to "does it return 200?" testing.

**Patient identifiers never touch the command line.** They live in one
JSON file beside the run and are read only when a report is produced. The
pipeline subprocess is never told them, so they cannot end up in a tool's
stdout, a manifest, or a crash trace.

---

## 8. Architecture patterns worth stealing

**Stages that hand off through a manifest.** Each stage writes a JSON file
naming its outputs; the next stage reads it. Not filename guessing — an
explicit, inspectable contract. It also means you can re-run stage 3
without redoing stage 1.

**Share the expensive thing.** The reference genome and its indices live
in one directory that every run resolves against. It used to be a download
and a 1–2 hour index build *per output directory*, for bytes that were
identical.

**Make resumption refuse on drift.** `--resume` reuses finished steps — but
it compares the previous run's parameters first and *aborts* if they
differ. Silently reusing a BQSR BAM built with different known-sites would
produce a call set assembled from two different analyses, and it would
finish successfully.

**Separate environments for conflicting dependencies.** GATK pins a JDK
range; STAR and Arriba do not care; PCGR needs R. Three conda
environments, and the job runner activates the right one per job. One
environment would mean a solver resolving unrelated constraints, and the
thing that gives is the pin nobody is watching.

**Make the expensive optional thing opt-in.** The RNA branch costs ~32 GB
and an hour of CPU. A DNA-only lab should not pay that by accident, so it
installs only with `--with-rna`.

**Idempotent installers.** Every step checks first and skips what is
complete. An interrupted install is resumed by running the same command.

**Configuration as data, with visible precedence.** Panel profiles are
JSON. The rule is `explicit flag > profile > default`, and every run prints
what the profile set *and what it wanted to set but wasn't allowed to*.
The second list is where surprising results usually come from.

**Write down what you did not build.** [CNV_SCOPE.md](CNV_SCOPE.md) and
[RNA_SCOPE.md](RNA_SCOPE.md) describe capabilities that are *designed and
deliberately not implemented*, with what each would cost. This is a real
deliverable: it stops the half-built version appearing later, and it tells
the next person which half was the hard one.

---

## 9. A checklist for your own pipeline

Before you ship something that produces clinical-looking output:

- [ ] For every step: **if this silently did nothing, how would anyone
      know?** Fix the ones with no answer.
- [ ] For every threshold and denominator: where did that number come
      from? Is it right for *this* assay, or for the author's?
- [ ] Does an empty result look different from an unmeasured one?
- [ ] Are expensive checks done *before* expensive work?
- [ ] Does each run record what produced it — settings, versions,
      reference files?
- [ ] Do your reports name what they could **not** measure?
- [ ] Does the interface make the dangerous action harder than the safe
      one?
- [ ] Have you written down what you deliberately did not build?
- [ ] Is there an orthogonal method to check your output against? If not,
      say so loudly in the output itself.

---

## 10. And the honest ending

This pipeline is validated on **one sample, one assay, tumour-only**. The
RNA branch has no validation at all. Every report it produces says so.

That is not modesty, it is the same principle as everything else in this
document: **the system must not claim more than it has earned**, because
the person reading the output cannot tell the difference — and in this
domain, they may act on it.

See the *Read before clinical use* section of [README.md](README.md), and
the scope documents, for exactly what is and is not established.
