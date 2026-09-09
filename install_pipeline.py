#!/usr/bin/env python3
"""
install_pipeline.py
===================
Install the cancer DNA pipeline and everything it needs on a fresh machine.

WHAT THIS IS FOR
----------------
  `install.md` describes the install for a human reader. This script performs
  it, and it exists because doing so by hand went wrong in six separate ways
  on the first attempt -- every one of them silent or misleading, and each is
  now encoded here as a check rather than a paragraph someone has to notice.

  Roughly 70 GB and, on a fast link, about two hours. Most of that is waiting
  on downloads; the bwa-mem2 index is the one long CPU step (~15 minutes).

WHAT IT INSTALLS
----------------
  1. conda environment `cancer_pipeline`  -- every analysis tool
  2. conda environments `pcgr` + `pcgrr`  -- the clinical reporter
  3. hg38 reference + .fai/.dict/bwa-mem2 index
  4. GATK resource VCFs (dbSNP, Mills, known indels, gnomAD, PoN, ExAC)
  5. SnpEff hg38 database
  6. PCGR reference bundle + Ensembl VEP cache
  7. MSIsensor2 models for tumour-only MSI scoring
  8. COSMIC, if you supply the file (see COSMIC below)
  9. the man page

THE TRAPS THIS ENCODES
----------------------
  Each of these cost real time to diagnose. They are checks, not comments.

  * JDK VERSION IS LOAD-BEARING AND `gatk --version` WILL NOT TELL YOU.
    `gatk4` declares only a floor, so a fresh solve installs the newest JDK.
    GATK's Spark tools cannot run on 24+ (JEP 486 made Subject.getSubject()
    throw, and Hadoop calls it while Spark starts), so MarkDuplicatesSpark
    dies AFTER alignment has already run. snpEff 5.4 needs >= 21. The window
    is 21-23 and this script verifies it rather than trusting the solve.

  * PCGR IS NOT ON BIOCONDA. `conda install -c bioconda pcgr` finds nothing.
    It ships from its own channel via version-pinned lock files, and needs
    TWO sibling environments: `pcgr` shells out to `pcgrr`, locating it as
    dirname($CONDA_PREFIX)/pcgrr.

  * TWO BROAD BUCKETS, AND THE FILES ARE SPLIT ACROSS THEM. dbSNP, Mills,
    known indels and the reference live in `gcp-public-data--broad-references`;
    gnomAD, the PoN and small_exac_common_3 live in `gatk-best-practices`.
    Asking either for the other's files returns 403 or 404. The older
    `genomics-public-data` path in circulation is dead to anonymous callers.

  * A TRUNCATED DOWNLOAD IS REUSED FOREVER if the only check is "does the
    file exist". Everything here downloads to a .partial, verifies the size
    against the server's Content-Length, and renames only then.

  * COSMIC SHIPS ENSEMBL CONTIG NAMES (`1`, `MT`) while hg38 here is UCSC
    (`chr1`, `chrM`). bcftools annotate matches on contig name, so the
    unrenamed file annotates NOTHING, silently. This script renames it.

USAGE
-----
  python3 install_pipeline.py --check          # what is missing, change nothing
  python3 install_pipeline.py --dry-run        # print every command, run none
  python3 install_pipeline.py                  # install what is missing
  python3 install_pipeline.py --only envs ref  # just those steps
  python3 install_pipeline.py --cosmic /path/to/Cosmic_*.vcf.gz

  Re-running is safe: every step checks first and skips what is already
  present and complete. An interrupted run is resumed by running it again.

COSMIC
------
  COSMIC cannot be downloaded unattended -- it needs a registered account at
  https://cancer.sanger.ac.uk/cosmic. Download it yourself, then pass the
  file with --cosmic and this script will bgzip it if needed, rename the
  contigs to UCSC style, index it, and verify the result.

REQUIREMENTS
------------
  Python 3.8+, conda or mamba on PATH, wget or curl, and ~70 GB free.
  Nothing else: this script uses only the standard library.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

# =============================================================================
# WHAT GETS INSTALLED
# =============================================================================

# The two Broad buckets. Which file lives in which is not guessable, and
# asking the wrong one gives 403/404 rather than a helpful error.
BROAD_REF = ("https://storage.googleapis.com/"
             "gcp-public-data--broad-references/hg38/v0")
GATK_BP = ("https://storage.googleapis.com/"
           "gatk-best-practices/somatic-hg38")

REFERENCE_FILES = [
    # (url, destination basename). The .fai and .dict are downloaded rather
    # than generated: they are authoritative, and building the dict locally
    # is where a fresh install first fails if GATK is misconfigured.
    (f"{BROAD_REF}/Homo_sapiens_assembly38.fasta",
     "Homo_sapiens_assembly38.fasta"),
    (f"{BROAD_REF}/Homo_sapiens_assembly38.fasta.fai",
     "Homo_sapiens_assembly38.fasta.fai"),
    (f"{BROAD_REF}/Homo_sapiens_assembly38.dict",
     "Homo_sapiens_assembly38.dict"),
]

RESOURCE_FILES = [
    (f"{BROAD_REF}/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
     "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),
    (f"{BROAD_REF}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi",
     "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"),
    (f"{BROAD_REF}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
     "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
    (f"{BROAD_REF}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
     "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"),
    (f"{BROAD_REF}/Homo_sapiens_assembly38.known_indels.vcf.gz",
     "Homo_sapiens_assembly38.known_indels.vcf.gz"),
    (f"{BROAD_REF}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",
     "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"),
    (f"{GATK_BP}/af-only-gnomad.hg38.vcf.gz",
     "af-only-gnomad.hg38.vcf.gz"),
    (f"{GATK_BP}/af-only-gnomad.hg38.vcf.gz.tbi",
     "af-only-gnomad.hg38.vcf.gz.tbi"),
    (f"{GATK_BP}/1000g_pon.hg38.vcf.gz", "1000g_pon.hg38.vcf.gz"),
    (f"{GATK_BP}/1000g_pon.hg38.vcf.gz.tbi", "1000g_pon.hg38.vcf.gz.tbi"),
    (f"{GATK_BP}/small_exac_common_3.hg38.vcf.gz",
     "small_exac_common_3.hg38.vcf.gz"),
    (f"{GATK_BP}/small_exac_common_3.hg38.vcf.gz.tbi",
     "small_exac_common_3.hg38.vcf.gz.tbi"),
]

PCGR_VERSION = "2.3.2"
PCGR_LOCK = (f"https://raw.githubusercontent.com/sigven/pcgr/"
             f"v{PCGR_VERSION}/conda/env/lock")
PCGR_BUNDLE_RELEASE = "20260620"
PCGR_BUNDLE_URL = (f"https://insilico.hpc.uio.no/pcgr/"
                   f"pcgr_ref_data.{PCGR_BUNDLE_RELEASE}.grch38.tgz")

VEP_RELEASE = "115"
VEP_CACHE_URL = (f"https://ftp.ensembl.org/pub/release-{VEP_RELEASE}/"
                 f"variation/indexed_vep_cache/"
                 f"homo_sapiens_vep_{VEP_RELEASE}_GRCh38.tar.gz")

# MSIsensor2's tumour-only mode needs a trained model directory, which the
# conda package does NOT ship -- only the binary. The models live in the
# upstream repository, so the whole tarball is fetched and one directory kept.
MSI_MODELS_URL = ("https://codeload.github.com/niu-lab/msisensor2/"
                  "tar.gz/refs/heads/master")
MSI_MODELS_MEMBER = "msisensor2-master/models_hg38"

PIPELINE_ENV = "cancer_pipeline"
PCGR_ENV = "pcgr"
PCGRR_ENV = "pcgrr"

# Bytes needed, with headroom for the tarballs that are extracted then removed.
DISK_REQUIRED = 75 * 1000 ** 3

# GATK's Spark tools break at 24; snpEff 5.4 needs 21. See module docstring.
JDK_MIN, JDK_MAX = 21, 23


# =============================================================================
# OUTPUT
# =============================================================================

class Log:
    """Console output. Quiet enough to read, loud enough to debug."""

    BOLD, GREEN, YELLOW, RED, DIM, OFF = (
        "\033[1m", "\033[32m", "\033[33m", "\033[31m", "\033[2m", "\033[0m")

    def __init__(self, colour=None):
        self.colour = sys.stdout.isatty() if colour is None else colour

    def _c(self, code, text):
        return f"{code}{text}{self.OFF}" if self.colour else text

    def step(self, text):
        print(f"\n{self._c(self.BOLD, '==> ' + text)}", flush=True)

    def ok(self, text):
        print(f"  {self._c(self.GREEN, 'ok')}    {text}", flush=True)

    def skip(self, text):
        print(f"  {self._c(self.DIM, 'skip')}  {text}", flush=True)

    def warn(self, text):
        print(f"  {self._c(self.YELLOW, 'warn')}  {text}", flush=True)

    def fail(self, text):
        print(f"  {self._c(self.RED, 'FAIL')}  {text}", flush=True)

    def info(self, text):
        print(f"        {text}", flush=True)


LOG = Log()


def human(n):
    """Bytes as a short string."""
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if abs(n) < 1000 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1000.0
    return f"{n:.1f} TB"


# =============================================================================
# SHELL AND DOWNLOAD
# =============================================================================

def run(cmd, dry_run=False, check=True, capture=False, env=None,
        timeout=None):
    """
    Run a command, echoing it first.

    Long steps stream their output rather than buffering it: an installer
    that appears to hang for fifteen minutes gets killed by its operator.
    """
    printable = " ".join(cmd) if isinstance(cmd, list) else cmd
    LOG.info(f"$ {printable}")
    if dry_run:
        return 0, ""
    try:
        if capture:
            res = subprocess.run(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True,
                                 shell=isinstance(cmd, str), env=env,
                                 timeout=timeout)
            if check and res.returncode != 0:
                LOG.fail(f"exit {res.returncode}")
                for line in (res.stdout or "").splitlines()[-15:]:
                    LOG.info(line)
            return res.returncode, res.stdout or ""
        res = subprocess.run(cmd, shell=isinstance(cmd, str), env=env,
                             timeout=timeout)
        if check and res.returncode != 0:
            LOG.fail(f"exit {res.returncode}")
        return res.returncode, ""
    except subprocess.TimeoutExpired:
        LOG.fail(f"timed out after {timeout}s")
        return 124, ""
    except FileNotFoundError as exc:
        LOG.fail(str(exc))
        return 127, ""


def remote_size(url, timeout=60):
    """
    Size the server reports, or None if it cannot be learned.

    None always means "unknown", never "mismatch": a mirror that omits
    Content-Length must not make the install impossible.
    """
    if shutil.which("curl"):
        cmd = ["curl", "-sIL", "--max-time", str(timeout), url]
    elif shutil.which("wget"):
        cmd = ["wget", "--spider", "-S", "--timeout", str(timeout), url]
    else:
        return None
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True,
                             timeout=timeout + 30)
    except Exception:
        return None
    size = None
    for line in res.stdout.splitlines():
        if line.lower().strip().startswith("content-length:"):
            try:
                size = int(line.split(":", 1)[1].strip())
            except ValueError:
                pass
    return size


def download(url, dest, dry_run=False):
    """
    Fetch `url` to `dest` so that `dest` exists only when it is complete.

    Writes to <dest>.partial and renames after the size checks out, because
    the alternative -- writing straight to the final path -- leaves a stump
    behind on any interruption, and every "is it installed?" check in this
    script would then accept it forever.
    """
    if os.path.exists(dest):
        expected = remote_size(url)
        actual = os.path.getsize(dest)
        if expected is None:
            LOG.skip(f"{os.path.basename(dest)} present "
                     f"({human(actual)}, unverified: server unreachable)")
            return True
        if actual == expected:
            LOG.skip(f"{os.path.basename(dest)} present ({human(actual)})")
            return True
        LOG.warn(f"{os.path.basename(dest)} is {human(actual)}, expected "
                 f"{human(expected)} -- re-downloading")

    if not dry_run:
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    partial = dest + ".partial"
    if shutil.which("wget"):
        cmd = ["wget", "-q", "--show-progress", "-c", "-O", partial, url]
    elif shutil.which("curl"):
        cmd = ["curl", "-fL", "--retry", "3", "-o", partial, url]
    else:
        LOG.fail("need wget or curl")
        return False

    LOG.info(f"downloading {os.path.basename(dest)}")
    code, _ = run(cmd, dry_run=dry_run, check=False)
    if dry_run:
        return True
    if code != 0:
        LOG.fail(f"download failed: {url}")
        _discard(partial)
        return False

    actual = os.path.getsize(partial) if os.path.exists(partial) else 0
    if actual == 0:
        LOG.fail(f"empty download: {url}")
        _discard(partial)
        return False
    expected = remote_size(url)
    if expected is not None and actual != expected:
        LOG.fail(f"{os.path.basename(dest)}: got {human(actual)}, server "
                 f"said {human(expected)} -- treating as truncated")
        _discard(partial)
        return False

    os.replace(partial, dest)
    LOG.ok(f"{os.path.basename(dest)} ({human(actual)})")
    return True


def _discard(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        LOG.warn(f"could not remove {path}; delete it by hand")


# =============================================================================
# CONDA
# =============================================================================

def conda_exe():
    """mamba if available -- it solves this environment far faster."""
    return shutil.which("mamba") or shutil.which("conda")


def conda_env_root():
    """Where named environments live, or None if conda cannot be found."""
    exe = conda_exe()
    if not exe:
        return None
    try:
        res = subprocess.run([exe, "info", "--json"], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, timeout=120)
        info = json.loads(res.stdout)
        dirs = info.get("envs_dirs") or []
        return dirs[0] if dirs else None
    except Exception:
        return None


def env_path(name):
    root = conda_env_root()
    return os.path.join(root, name) if root else None


def env_exists(name):
    path = env_path(name)
    return bool(path and os.path.isdir(os.path.join(path, "bin")))


def env_run(name, cmd, **kwargs):
    """
    Run a command inside an environment.

    A properly activated one: bin/ on PATH and CONDA_PREFIX set. Using the
    environment's interpreter alone is not enough -- the tools are found
    with which(), and PCGR reads CONDA_PREFIX to locate its VEP plugins.
    """
    path = env_path(name)
    env = dict(os.environ)
    if path:
        env["PATH"] = os.path.join(path, "bin") + os.pathsep + env.get("PATH", "")
        env["CONDA_PREFIX"] = path
        env["CONDA_DEFAULT_ENV"] = name
        env.pop("PYTHONHOME", None)
        env.pop("PYTHONPATH", None)
    return run(cmd, env=env, **kwargs)


# =============================================================================
# STEPS
# =============================================================================

def step_preflight(args):
    """Check the machine can do this before spending two hours finding out."""
    ok = True

    exe = conda_exe()
    if exe:
        LOG.ok(f"conda: {exe}")
    else:
        LOG.fail("neither mamba nor conda on PATH -- install Miniforge or "
                 "Miniconda first: https://conda-forge.org/download/")
        ok = False

    if shutil.which("wget") or shutil.which("curl"):
        LOG.ok("downloader: " + (shutil.which("wget") or shutil.which("curl")))
    else:
        LOG.fail("need wget or curl")
        ok = False

    if sys.version_info < (3, 8):
        LOG.fail(f"Python {sys.version_info.major}.{sys.version_info.minor}; "
                 f"the pipeline needs 3.8+")
        ok = False
    else:
        LOG.ok(f"python {sys.version_info.major}.{sys.version_info.minor}")

    try:
        free = shutil.disk_usage(args.data_dir if os.path.isdir(args.data_dir)
                                 else os.path.dirname(args.data_dir) or "/").free
        if free < DISK_REQUIRED:
            LOG.fail(f"{human(free)} free, need about "
                     f"{human(DISK_REQUIRED)}")
            ok = False
        else:
            LOG.ok(f"disk: {human(free)} free")
    except OSError as exc:
        LOG.warn(f"could not check free space: {exc}")

    # Reaching one bucket is enough to know the network works; both are
    # checked because the files are split across them and only one being
    # blocked is a confusing failure much later.
    for label, url in (("broad-references", REFERENCE_FILES[1][0]),
                       ("gatk-best-practices", RESOURCE_FILES[-1][0])):
        if remote_size(url):
            LOG.ok(f"reachable: {label}")
        else:
            LOG.warn(f"cannot reach {label} -- downloads may fail")
    return ok


def step_envs(args):
    """The analysis environment, from the repo's environment.yml."""
    if env_exists(PIPELINE_ENV) and not args.force:
        LOG.skip(f"environment '{PIPELINE_ENV}' exists")
        return verify_jdk(args)

    yml = os.path.join(args.repo, "environment.yml")
    if not os.path.exists(yml):
        LOG.fail(f"environment.yml not found at {yml}; use --repo")
        return False

    code, _ = run([conda_exe(), "env", "create", "-f", yml],
                  dry_run=args.dry_run, check=True)
    if code != 0 and not args.dry_run:
        return False
    return verify_jdk(args)


def verify_jdk(args):
    """
    The check that `gatk --version` cannot do for you.

    GATK's Spark tools fail on JDK 24+ AFTER alignment has run, which is an
    expensive way to discover a solve picked the wrong JDK.
    """
    if args.dry_run:
        LOG.info("(dry run) would verify the JDK version")
        return True
    code, out = env_run(PIPELINE_ENV, ["java", "-version"], capture=True,
                        check=False)
    if code != 0:
        LOG.fail("java not found in the environment")
        return False
    first = (out or "").splitlines()[0] if out else ""
    major = None
    for token in first.replace('"', " ").split():
        head = token.split(".")[0]
        if head.isdigit():
            major = int(head)
            break
    if major is None:
        LOG.warn(f"could not parse the Java version from: {first.strip()}")
        return True
    if JDK_MIN <= major <= JDK_MAX:
        LOG.ok(f"JDK {major} (supported range {JDK_MIN}-{JDK_MAX})")
        return True

    LOG.fail(f"JDK {major} is outside the supported range "
             f"{JDK_MIN}-{JDK_MAX}")
    if major > JDK_MAX:
        LOG.info("GATK's Spark tools cannot run on 24+: MarkDuplicatesSpark")
        LOG.info("dies with 'getSubject is not supported' AFTER alignment.")
    else:
        LOG.info("snpEff 5.4 requires JDK 21 or newer.")
    LOG.info(f"fix: conda install -n {PIPELINE_ENV} -c conda-forge "
             f"'openjdk>=21,<24'")
    return False


def step_pcgr_envs(args):
    """PCGR's two environments, from upstream's pinned lock files."""
    ok = True
    for name, lock in ((PCGR_ENV, f"{PCGR_LOCK}/pcgr-linux-64.lock"),
                       (PCGRR_ENV, f"{PCGR_LOCK}/pcgrr-linux-64.lock")):
        if env_exists(name) and not args.force:
            LOG.skip(f"environment '{name}' exists")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            local = os.path.join(tmp, os.path.basename(lock))
            if not download(lock, local, dry_run=args.dry_run):
                LOG.fail(f"could not fetch the {name} lock file")
                ok = False
                continue
            code, _ = run([conda_exe(), "create", "-y", "-n", name,
                           "--file", local], dry_run=args.dry_run)
            if code != 0 and not args.dry_run:
                ok = False

    # Both must be siblings: pcgr finds pcgrr at dirname($CONDA_PREFIX)/pcgrr.
    if not args.dry_run and env_exists(PCGR_ENV) and env_exists(PCGRR_ENV):
        if os.path.dirname(env_path(PCGR_ENV)) != \
                os.path.dirname(env_path(PCGRR_ENV)):
            LOG.fail("'pcgr' and 'pcgrr' are not sibling directories; PCGR "
                     "resolves pcgrr as dirname($CONDA_PREFIX)/pcgrr")
            ok = False
        else:
            LOG.ok("pcgr and pcgrr are siblings")
    return ok


def step_reference(args):
    """hg38 plus its indices. The bwa-mem2 index is the long CPU step."""
    ref_dir = os.path.join(args.data_dir, "references", "hg38")
    os.makedirs(ref_dir, exist_ok=True) if not args.dry_run else None

    for url, name in REFERENCE_FILES:
        if not download(url, os.path.join(ref_dir, name),
                        dry_run=args.dry_run):
            return False

    fasta = os.path.join(ref_dir, "Homo_sapiens_assembly38.fasta")
    bwt = fasta + ".bwt.2bit.64"
    if os.path.exists(bwt) and not args.force:
        LOG.skip(f"bwa-mem2 index present ({human(os.path.getsize(bwt))})")
        return True

    LOG.info("building the bwa-mem2 index -- about 15 minutes, ~20 GB")
    code, _ = env_run(PIPELINE_ENV, ["bwa-mem2", "index", fasta],
                      dry_run=args.dry_run)
    if args.dry_run:
        return True
    if code != 0:
        LOG.fail("bwa-mem2 index failed")
        return False
    LOG.ok("bwa-mem2 index built")
    return True


def step_resources(args):
    """The GATK resource VCFs, from whichever bucket actually holds them."""
    res_dir = os.path.join(args.data_dir, "resources", "hg38")
    if not args.dry_run:
        os.makedirs(res_dir, exist_ok=True)
    ok = True
    for url, name in RESOURCE_FILES:
        if not download(url, os.path.join(res_dir, name),
                        dry_run=args.dry_run):
            ok = False
    return ok


def step_snpeff(args):
    """
    The SnpEff hg38 database.

    Its location is decided by snpEff's own config, so the database is
    located by asking snpEff rather than by guessing a path that differs
    between conda-forge and bioconda builds.
    """
    if args.dry_run:
        LOG.info("(dry run) would run: snpEff download hg38")
        return True
    code, out = env_run(PIPELINE_ENV, ["snpEff", "databases"], capture=True,
                        check=False, timeout=600)
    installed = False
    if code == 0:
        share = os.path.join(env_path(PIPELINE_ENV) or "", "share")
        for root, dirs, _files in os.walk(share):
            if os.path.basename(root) == "hg38" and \
                    os.path.basename(os.path.dirname(root)) == "data":
                LOG.skip(f"SnpEff hg38 database present ({root})")
                installed = True
                break
    if installed and not args.force:
        return True

    LOG.info("downloading the SnpEff hg38 database (~450 MB)")
    code, _ = env_run(PIPELINE_ENV, ["snpEff", "download", "hg38"],
                      check=False, timeout=3600)
    if code != 0:
        LOG.fail("snpEff download hg38 failed")
        return False
    LOG.ok("SnpEff hg38 database installed")
    return True


def _extract_tgz(url, dest_dir, marker, label, args):
    """
    Download a tarball, extract it, and delete the archive.

    Streamed through tar rather than read into memory: these are 7-24 GB.
    The marker is what proves the extraction finished, so an interrupted
    run redoes it instead of leaving a half-tree that looks installed.
    """
    if os.path.exists(marker) and not args.force:
        LOG.skip(f"{label} present ({marker})")
        return True
    if not args.dry_run:
        os.makedirs(dest_dir, exist_ok=True)
    archive = os.path.join(dest_dir, os.path.basename(url))
    if not download(url, archive, dry_run=args.dry_run):
        return False

    LOG.info(f"extracting {label} -- this takes a while and needs the space "
             f"twice over until the archive is removed")
    code, _ = run(f"gzip -dc {archive!r} | tar xf - -C {dest_dir!r}",
                  dry_run=args.dry_run, check=True)
    if code != 0 and not args.dry_run:
        LOG.fail(f"could not extract {label}")
        return False
    if args.dry_run:
        return True
    _discard(archive)
    if not os.path.exists(marker):
        LOG.fail(f"{label} extracted but {marker} is missing")
        return False
    LOG.ok(f"{label} installed")
    return True


def step_pcgr_data(args):
    """PCGR's reference bundle and the Ensembl VEP cache."""
    pcgr_root = os.path.join(args.data_dir, "pcgr")
    release_dir = os.path.join(pcgr_root, PCGR_BUNDLE_RELEASE)
    marker = os.path.join(release_dir, "data", "grch38")

    ok = True
    if os.path.exists(marker) and not args.force:
        LOG.skip(f"PCGR bundle present ({marker})")
    else:
        # The tarball unpacks to ./data, which must end up under a directory
        # named for the release -- that name is how PCGR is pointed at it.
        if not args.dry_run:
            os.makedirs(release_dir, exist_ok=True)
        if not _extract_tgz(PCGR_BUNDLE_URL, release_dir,
                            marker, "PCGR reference bundle", args):
            ok = False

    vep_dir = os.path.join(args.data_dir, "vep_cache")
    vep_marker = os.path.join(vep_dir, "homo_sapiens",
                              f"{VEP_RELEASE}_GRCh38")
    if not _extract_tgz(VEP_CACHE_URL, vep_dir, vep_marker,
                        f"Ensembl VEP {VEP_RELEASE} cache", args):
        ok = False
    return ok


def step_msi_models(args):
    """
    Install the MSIsensor2 models for tumour-only MSI scoring.

    The conda package ships the binary only. Without these models the
    tumour-only mode cannot run at all, and PCGR does not cover the gap:
    it restricts MSI to WGS/WES tumour-control runs and omits the section
    silently on anything else.

    The upstream project publishes the models inside its git repository
    rather than as a release asset, so the source tarball is fetched and a
    single directory extracted from it. About 310 MB down, 251 MB kept.
    """
    dest = os.path.join(args.data_dir, "msisensor2")
    marker = os.path.join(dest, "models_hg38")
    if os.path.isdir(marker) and not args.force:
        count = len(os.listdir(marker))
        LOG.skip(f"MSI models present ({count} files)")
        return True
    if args.dry_run:
        LOG.info(f"(dry run) would fetch {MSI_MODELS_URL}")
        LOG.info(f"(dry run) would extract {MSI_MODELS_MEMBER} -> {marker}")
        return True

    os.makedirs(dest, exist_ok=True)
    archive = os.path.join(dest, "msisensor2-master.tar.gz")
    if not download(MSI_MODELS_URL, archive):
        LOG.fail("could not fetch the MSIsensor2 models")
        return False

    LOG.info("extracting the hg38 models")
    code, _ = run(["tar", "xzf", archive, "-C", dest,
                   "--strip-components=1", MSI_MODELS_MEMBER], check=True)
    _discard(archive)
    if code != 0 or not os.path.isdir(marker):
        LOG.fail("could not extract the MSIsensor2 models")
        return False
    LOG.ok(f"MSI models installed ({len(os.listdir(marker))} files)")
    return True


def step_cosmic(args):
    """
    Import a COSMIC VCF the operator downloaded, and fix its contig names.

    COSMIC needs a registered account, so it cannot be fetched here. What
    this does is the part that is easy to get wrong: COSMIC uses Ensembl
    contig names (1, MT) while this reference is UCSC (chr1, chrM), and
    bcftools annotate matches on the name. Unrenamed, COSMIC annotation
    silently matches nothing at all and the run looks entirely successful.
    """
    if not args.cosmic:
        LOG.skip("no --cosmic given; COSMIC annotation will be unavailable")
        LOG.info("register at https://cancer.sanger.ac.uk/cosmic, download "
                 "the GRCh38 VCF, then re-run with --cosmic <file>")
        return True

    src = os.path.abspath(args.cosmic)
    if not os.path.exists(src):
        LOG.fail(f"--cosmic file not found: {src}")
        return False

    res_dir = os.path.join(args.data_dir, "resources", "hg38")
    if not args.dry_run:
        os.makedirs(res_dir, exist_ok=True)
    base = os.path.basename(src)
    for suffix in (".vcf.gz", ".vcf"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    out = os.path.join(res_dir, base + ".chr.vcf.gz")

    if os.path.exists(out) and not args.force:
        LOG.skip(f"renamed COSMIC present ({out})")
        return True

    if args.dry_run:
        LOG.info(f"(dry run) would rename contigs of {src} -> {out}")
        return True

    # 1..22, X, Y and MT->chrM. Contigs absent from the map are dropped by
    # bcftools, which is what we want: COSMIC's scaffolds have no counterpart
    # in this reference anyway.
    chr_map = os.path.join(res_dir, ".cosmic_chr_map.txt")
    with open(chr_map, "w", encoding="utf-8") as fh:
        for c in list(range(1, 23)) + ["X", "Y"]:
            fh.write(f"{c} chr{c}\n")
        fh.write("MT chrM\n")

    LOG.info("renaming COSMIC contigs to UCSC style")
    code, _ = env_run(PIPELINE_ENV,
                      ["bcftools", "annotate", "--rename-chrs", chr_map,
                       "--threads", str(args.threads), "-Oz", "-o", out, src],
                      check=True)
    if code != 0:
        LOG.fail("bcftools annotate failed")
        return False
    code, _ = env_run(PIPELINE_ENV, ["tabix", "-f", "-p", "vcf", out])
    if code != 0:
        LOG.fail("tabix failed")
        return False

    # Prove the rename worked rather than assuming: the whole point is that
    # the failure mode is silent.
    code, out_txt = env_run(PIPELINE_ENV,
                            ["bcftools", "view", "-h", out],
                            capture=True, check=False)
    if code == 0 and "##contig=<ID=chr" not in (out_txt or ""):
        LOG.fail("renamed COSMIC has no chr-prefixed contigs; check the input")
        return False
    LOG.ok(f"COSMIC installed and renamed ({out})")
    return True


def step_manpage(args):
    """The man page, into the user's own man path."""
    src = os.path.join(args.repo, "cancer-dna-pipeline.1")
    if not os.path.exists(src):
        LOG.skip("no man page in the repo")
        return True
    dest_dir = os.path.expanduser("~/.local/share/man/man1")
    dest = os.path.join(dest_dir, "cancer-dna-pipeline.1")
    if args.dry_run:
        LOG.info(f"(dry run) would install {src} -> {dest}")
        return True
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copyfile(src, dest)
    os.chmod(dest, 0o644)
    LOG.ok(f"man page installed ({dest})")
    LOG.info("read it with: man cancer-dna-pipeline")
    return True


# =============================================================================
# VERIFICATION
# =============================================================================

def step_verify(args):
    """
    Prove the install, rather than reporting that the steps ran.

    Checks what the pipeline actually needs at runtime: the tools resolve
    inside their environment, the JDK is in range, indices sit beside their
    data, and the two PCGR environments are siblings.
    """
    if args.dry_run:
        LOG.info("(dry run) would verify the installation")
        return True

    ok = True
    tools = ["fastp", "bwa-mem2", "samtools", "gatk", "bcftools", "snpEff",
             "msisensor2"]
    missing = []
    for tool in tools:
        code, _ = env_run(PIPELINE_ENV, ["bash", "-lc",
                                         f"command -v {tool}"],
                          capture=True, check=False)
        if code != 0:
            missing.append(tool)
    if missing:
        LOG.fail(f"not on PATH in '{PIPELINE_ENV}': {', '.join(missing)}")
        ok = False
    else:
        LOG.ok(f"all {len(tools)} analysis tools resolve")

    if not verify_jdk(args):
        ok = False

    ref_dir = os.path.join(args.data_dir, "references", "hg38")
    fasta = os.path.join(ref_dir, "Homo_sapiens_assembly38.fasta")
    for suffix, what in ((".fai", "samtools index"),
                         (".bwt.2bit.64", "bwa-mem2 index")):
        if os.path.exists(fasta + suffix):
            LOG.ok(f"{what} present")
        else:
            LOG.fail(f"{what} missing ({fasta + suffix})")
            ok = False
    dict_file = os.path.splitext(fasta)[0] + ".dict"
    if os.path.exists(dict_file):
        LOG.ok("GATK sequence dictionary present")
    else:
        LOG.fail(f"sequence dictionary missing ({dict_file})")
        ok = False

    res_dir = os.path.join(args.data_dir, "resources", "hg38")
    for _url, name in RESOURCE_FILES:
        if name.endswith(".tbi"):
            continue
        path = os.path.join(res_dir, name)
        if not os.path.exists(path):
            LOG.fail(f"missing resource: {name}")
            ok = False
        elif not os.path.exists(path + ".tbi"):
            LOG.fail(f"{name} has no .tbi beside it -- GATK and bcftools "
                     f"never search for indexes elsewhere")
            ok = False
    if ok:
        LOG.ok("resource VCFs present and indexed")

    # The clinical report is the deliverable, not an optional extra: a run
    # that ends at an annotated VCF is not a finished run. Anything PCGR
    # needs is therefore a failure here, not a warning.
    for name in (PCGR_ENV, PCGRR_ENV):
        if env_exists(name):
            LOG.ok(f"environment '{name}' present")
        else:
            LOG.fail(f"environment '{name}' missing -- no clinical reports")
            ok = False

    bundle = os.path.join(args.data_dir, "pcgr", PCGR_BUNDLE_RELEASE,
                          "data", "grch38")
    vep = os.path.join(args.data_dir, "vep_cache", "homo_sapiens",
                       f"{VEP_RELEASE}_GRCh38")
    for path, label in ((bundle, "PCGR bundle"), (vep, "VEP cache")):
        if os.path.isdir(path):
            LOG.ok(f"{label} present")
        else:
            LOG.fail(f"{label} missing ({path}) -- no clinical reports")
            ok = False

    msi_models = os.path.join(args.data_dir, "msisensor2", "models_hg38")
    if os.path.isdir(msi_models):
        LOG.ok("MSIsensor2 models present")
    else:
        LOG.warn(f"MSI models missing ({msi_models}) -- step 8 will be "
                 f"skipped. PCGR does not cover MSI on a panel.")

    # COSMIC stays a warning: it cannot be downloaded without an account,
    # so its absence is a licensing fact rather than a broken install.
    cosmic = [f for f in os.listdir(res_dir)
              if f.endswith(".chr.vcf.gz")] if os.path.isdir(res_dir) else []
    if cosmic:
        LOG.ok(f"COSMIC present ({cosmic[0]})")
    else:
        LOG.warn("no chr-renamed COSMIC -- step 10 will be skipped. Supply "
                 "it with --cosmic <file>")

    # A dry run of the pipeline itself is the end-to-end check: it builds
    # every command line without touching data or the network.
    script = os.path.join(args.repo, "comprehensive_variant_calling.py")
    if os.path.exists(script):
        code, _ = env_run(PIPELINE_ENV, ["python", script, "--help"],
                          capture=True, check=False)
        if code == 0:
            LOG.ok("pipeline script runs")
        else:
            LOG.fail("the pipeline script failed to start")
            ok = False

    # The pipeline imports these from beside itself, late in a run: PCGR at
    # step 13 and the coverage check at step 12. A bundle missing one fails
    # after the hours, not before them.
    for helper, what in (("pcgr_report.py", "clinical report"),
                         ("coverage_report.py", "coverage check")):
        path = os.path.join(args.repo, helper)
        if not os.path.exists(path):
            LOG.fail(f"{helper} missing -- the {what} cannot run")
            ok = False
            continue
        code, _ = env_run(PIPELINE_ENV, ["python", path, "--help"],
                          capture=True, check=False)
        if code == 0:
            LOG.ok(f"{helper} runs")
        else:
            LOG.fail(f"{helper} failed to start")
            ok = False

    # The pipeline looks for a named genome under ~/data/references unless
    # told otherwise. With a moved --data-dir it would not find the one just
    # installed, and would download a second copy instead.
    default_ref_dir = os.path.expanduser("~/data/references")
    installed_ref_dir = os.path.join(args.data_dir, "references")
    if os.path.abspath(installed_ref_dir) != os.path.abspath(default_ref_dir):
        LOG.warn(
            f"the genome is in {installed_ref_dir}, but the pipeline looks in "
            f"{default_ref_dir} by default. Export "
            f"PIPELINE_REFERENCE_DIR={installed_ref_dir} (or pass "
            f"--reference-dir) so '--reference hg38' finds this copy instead "
            f"of downloading and indexing another.")
    else:
        LOG.ok("the pipeline will find this genome by name (--reference hg38)")
    return ok


# =============================================================================
# DRIVER
# =============================================================================

STEPS = [
    ("preflight", "Check conda, downloader, disk and network", step_preflight),
    ("envs", "Create the cancer_pipeline environment", step_envs),
    ("pcgr-envs", "Create the pcgr and pcgrr environments", step_pcgr_envs),
    ("ref", "hg38 reference and indices", step_reference),
    ("resources", "GATK resource VCFs", step_resources),
    ("snpeff", "SnpEff hg38 database", step_snpeff),
    ("pcgr-data", "PCGR bundle and VEP cache", step_pcgr_data),
    ("msi-models", "MSIsensor2 models for tumour-only MSI", step_msi_models),
    ("cosmic", "Import and rename a COSMIC VCF", step_cosmic),
    ("man", "Install the man page", step_manpage),
    ("verify", "Verify the installation", step_verify),
]


def build_parser():
    names = ", ".join(name for name, _d, _f in STEPS)
    parser = argparse.ArgumentParser(
        description="Install the cancer DNA pipeline and its dependencies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Steps, in order: {names}\n\n"
               f"Re-running is safe: each step checks first and skips what is "
               f"already complete.\n")
    parser.add_argument("--data-dir", default=os.path.expanduser("~/data"),
                        help="Where references and databases go "
                             "(default: %(default)s).")
    parser.add_argument("--repo", default=os.path.dirname(
        os.path.abspath(__file__)),
        help="The pipeline checkout, for environment.yml and the man page "
             "(default: this script's directory).")
    parser.add_argument("--cosmic", default=None, metavar="VCF",
                        help="A COSMIC GRCh38 VCF you downloaded. It is "
                             "copied in, its contigs renamed to UCSC style "
                             "and indexed. Without this, COSMIC annotation "
                             "is unavailable.")
    parser.add_argument("--only", nargs="+", metavar="STEP", default=None,
                        help="Run only these steps.")
    parser.add_argument("--skip", nargs="+", metavar="STEP", default=[],
                        help="Skip these steps.")
    parser.add_argument("--check", action="store_true",
                        help="Report what is missing and change nothing. "
                             "Equivalent to --only verify.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print every command without running any.")
    parser.add_argument("--force", action="store_true",
                        help="Redo steps whose output already exists. Will "
                             "re-download tens of gigabytes.")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads for indexing and COSMIC conversion "
                             "(default: %(default)s).")
    parser.add_argument("--no-colour", action="store_true",
                        help="Plain output, for logs and CI.")
    return parser


def main():
    args = build_parser().parse_args()
    if args.no_colour:
        global LOG
        LOG = Log(colour=False)

    if args.check:
        args.only = ["verify"]

    valid = {name for name, _d, _f in STEPS}
    for group in (args.only or [], args.skip):
        unknown = set(group) - valid
        if unknown:
            print(f"Unknown step(s): {', '.join(sorted(unknown))}\n"
                  f"Valid: {', '.join(sorted(valid))}", file=sys.stderr)
            return 2

    selected = [s for s in STEPS
                if (args.only is None or s[0] in args.only)
                and s[0] not in args.skip]

    print(f"{LOG._c(Log.BOLD, 'Cancer DNA pipeline installer')}")
    print(f"  repo     : {args.repo}")
    print(f"  data dir : {args.data_dir}")
    print(f"  steps    : {', '.join(s[0] for s in selected)}")
    if args.dry_run:
        print(f"  {LOG._c(Log.YELLOW, 'DRY RUN -- nothing will change')}")

    started = time.time()
    failed = []
    for name, description, func in selected:
        LOG.step(f"{name}: {description}")
        try:
            if not func(args):
                failed.append(name)
                # Everything downstream needs the environment, so there is
                # no value in continuing without it.
                if name in ("preflight", "envs"):
                    LOG.fail(f"'{name}' is required; stopping here")
                    break
        except KeyboardInterrupt:
            print("\nInterrupted. Re-run to resume: completed steps are "
                  "skipped.", file=sys.stderr)
            return 130
        except Exception as exc:                # noqa: BLE001 - reported
            LOG.fail(f"{type(exc).__name__}: {exc}")
            failed.append(name)

    mins = (time.time() - started) / 60
    print()
    if failed:
        print(LOG._c(Log.RED, f"Failed: {', '.join(failed)}  ({mins:.1f} min)"))
        print("Fix the causes above and re-run; finished steps are skipped.")
        return 1
    print(LOG._c(Log.GREEN, f"All steps completed ({mins:.1f} min)"))
    if not args.dry_run and (args.only is None or "verify" in (args.only or [])):
        print("\nNext:")
        print(f"  conda activate {PIPELINE_ENV}")
        if os.path.abspath(args.data_dir) != os.path.expanduser("~/data"):
            print(f"  export PIPELINE_REFERENCE_DIR="
                  f"{os.path.join(args.data_dir, 'references')}"
                  f"   # or the pipeline looks in ~/data/references")
        print(f"  python {os.path.join(args.repo, 'webapp', 'app.py')} "
              f"--allow-root {args.data_dir} --allow-root /path/to/fastqs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
