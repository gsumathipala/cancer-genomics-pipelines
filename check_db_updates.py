#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
check_db_updates.py
===================
Ask whether any installed reference database has a newer release, and write
the answer to a file the web interface reads.

WHAT IT WILL NOT DO
-------------------
  This is a notifier, not an updater. It is deliberately incapable of
  changing anything a run depends on:

    * It never writes to --data-dir. Every check is an HTTP HEAD or a
      read-only API call; nothing is downloaded and nothing is replaced.
    * It never touches a conda environment.
    * Its only output is one small JSON status file, outside the data tree.
    * It exits 0 even when checks fail, so a cron entry cannot start
      mailing you about a flaky network.

  Upgrading a database mid-project is a decision with consequences -- a new
  VEP cache changes annotations, a new COSMIC changes IDs, and either can
  make today's report disagree with last week's for reasons that have
  nothing to do with the sample. So this tells you, and stops.

  Acting on it is install_pipeline.py's job. The web interface runs this
  script at startup and on demand, and offers a per-database Update button
  that runs the installer for that one step -- but the decision is still
  taken by a person, in front of a warning, and never while a run is in
  progress. Nothing in this file changes as a result: it still only reports.

WHAT IT CHECKS
--------------
  Ensembl VEP cache      release directory listing on ftp.ensembl.org
  PCGR                   latest release tag on GitHub
  MSIsensor2 models      last commit touching models_hg38
  SnpEff                 latest version on the bioconda channel
  GATK resource VCFs     remote size vs the local file (drift/truncation)
  COSMIC                 reported as MANUAL: it needs an account, so no
                         unattended check is possible. The installed
                         version is shown next to the release-notes URL.

  A source that cannot be reached is reported as "unknown", never as
  "up to date". Silence and good news must not look the same.

USAGE
-----
  python3 check_db_updates.py                  # check, write the status file
  python3 check_db_updates.py --print          # human-readable, to stdout
  python3 check_db_updates.py --json           # machine-readable, to stdout

  Weekly, via cron -- the leading hyphen keeps it quiet on success:
    17 6 * * 1  /usr/bin/python3 /path/to/check_db_updates.py

REQUIREMENTS
------------
  Python 3.8+, and wget or curl. Standard library only.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

# The status file lives outside the data tree, so this tool has no reason to
# hold a write handle anywhere the pipeline reads.
DEFAULT_STATUS = os.path.expanduser("~/.cache/cancer_pipeline/db_updates.json")
DEFAULT_DATA = os.path.expanduser("~/data")

# How long before a status file is treated as stale rather than trusted.
STALE_AFTER_DAYS = 30

HTTP_TIMEOUT = 30


# =============================================================================
# FETCHING
# =============================================================================

def fetch(url, timeout=HTTP_TIMEOUT):
    """
    GET a URL and return its body as text, or None.

    None means "could not ask", which callers must not confuse with "nothing
    has changed" -- an unreachable source is reported as unknown.
    """
    if shutil.which("curl"):
        cmd = ["curl", "-sSL", "--max-time", str(timeout), url]
    elif shutil.which("wget"):
        cmd = ["wget", "-qO-", "--timeout", str(timeout), url]
    else:
        return None
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True,
                             timeout=timeout + 15)
        return res.stdout if res.returncode == 0 else None
    except Exception:
        return None


def fetch_json(url, timeout=HTTP_TIMEOUT):
    """GET and parse JSON, or None."""
    body = fetch(url, timeout)
    if not body:
        return None
    try:
        return json.loads(body)
    except ValueError:
        return None


def remote_size(url, timeout=HTTP_TIMEOUT):
    """Content-Length for a URL, or None when the server will not say."""
    if shutil.which("curl"):
        cmd = ["curl", "-sIL", "--max-time", str(timeout), url]
    elif shutil.which("wget"):
        cmd = ["wget", "--spider", "-S", "--timeout", str(timeout), url]
    else:
        return None
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True,
                             timeout=timeout + 15)
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


# =============================================================================
# ONE RESULT PER DATABASE
# =============================================================================

def result(name, installed, latest=None, state="unknown", note="", url=""):
    """
    One database's answer.

    state is one of:
      current  -- installed matches the newest available
      update   -- something newer exists
      manual   -- cannot be checked unattended; a human must look
      unknown  -- the source could not be reached
      changed  -- the remote file differs from the local one (see GATK)
    """
    return {"name": name, "installed": installed, "latest": latest,
            "state": state, "note": note, "url": url}


def check_vep(data_dir):
    """
    Ensembl VEP cache: the release directories on the FTP server.

    Upgrading is a real decision, not housekeeping -- a new VEP release can
    change transcript sets and therefore consequence calls, so a report
    regenerated after an upgrade may legitimately differ from one issued
    before it.
    """
    cache = os.path.join(data_dir, "vep_cache", "homo_sapiens")
    installed = None
    if os.path.isdir(cache):
        rels = [int(m.group(1)) for d in os.listdir(cache)
                for m in [re.match(r"(\d+)_GRCh38$", d)] if m]
        installed = max(rels) if rels else None

    body = fetch("https://ftp.ensembl.org/pub/")
    if not body:
        return result("Ensembl VEP cache", installed, state="unknown",
                      note="ftp.ensembl.org unreachable",
                      url="https://ftp.ensembl.org/pub/")
    found = sorted({int(m) for m in re.findall(r"release-(\d+)", body)})
    latest = found[-1] if found else None
    if installed is None:
        return result("Ensembl VEP cache", "not installed", latest,
                      "update", "no cache found under vep_cache/",
                      "https://ftp.ensembl.org/pub/")
    if latest and latest > installed:
        return result(
            "Ensembl VEP cache", installed, latest, "update",
            f"release {latest} available. Upgrading changes the transcript "
            f"set, so annotations may shift; re-run affected samples rather "
            f"than mixing releases in one cohort.",
            f"https://ftp.ensembl.org/pub/release-{latest}/variation/"
            f"indexed_vep_cache/")
    return result("Ensembl VEP cache", installed, latest, "current")


def check_pcgr(data_dir):
    """
    PCGR, via its GitHub release tag.

    The reference bundle itself cannot be checked: its directory is not
    listable, and the filename carries a build date with no index. A new
    PCGR release is the usable proxy, since the two are versioned together
    and a mismatched pair is the classic failed first run.
    """
    installed_env = None
    envs = os.path.expanduser("~/miniconda3/envs/pcgr/bin/pcgr")
    if os.path.exists(envs):
        try:
            out = subprocess.run([envs, "--version"], stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True,
                                 timeout=120).stdout
            m = re.search(r"(\d+\.\d+\.\d+)", out or "")
            installed_env = m.group(1) if m else None
        except Exception:
            pass

    bundles = []
    root = os.path.join(data_dir, "pcgr")
    if os.path.isdir(root):
        bundles = sorted(d for d in os.listdir(root) if re.fullmatch(r"\d{8}", d))
    bundle = bundles[-1] if bundles else "none"

    data = fetch_json("https://api.github.com/repos/sigven/pcgr/releases/latest")
    latest = (data or {}).get("tag_name", "").lstrip("v") or None
    installed = f"{installed_env or '?'} (bundle {bundle})"
    if not latest:
        return result("PCGR", installed, state="unknown",
                      note="GitHub API unreachable",
                      url="https://github.com/sigven/pcgr/releases")
    if installed_env and latest != installed_env:
        return result(
            "PCGR", installed, latest, "update",
            f"PCGR {latest} released. The reference bundle is versioned with "
            f"the software -- upgrade both together, or the run fails on a "
            f"version mismatch.",
            "https://github.com/sigven/pcgr/releases")
    return result("PCGR", installed, latest, "current",
                  "bundle release cannot be checked automatically; it is "
                  "published without a directory listing")


def check_msi_models(data_dir):
    """MSIsensor2 models: the last commit that touched models_hg38."""
    local = os.path.join(data_dir, "msisensor2", "models_hg38")
    installed = f"{len(os.listdir(local))} files" if os.path.isdir(local) \
        else "not installed"
    data = fetch_json("https://api.github.com/repos/niu-lab/msisensor2/"
                      "commits?path=models_hg38&per_page=1")
    if not data:
        return result("MSIsensor2 models", installed, state="unknown",
                      note="GitHub API unreachable",
                      url="https://github.com/niu-lab/msisensor2")
    try:
        sha = data[0]["sha"][:12]
        date = data[0]["commit"]["committer"]["date"][:10]
    except (KeyError, IndexError, TypeError):
        return result("MSIsensor2 models", installed, state="unknown",
                      note="unexpected API response")
    if installed == "not installed":
        return result("MSIsensor2 models", installed, f"{sha} ({date})",
                      "update", "no models installed; MSI scoring is skipped",
                      "https://github.com/niu-lab/msisensor2")
    return result("MSIsensor2 models", installed, f"{sha} ({date})", "current",
                  f"upstream models last changed {date}")


def check_snpeff():
    """SnpEff, via the bioconda channel."""
    env = os.path.expanduser("~/miniconda3/envs/cancer_pipeline/bin/snpEff")
    installed = None
    if os.path.exists(env):
        try:
            out = subprocess.run([env, "-version"], stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True,
                                 timeout=180).stdout
            m = re.search(r"([0-9]+\.[0-9]+[a-z]?)", out or "")
            installed = m.group(1) if m else None
        except Exception:
            pass
    data = fetch_json("https://api.anaconda.org/package/bioconda/snpeff")
    latest = (data or {}).get("latest_version")
    if not latest:
        return result("SnpEff", installed or "?", state="unknown",
                      note="anaconda.org unreachable")
    # 5.4c and 5.4.0c are the same build numbered differently by channel.
    norm = lambda v: re.sub(r"[^0-9a-z]", "", (v or "").lower())
    if installed and norm(latest).startswith(norm(installed)[:3]) and \
            norm(installed) not in norm(latest) and \
            norm(latest) not in norm(installed):
        state = "update"
    else:
        state = "current"
    return result("SnpEff", installed or "?", latest, state,
                  "" if state == "current" else
                  "a newer build exists; the hg38 database must be "
                  "re-downloaded after upgrading",
                  "https://anaconda.org/bioconda/snpeff")


def check_gatk_resources(data_dir):
    """
    The GATK resource VCFs.

    These are effectively frozen, so this is not really an update check --
    it is a drift check. A local file whose size no longer matches the
    bucket means one of them changed underneath you, or the local copy is
    truncated. Either is worth knowing before a run trusts it.
    """
    broad = ("https://storage.googleapis.com/"
             "gcp-public-data--broad-references/hg38/v0")
    bp = "https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
    files = [
        (f"{broad}/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
         "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),
        (f"{bp}/af-only-gnomad.hg38.vcf.gz", "af-only-gnomad.hg38.vcf.gz"),
        (f"{bp}/1000g_pon.hg38.vcf.gz", "1000g_pon.hg38.vcf.gz"),
        (f"{bp}/small_exac_common_3.hg38.vcf.gz",
         "small_exac_common_3.hg38.vcf.gz"),
    ]
    res_dir = os.path.join(data_dir, "resources", "hg38")
    checked, drifted, unknown = 0, [], 0
    for url, name in files:
        path = os.path.join(res_dir, name)
        if not os.path.exists(path):
            drifted.append(f"{name}: missing locally")
            continue
        remote = remote_size(url)
        if remote is None:
            unknown += 1
            continue
        checked += 1
        local = os.path.getsize(path)
        if local != remote:
            drifted.append(f"{name}: local {local:,} vs remote {remote:,}")
    if drifted:
        return result("GATK resource VCFs", f"{checked} verified", None,
                      "changed", "; ".join(drifted))
    if checked == 0:
        return result("GATK resource VCFs", "?", None, "unknown",
                      "no file could be checked")
    return result("GATK resource VCFs", f"{checked} of {len(files)} verified",
                  None, "current",
                  f"{unknown} could not be reached" if unknown else "")


def check_cosmic(data_dir):
    """
    COSMIC: reported as manual, because it genuinely is.

    Downloads need a registered account, so there is no unattended check
    that would not amount to scraping a login-walled page and guessing.
    Saying "manual" is more useful than a check that quietly fails open.
    """
    res_dir = os.path.join(data_dir, "resources", "hg38")
    installed = "not installed"
    if os.path.isdir(res_dir):
        found = sorted(f for f in os.listdir(res_dir)
                       if f.startswith("Cosmic") and f.endswith(".chr.vcf.gz"))
        if found:
            m = re.search(r"_v(\d+)_", found[-1])
            installed = f"v{m.group(1)}" if m else found[-1]
    return result("COSMIC", installed, None, "manual",
                  "needs a registered account; check the release notes and "
                  "re-import with install_pipeline.py --cosmic",
                  "https://cancer.sanger.ac.uk/cosmic/release_notes")


# =============================================================================
# =============================================================================
# THE RNA BRANCH'S REFERENCE DATA
# =============================================================================
# Only reported when it is installed. A DNA-only laboratory has not chosen
# to have these and should not be told they are out of date.

def check_gencode(data_dir):
    """
    GENCODE annotation, via the release directories on the EBI FTP server.

    Upgrading is the same kind of decision a VEP upgrade is, with one extra
    consequence: the annotation is BAKED INTO the STAR index. Changing
    GENCODE without rebuilding the index leaves the aligner finding
    junctions from one annotation while the fusion caller names genes from
    another -- see check_star_index() below, which is the check that
    actually catches it.
    """
    directory = os.path.join(data_dir, "references", "gencode")
    installed = None
    if os.path.isdir(directory):
        releases = [int(m.group(1)) for f in os.listdir(directory)
                    for m in [re.match(r"gencode\.v(\d+)\.", f)] if m]
        installed = max(releases) if releases else None
    if installed is None:
        return None                       # the RNA branch is not installed

    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    body = fetch(url)
    if not body:
        return result("GENCODE annotation (RNA)", installed, state="unknown",
                      note="ftp.ebi.ac.uk unreachable", url=url)
    found = sorted({int(m) for m in re.findall(r"release_(\d+)", body)})
    latest = found[-1] if found else None
    if latest and latest > installed:
        return result(
            "GENCODE annotation (RNA)", installed, latest, "update",
            f"release {latest} available. The annotation is baked into the "
            f"STAR index, so upgrading means rebuilding it "
            f"(install_pipeline.py --only gencode star-index --force) -- "
            f"about an hour. Gene names and transcript sets shift between "
            f"releases, so do not mix them within a cohort.",
            f"{url}release_{latest}/")
    return result("GENCODE annotation (RNA)", installed, latest, "current")


def check_star_index(data_dir):
    """
    Is each STAR index still consistent with the installed annotation?

    This is the check that matters, and nothing else performs it. An index
    outliving the GTF it was built from is completely silent: STAR runs,
    the mapping rate looks normal, and the junctions it knows about are the
    old ones. align_rna.py warns at run time when it can compare the two,
    but only a run does that -- this asks the question without one.
    """
    references = os.path.join(data_dir, "references")
    if not os.path.isdir(references):
        return None
    indexes = [d for d in sorted(os.listdir(references))
               if d.startswith("star_")
               and os.path.exists(os.path.join(references, d, "SAindex"))]
    if not indexes:
        return None                       # the RNA branch is not installed

    gencode_dir = os.path.join(references, "gencode")
    installed_gtfs = set()
    if os.path.isdir(gencode_dir):
        installed_gtfs = {f for f in os.listdir(gencode_dir)
                          if f.endswith(".gtf")}

    stale = []
    unknown = []
    for name in indexes:
        record_path = os.path.join(references, name, "pipeline_index.json")
        try:
            with open(record_path, encoding="utf-8") as fh:
                record = json.load(fh)
        except (OSError, ValueError):
            unknown.append(name)
            continue
        built_from = os.path.basename(record.get("gtf") or "")
        if installed_gtfs and built_from and built_from not in installed_gtfs:
            stale.append(f"{name} (built from {built_from})")

    summary = f"{len(indexes)} index(es): {', '.join(indexes)}"
    if stale:
        return result(
            "STAR index (RNA)", summary, "rebuild needed", "changed",
            f"built from an annotation that is no longer installed: "
            f"{'; '.join(stale)}. STAR would find junctions from one "
            f"annotation while the fusion caller names genes from another. "
            f"Rebuild with install_pipeline.py --only star-index --force "
            f"--rna-read-length <N>.",
            "")
    if unknown:
        return result(
            "STAR index (RNA)", summary, None, "manual",
            f"no build record in {', '.join(unknown)}, so the annotation "
            f"and read length they were built for cannot be checked. An "
            f"index built elsewhere carries none; confirm it matches the "
            f"installed GENCODE.", "")
    return result("STAR index (RNA)", summary, "consistent", "current",
                  "built from the installed annotation")


# =============================================================================
# DRIVER
# =============================================================================

def run_checks(data_dir):
    """Every check, in report order. Never raises."""
    checks = [
        lambda: check_vep(data_dir),
        lambda: check_pcgr(data_dir),
        lambda: check_msi_models(data_dir),
        check_snpeff,
        lambda: check_gatk_resources(data_dir),
        lambda: check_cosmic(data_dir),
        # The RNA branch. Both return None when it is not installed, so a
        # DNA-only machine's report is unchanged.
        lambda: check_gencode(data_dir),
        lambda: check_star_index(data_dir),
    ]
    out = []
    for fn in checks:
        try:
            answer = fn()
            if answer is None:
                continue                  # not installed; not applicable
            out.append(answer)
        except Exception as exc:              # noqa: BLE001 - reported inline
            out.append(result("(check failed)", "?", None, "unknown",
                              f"{type(exc).__name__}: {exc}"))
    return out


def build_status(data_dir):
    results = run_checks(data_dir)
    counts = {}
    for r in results:
        counts[r["state"]] = counts.get(r["state"], 0) + 1
    return {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "data_dir": data_dir,
        "counts": counts,
        "updates_available": counts.get("update", 0),
        "needs_attention": counts.get("update", 0) + counts.get("changed", 0),
        "results": results,
    }


def write_status(status, path):
    """Write the status file atomically, so a reader never sees half of it."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(status, fh, indent=2)
    os.replace(tmp, path)


def render(status):
    """Human-readable summary."""
    mark = {"current": "  ok  ", "update": "UPDATE", "manual": "manual",
            "unknown": "  ?   ", "changed": "CHANGED"}
    lines = [f"Database check — {status['checked_utc'][:19]}Z",
             f"data-dir: {status['data_dir']}", ""]
    for r in status["results"]:
        lines.append(f"[{mark.get(r['state'], '?'):^7}] {r['name']}")
        lines.append(f"          installed: {r['installed']}")
        if r.get("latest"):
            lines.append(f"          latest   : {r['latest']}")
        if r.get("note"):
            lines.append(f"          {r['note']}")
        if r.get("url") and r["state"] in ("update", "manual", "changed"):
            lines.append(f"          {r['url']}")
        lines.append("")
    n = status["needs_attention"]
    lines.append("Nothing to do." if not n else
                 f"{n} item(s) need attention. Nothing was changed — this "
                 f"tool only reports.")
    return "\n".join(lines)


def load_status(path=DEFAULT_STATUS):
    """
    Read a status file, for the web interface.

    Returns None when there is no readable file. Adds `stale` and `age_days`
    so a caller can distinguish "checked recently, all well" from "nobody
    has looked in three months".
    """
    try:
        with open(path, encoding="utf-8") as fh:
            status = json.load(fh)
    except (OSError, ValueError):
        return None
    try:
        when = datetime.fromisoformat(status["checked_utc"])
        age = (datetime.now(timezone.utc) - when).days
        status["age_days"] = age
        status["stale"] = age > STALE_AFTER_DAYS
    except (KeyError, ValueError):
        status["age_days"] = None
        status["stale"] = True
    return status


def main():
    parser = argparse.ArgumentParser(
        description="Check for reference-database updates and record the "
                    "result. Reports only: it never downloads, replaces or "
                    "modifies anything a pipeline run depends on.")
    parser.add_argument("--data-dir", default=DEFAULT_DATA,
                        help="Where the databases live (default: %(default)s).")
    parser.add_argument("--status-file", default=DEFAULT_STATUS,
                        help="Where to write the result. Deliberately outside "
                             "the data directory (default: %(default)s).")
    parser.add_argument("--print", dest="show", action="store_true",
                        help="Print a readable summary as well as writing it.")
    parser.add_argument("--json", action="store_true",
                        help="Print the status as JSON on stdout.")
    parser.add_argument("--no-write", action="store_true",
                        help="Do not write the status file.")
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 when something needs attention. Off by "
                             "default so a cron entry stays quiet.")
    args = parser.parse_args()

    status = build_status(args.data_dir)
    if not args.no_write:
        try:
            write_status(status, args.status_file)
        except OSError as exc:
            print(f"could not write {args.status_file}: {exc}",
                  file=sys.stderr)

    if args.json:
        print(json.dumps(status, indent=2))
    elif args.show or sys.stdout.isatty():
        print(render(status))

    return 1 if (args.strict and status["needs_attention"]) else 0


if __name__ == "__main__":
    sys.exit(main())
