#!/usr/bin/env python3
"""
jobs.py
=======
Background execution of pipeline runs for the web interface.

WHY A JOB RUNNER AT ALL
-----------------------
  A real run takes hours: BWA-MEM2 alone is 1-2 hours per sample, and
  building the reference index the first time is another 1-2. A Flask
  request cannot hold that open, so submitting a run starts a thread which
  spawns the pipeline as a subprocess and streams its output to a log file
  on disk. The browser then polls a status endpoint.

DESIGN CONSTRAINTS THIS ENCODES
-------------------------------
  * ONE RUN AT A TIME. The pipeline is already parallel internally (it is
    handed --threads) and a second concurrent run would contend for the
    same CPUs and RAM -- BWA-MEM2 wants ~32 GB for hg38 on its own. The
    queue is therefore serial, and submitting while busy queues rather
    than oversubscribing the machine.
  * THE LOG IS THE SOURCE OF TRUTH. Everything the pipeline prints goes to
    <run>/webapp_run.log verbatim. The UI's progress bar is derived from
    that text, never from a separate bookkeeping channel that could drift
    away from what actually happened.
  * PATIENT DETAILS NEVER ENTER THE LOG. They are written to patient.json
    beside the run and read back only when a report is generated. The
    pipeline subprocess is never told them, so they cannot end up in a
    tool's stdout, in a manifest, or in a crash trace.
  * NOTHING IS KILLED IMPLICITLY. Cancelling terminates the process group
    so GATK's Java children die too, rather than being orphaned to carry
    on consuming the machine.
"""

import glob
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import uuid
from collections import deque
from datetime import datetime, timezone

# Progress is inferred from the pipeline's own step banners, e.g.
# "# STEP 4/11: Base Quality Score Recalibration (BQSR)".
STEP_RE = re.compile(r"#\s*STEP\s+(\d+)\s*/\s*(\d+)\s*:\s*(.+?)\s*$")

# Tail of the log kept in memory for the live view, so the UI does not
# re-read a multi-megabyte file on every poll.
TAIL_LINES = 400


class Job:
    """One pipeline run: its parameters, its process, and its progress."""

    def __init__(self, job_id, run_dir, script, argv, patient, meta,
                 pcgr_form=None):
        self.id = job_id
        self.run_dir = run_dir
        self.script = script
        self.argv = argv
        self.patient = patient
        self.meta = meta
        # Settings for the PCGR follow-up, or None to skip it. Kept as the
        # raw form so build_pcgr_argv() stays the single place that knows
        # how the report is invoked.
        self.pcgr_form = pcgr_form
        self.pcgr_status = None

        self.status = "queued"      # queued|running|finished|failed|cancelled
        self.returncode = None
        self.step = 0
        self.total_steps = 0
        self.step_name = ""
        self.started = None
        self.finished = None
        self.error = None

        self.log_path = os.path.join(run_dir, "webapp_run.log")
        self._tail = deque(maxlen=TAIL_LINES)
        self._proc = None
        self._lock = threading.Lock()

    # -- state for the UI -------------------------------------------------
    def snapshot(self):
        with self._lock:
            pct = 0
            if self.total_steps:
                pct = int(100 * self.step / self.total_steps)
            elif self.status == "finished":
                pct = 100
            return {
                "id": self.id,
                "status": self.status,
                "returncode": self.returncode,
                "step": self.step,
                "total_steps": self.total_steps,
                "step_name": self.step_name,
                "percent": pct,
                "started": self.started,
                "finished": self.finished,
                "error": self.error,
                "run_dir": self.run_dir,
                "command": " ".join(self.argv),
                "sample": self.meta.get("tumour_sample"),
                "patient_id": self.patient.get("patient_id", ""),
                "output_dir": self.meta.get("output_dir"),
                "pcgr_status": self.pcgr_status,
            }

    def tail(self, limit=TAIL_LINES):
        with self._lock:
            return list(self._tail)[-limit:]

    # -- execution --------------------------------------------------------
    def run(self):
        with self._lock:
            self.status = "running"
            self.started = datetime.now(timezone.utc).isoformat()

        try:
            with open(self.log_path, "a", encoding="utf-8") as log:
                log.write(f"\n### webapp job {self.id}\n"
                          f"### {' '.join(self.argv)}\n")
                log.flush()
                # start_new_session puts the pipeline in its own process
                # group; cancel() then signals the whole group so GATK's
                # Java child processes go too.
                pipeline_env = find_conda_env(PIPELINE_ENV_NAME)
                self._proc = subprocess.Popen(
                    self.argv,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    cwd=os.path.dirname(self.script) or ".",
                    env=(activated_env(pipeline_env) if pipeline_env
                         else None),
                    start_new_session=True,
                )
                for line in self._proc.stdout:
                    line = line.rstrip("\n")
                    log.write(line + "\n")
                    log.flush()
                    self._observe(line)
                self._proc.wait()

            with self._lock:
                self.returncode = self._proc.returncode
                if self.status == "cancelled":
                    pass                       # keep the cancelled status
                elif self.returncode == 0:
                    self.status = "finished"
                    self.step = self.total_steps or self.step
                else:
                    self.status = "failed"
                    self.error = f"pipeline exited {self.returncode}"

            # The clinical report runs only after a clean pipeline, and its
            # failure never fails the run: the call set is already on disk
            # and is the thing that took the hours.
            if self.status == "finished" and self.pcgr_form:
                self._run_pcgr()
        except Exception as exc:                # noqa: BLE001 - reported to UI
            with self._lock:
                self.status = "failed"
                self.error = f"{type(exc).__name__}: {exc}"
        finally:
            with self._lock:
                self.finished = datetime.now(timezone.utc).isoformat()
            self._write_state()

    def _run_pcgr(self):
        """
        Generate the clinical report, in PCGR's own conda environment.

        Runs only after a clean pipeline. Every failure here is recorded and
        then swallowed: the variant calls are already on disk and represent
        hours of compute, so a missing report must not turn a finished run
        into a failed one. The outcome lands in self.pcgr_status so the UI
        can say which of the two happened instead of silently showing a run
        with no report.
        """
        output_dir = self.meta.get("output_dir")
        if not output_dir:
            self.pcgr_status = "skipped: no output directory recorded"
            return

        python_exe = find_pcgr_python(self.meta.get("pcgr_python"))
        if not python_exe:
            self.pcgr_status = ("skipped: no PCGR environment found "
                                "(set PCGR_PYTHON)")
            self._log_pcgr(
                "[PCGR] Skipped: could not find a Python interpreter with "
                "PCGR installed.\n"
                "[PCGR] PCGR lives in its own conda environment; point at it "
                "with the PCGR_PYTHON environment variable, e.g.\n"
                "[PCGR]   export PCGR_PYTHON=~/miniconda3/envs/pcgr/bin/python\n")
            return

        manifest = latest_manifest(output_dir)
        if manifest is None:
            self.pcgr_status = "skipped: no run manifest found"
            self._log_pcgr("[PCGR] Skipped: no pipeline manifest in "
                           f"{output_dir}.\n")
            return

        argv, note = build_pcgr_argv(
            python_exe, os.path.dirname(self.script) or ".",
            manifest, self.pcgr_form, output_dir)
        if argv is None:
            self.pcgr_status = f"skipped: {note}"
            self._log_pcgr(f"[PCGR] Skipped: {note}.\n")
            return

        with self._lock:
            self.step_name = "Clinical report (PCGR)"

        self._log_pcgr(f"\n### clinical report\n### {' '.join(argv)}\n")
        try:
            with open(self.log_path, "a", encoding="utf-8") as log:
                pcgr_env = find_conda_env(PCGR_ENV_NAME)
                proc = subprocess.Popen(
                    argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1,
                    cwd=os.path.dirname(self.script) or ".",
                    # CONDA_PREFIX must point at the PCGR env: PCGR reads it
                    # to find its VEP plugin directory, and a PATH-only fix
                    # fails later with "No ensembl-vep directories found".
                    env=activated_env(pcgr_env) if pcgr_env else None,
                    start_new_session=True,
                )
                with self._lock:
                    self._proc = proc          # so cancel() reaches it too
                for line in proc.stdout:
                    line = line.rstrip("\n")
                    log.write(line + "\n")
                    log.flush()
                    with self._lock:
                        self._tail.append(line)
                proc.wait()
            if proc.returncode == 0:
                reports = glob.glob(os.path.join(output_dir, "pcgr", "*.html"))
                self.pcgr_status = ("report generated" if reports
                                    else "ran, but produced no HTML report")
            else:
                self.pcgr_status = f"failed (exit {proc.returncode})"
        except Exception as exc:               # noqa: BLE001 - reported to UI
            self.pcgr_status = f"failed ({type(exc).__name__}: {exc})"
        self._write_state()

    def _log_pcgr(self, text):
        """Append a note to the run log, ignoring log-write failures."""
        try:
            with open(self.log_path, "a", encoding="utf-8") as log:
                log.write(text)
        except OSError:
            pass

    def _observe(self, line):
        """Update progress from one line of pipeline output."""
        with self._lock:
            self._tail.append(line)
            match = STEP_RE.search(line)
            if match:
                self.step = int(match.group(1))
                self.total_steps = int(match.group(2))
                self.step_name = match.group(3)

    def cancel(self):
        """Terminate the run and everything it spawned."""
        with self._lock:
            if self.status not in ("queued", "running"):
                return False
            self.status = "cancelled"
            self.error = "cancelled by user"
            proc = self._proc
        if proc and proc.poll() is None:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except (ProcessLookupError, PermissionError, OSError):
                proc.terminate()
        return True

    def _write_state(self):
        """Persist the job record so it survives a webapp restart."""
        try:
            with open(os.path.join(self.run_dir, "webapp_job.json"), "w",
                      encoding="utf-8") as fh:
                json.dump({"job": self.snapshot(), "meta": self.meta},
                          fh, indent=2, default=str)
        except OSError:
            pass


class JobManager:
    """A serial queue of pipeline runs."""

    def __init__(self, runs_root):
        self.runs_root = os.path.abspath(runs_root)
        os.makedirs(self.runs_root, exist_ok=True)
        self._jobs = {}
        self._order = []
        self._lock = threading.Lock()
        self._worker = None
        self._queue = deque()

    def submit(self, script, argv, patient, meta, pcgr_form=None):
        """Create a run directory, record the patient, and queue the job."""
        job_id = datetime.now().strftime("%Y%m%dT%H%M%S") + "-" + \
            uuid.uuid4().hex[:6]
        run_dir = os.path.join(self.runs_root, job_id)
        os.makedirs(run_dir, exist_ok=True)

        # Patient details live here and ONLY here -- never on the pipeline's
        # command line, and never in the log.
        with open(os.path.join(run_dir, "patient.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(patient, fh, indent=2)

        job = Job(job_id, run_dir, script, argv, patient, meta,
                  pcgr_form=pcgr_form)
        with self._lock:
            self._jobs[job_id] = job
            self._order.append(job_id)
            self._queue.append(job)
            if self._worker is None or not self._worker.is_alive():
                self._worker = threading.Thread(target=self._drain,
                                                daemon=True)
                self._worker.start()
        return job

    def _drain(self):
        """Run queued jobs one at a time."""
        while True:
            with self._lock:
                if not self._queue:
                    self._worker = None
                    return
                job = self._queue.popleft()
            job.run()

    def get(self, job_id):
        with self._lock:
            return self._jobs.get(job_id)

    def all(self):
        with self._lock:
            return [self._jobs[i] for i in reversed(self._order)]

    def busy(self):
        with self._lock:
            return any(j.status in ("running", "queued")
                       for j in self._jobs.values())

    def load_existing(self):
        """
        Re-attach to runs from previous webapp sessions.

        Their processes are gone, so they are listed with whatever status
        was last persisted; this exists so a restart does not lose the
        history (and therefore the ability to generate a report).
        """
        for entry in sorted(os.listdir(self.runs_root)):
            run_dir = os.path.join(self.runs_root, entry)
            state = os.path.join(run_dir, "webapp_job.json")
            if entry in self._jobs or not os.path.isfile(state):
                continue
            try:
                with open(state, encoding="utf-8") as fh:
                    saved = json.load(fh)
                patient = {}
                pj = os.path.join(run_dir, "patient.json")
                if os.path.isfile(pj):
                    with open(pj, encoding="utf-8") as fh:
                        patient = json.load(fh)
            except (OSError, json.JSONDecodeError):
                continue
            snap = saved.get("job", {})
            job = Job(entry, run_dir, "", [], patient, saved.get("meta", {}))
            job.status = snap.get("status", "finished")
            job.returncode = snap.get("returncode")
            job.step = snap.get("step", 0)
            job.total_steps = snap.get("total_steps", 0)
            job.started = snap.get("started")
            job.finished = snap.get("finished")
            job.error = snap.get("error")
            # An interrupted webapp leaves jobs stuck as "running"; they are
            # not, because the process died with the server.
            if job.status in ("running", "queued"):
                job.status = "interrupted"
                job.error = "webapp restarted while this run was in progress"
            with self._lock:
                self._jobs[entry] = job
                self._order.append(entry)


# =============================================================================
# INTERMEDIATE FILE CLEANUP
# =============================================================================
# A run leaves about 5 GB of intermediates against 70 MB worth keeping: the
# cleaned FASTQs, the aligned, duplicate-marked and recalibrated BAMs. Each
# is reproducible from the inputs, and none is needed once the VCFs and the
# report exist. Nothing is deleted automatically -- reprocessing costs an
# hour, and that is the user's call to make, not this module's.

# Reproducible from the FASTQs; safe to remove after a finished run.
INTERMEDIATE_DIRS = ("cleaned_fastq", "aligned", "dedup", "bqsr", "tmp")

# Never touched: the call sets, the report, and everything needed to say how
# they were produced.
PRESERVED_DIRS = ("mutect2", "annotated", "pcgr", "metrics", "stats",
                  "logs", "fastp_reports", "reference")


def _tree_size(path):
    """Bytes used by a directory tree, ignoring anything unreadable."""
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total


def survey_intermediates(output_dir):
    """
    Report what cleanup would remove, without removing anything.

    Returns (entries, total_bytes) where entries is a list of
    (name, bytes) for the intermediate directories that exist. The UI
    shows this before asking for confirmation, so the decision is made
    against real numbers rather than a vague promise of space.
    """
    entries = []
    if not output_dir or not os.path.isdir(output_dir):
        return entries, 0
    for name in INTERMEDIATE_DIRS:
        path = os.path.join(output_dir, name)
        if os.path.isdir(path):
            entries.append((name, _tree_size(path)))
    return entries, sum(size for _n, size in entries)


def cleanup_intermediates(output_dir, keep_bam=False):
    """
    Delete the reproducible intermediates from a finished run.

    Only the directories this pipeline creates are touched, by name, and
    only inside output_dir -- a symlinked entry is unlinked rather than
    followed, so a link pointing at shared data cannot take that data with
    it. The VCFs, the report, the metrics and the logs are never touched:
    what remains is enough to read the result and to say how it was made.

    ARGS:
        output_dir: The run's --output-dir.
        keep_bam:   Keep bqsr/, the analysis-ready BAM. Choose this when
                    the calls may be revisited: re-calling from the BAM
                    takes minutes, while rebuilding it takes about an hour.

    RETURNS:
        (freed_bytes, removed_names, errors)
    """
    freed, removed, errors = 0, [], []
    if not output_dir or not os.path.isdir(output_dir):
        return 0, [], ["output directory not found: %s" % output_dir]

    root = os.path.realpath(output_dir)
    for name in INTERMEDIATE_DIRS:
        if keep_bam and name == "bqsr":
            continue
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        # Refuse anything that is not a real directory inside this run.
        if os.path.islink(path):
            errors.append("%s is a symlink; left alone" % name)
            continue
        if os.path.commonpath([root, os.path.realpath(path)]) != root:
            errors.append("%s resolves outside the run; left alone" % name)
            continue
        size = _tree_size(path)
        try:
            shutil.rmtree(path)
            freed += size
            removed.append(name)
        except OSError as exc:
            errors.append("%s: %s" % (name, exc))
    return freed, removed, errors


def human_bytes(n):
    """Size as a short human string (1.8 GB, 240 MB)."""
    step = 1000.0
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if abs(n) < step or unit == "TB":
            return "%.0f %s" % (n, unit) if unit == "B" else \
                   "%.1f %s" % (n, unit)
        n /= step
    return "%.1f TB" % n


# =============================================================================
# CONDA ENVIRONMENTS
# =============================================================================
# The pipeline and PCGR live in DIFFERENT conda environments and neither can
# run in the other. Launching a subprocess with "<env>/bin/python" is not
# enough on its own: the scripts locate their tools with shutil.which(), and
# PCGR additionally resolves its VEP plugin directory from $CONDA_PREFIX. So
# each child gets a properly activated environment -- its bin/ on PATH and
# CONDA_PREFIX pointed at it -- rather than inheriting whatever shell started
# the web server. Without this the webapp works only when launched from an
# already-activated env, and fails with "Required tools not found" otherwise.

PIPELINE_ENV_NAME = "cancer_pipeline"
PCGR_ENV_NAME = "pcgr"


def find_conda_env(name, explicit=None):
    """
    Locate a conda environment directory by name.

    Order: an explicit path, then the sibling of the env this process runs
    from (the usual case -- envs live together under one envs/ directory),
    then $CONDA_ROOT and the common miniconda/anaconda locations under the
    home directory. Returns None when nothing matches, which callers report
    rather than guessing.
    """
    candidates = []
    if explicit:
        candidates.append(explicit)

    here = os.path.dirname(os.path.dirname(os.path.abspath(sys.executable)))
    envs = os.path.dirname(here)
    if os.path.basename(envs) == "envs":
        candidates.append(os.path.join(envs, name))
    # Running from the base env: envs/ sits inside it.
    candidates.append(os.path.join(here, "envs", name))

    for root in (os.environ.get("CONDA_ROOT"),
                 os.path.expanduser("~/miniconda3"),
                 os.path.expanduser("~/anaconda3"),
                 os.path.expanduser("~/miniforge3")):
        if root:
            candidates.append(os.path.join(root, "envs", name))

    for path in candidates:
        if path and os.path.isdir(os.path.join(path, "bin")):
            return path
    return None


def activated_env(env_dir):
    """
    An os.environ copy equivalent to "conda activate <env_dir>".

    Prepends the environment's bin/ to PATH so shutil.which() finds its
    tools, and sets CONDA_PREFIX, which PCGR reads to locate its VEP plugin
    directory. PYTHONHOME and PYTHONPATH are removed: inherited from another
    interpreter they make the child import the wrong standard library.
    """
    env = dict(os.environ)
    bin_dir = os.path.join(env_dir, "bin")
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = env_dir
    env["CONDA_DEFAULT_ENV"] = os.path.basename(env_dir)
    env.pop("PYTHONHOME", None)
    env.pop("PYTHONPATH", None)
    return env


def env_python(env_dir):
    """Path to an environment's interpreter, or None."""
    if not env_dir:
        return None
    exe = os.path.join(env_dir, "bin", "python")
    return exe if os.path.isfile(exe) and os.access(exe, os.X_OK) else None


# =============================================================================
# PCGR FOLLOW-UP
# =============================================================================
# PCGR cannot run inside the pipeline's own conda environment: it resolves its
# VEP plugin directory from $CONDA_PREFIX/share, so a PATH tweak is not enough
# and step 11 always reports itself as skipped. The webapp therefore runs it as
# a SECOND process, in PCGR's environment, once the pipeline has succeeded.
# Without this the web interface can never produce a clinical report, which is
# the single most visible thing it is supposed to deliver.

def find_pcgr_python(explicit=None):
    """
    Interpreter for the PCGR environment, or None.

    $PCGR_PYTHON wins when set, so an unusual layout can be pointed at
    directly; otherwise the environment is found by name.
    """
    for path in (explicit, os.environ.get("PCGR_PYTHON")):
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return env_python(find_conda_env(PCGR_ENV_NAME))


def latest_manifest(output_dir):
    """The newest pipeline manifest in output_dir, parsed, or None."""
    paths = sorted(glob.glob(os.path.join(output_dir,
                                          "pipeline_manifest_*.json")))
    if not paths:
        return None
    try:
        with open(paths[-1], encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def build_pcgr_argv(python_exe, script_dir, manifest, form, output_dir):
    """
    Build the pcgr_report.py command for a finished run.

    Returns (argv, note). argv is None when the report cannot be produced,
    and note then explains why in one line for the run log.
    """
    sample = manifest.get("tumour") or \
        (manifest.get("result") or {}).get("tumour")
    if not sample:
        return None, "no tumour sample recorded in the manifest"

    # COSMIC-annotated if the step ran, otherwise the filtered calls. Never
    # the SnpEff output: PCGR runs its own VEP and validates INFO fields.
    vcf = manifest.get("cosmic_vcf") or manifest.get("filtered_vcf")
    if not vcf or not os.path.exists(vcf):
        return None, f"no filtered VCF on disk to report on ({vcf})"

    argv = [
        python_exe, os.path.join(script_dir, "pcgr_report.py"),
        "--input-vcf", vcf,
        "--output-dir", os.path.join(output_dir, "pcgr"),
        "--sample-id", sample,
        "--pcgr-refdata-dir", form["pcgr_refdata_dir"],
        "--lift-tags",
    ]
    if form.get("vep_dir"):
        argv += ["--vep-dir", form["vep_dir"]]
    if form.get("pcgr_assay"):
        argv += ["--assay", form["pcgr_assay"]]
    if form.get("pcgr_tumour_site"):
        argv += ["--tumour-site", str(form["pcgr_tumour_site"])]
    if form.get("pcgr_target_size_mb"):
        argv += ["--effective-target-size-mb",
                 str(form["pcgr_target_size_mb"])]
    for flag, key in (("--estimate-tmb", "pcgr_estimate_tmb"),
                      ("--estimate-msi", "pcgr_estimate_msi"),
                      ("--estimate-signatures", "pcgr_estimate_signatures")):
        if form.get(key):
            argv.append(flag)
    # Tumour-only unless the run had a matched normal.
    if not (manifest.get("normal") or
            (manifest.get("result") or {}).get("normal")):
        argv.append("--tumour-only")
    return argv, None


def build_pipeline_argv(python_exe, script, form):
    """
    Turn the submitted form into a comprehensive_variant_calling.py argv.

    Only options the form exposes are emitted; anything absent is left to
    the script's own defaults, so this wrapper never silently changes the
    pipeline's behaviour.
    """
    argv = [python_exe, script,
            "--output-dir", form["output_dir"],
            "--reference", form["reference"]]

    if form.get("manifest"):
        argv += ["--manifest", form["manifest"]]
    else:
        argv += ["--input-dir", form["input_dir"]]
        if form.get("auto_discover"):
            argv += ["--auto-discover"]
        else:
            argv += ["--tumour-sample", form["tumour_sample"],
                     "--tumour-r1", form["tumour_r1"],
                     "--tumour-r2", form["tumour_r2"]]
            if form.get("normal_sample"):
                argv += ["--normal-sample", form["normal_sample"],
                         "--normal-r1", form["normal_r1"],
                         "--normal-r2", form["normal_r2"]]

    for flag, key in (("--cosmic", "cosmic"),
                      ("--dbsnp", "dbsnp"),
                      ("--germline-resource", "germline_resource"),
                      ("--panel-of-normals", "panel_of_normals"),
                      ("--contamination-resource", "contamination_resource"),
                      ("--intervals", "intervals")):
        if form.get(key):
            argv += [flag, form[key]]

    # --known-indels takes a LIST, so it cannot join the loop above. It was
    # missing entirely until 2026-09-08, which meant BQSR ran on dbSNP alone
    # however the form was filled in -- indel recalibration silently absent.
    if form.get("known_indels"):
        paths = [p for p in str(form["known_indels"]).split() if p]
        if paths:
            argv += ["--known-indels"] + paths

    # Numeric options: 0 / blank means "off", so they are only passed when
    # actually set. --interval-padding is meaningless without --intervals.
    if form.get("intervals") and form.get("interval_padding"):
        argv += ["--interval-padding", str(form["interval_padding"])]
    if form.get("min_depth"):
        argv += ["--min-depth", str(form["min_depth"])]
    if form.get("min_allele_fraction"):
        argv += ["--min-allele-fraction", str(form["min_allele_fraction"])]
    if form.get("pcgr_target_size_mb"):
        argv += ["--pcgr-target-size-mb", str(form["pcgr_target_size_mb"])]

    if form.get("threads"):
        argv += ["--threads", str(form["threads"])]
    if form.get("two_pass"):
        argv += ["--two-pass"]
    if form.get("resume"):
        argv += ["--resume"]
    if form.get("dry_run"):
        argv += ["--dry-run"]

    if form.get("pcgr_refdata_dir"):
        argv += ["--pcgr-refdata-dir", form["pcgr_refdata_dir"]]
        if form.get("vep_dir"):
            argv += ["--vep-dir", form["vep_dir"]]
        if form.get("pcgr_assay"):
            argv += ["--pcgr-assay", form["pcgr_assay"]]
        if form.get("pcgr_lift_tags"):
            argv += ["--pcgr-lift-tags"]

    if form.get("skip_steps"):
        argv += ["--skip-steps"] + form["skip_steps"].split()

    return argv
