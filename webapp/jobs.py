#!/usr/bin/env python3
# Created by Brainstorm, 2026.
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
  * A RUN IS NOT FINISHED UNTIL EVERYTHING IT COMPRISES IS. The clinical
    report runs after the pipeline exits, in PCGR's own conda environment,
    and takes up to an hour more. That phase has its own non-terminal
    status, "reporting": the page keeps polling, keeps offering Cancel and
    withholds the report buttons until it ends. Calling the run finished at
    the pipeline's exit -- which is what this used to do -- meant the UI
    announced a completed run, and a PDF generated in that window recorded
    "PCGR did not run" about a PCGR that was running as it was written.

STATUSES
--------
  queued -> running -> [reporting ->] finished
                    +-> failed | cancelled
  interrupted is assigned on startup to a run whose webapp died under it;
  the process went with the server, so whatever it was doing is over.
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
                 pcgr_form=None, assay="dna"):
        self.id = job_id
        self.run_dir = run_dir
        self.script = script
        self.argv = argv
        self.patient = patient
        self.meta = meta
        # "dna" or "rna". Decides which conda environment the job runs in,
        # whether the PCGR follow-up applies at all (it does not on RNA --
        # PCGR is a somatic SNV reporter and has nothing to say about a
        # fusion), and which report the UI offers.
        self.assay = assay if assay in ("dna", "rna") else "dna"
        # Settings for the PCGR follow-up, or None to skip it. Kept as the
        # raw form so build_pcgr_argv() stays the single place that knows
        # how the report is invoked.
        self.pcgr_form = pcgr_form
        self.pcgr_status = None

        # queued|running|reporting|finished|failed|cancelled|interrupted.
        # "reporting" is the PCGR follow-up: the pipeline has exited but the
        # run is NOT done, and calling it finished there is what made the
        # page offer a PDF that said "PCGR did not run" while PCGR was in
        # fact running. Only the last thing this job does may say finished.
        self.status = "queued"
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
                "assay": self.assay,
                # The report phase has no step number of its own: the
                # pipeline's count is what the log shows, and the UI shows
                # the bar working rather than a step that does not exist.
                "indeterminate": self.status == "reporting",
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
        # Written up front so a webapp killed mid-run still has a record to
        # re-attach to; load_existing() skips a run directory with no state
        # file, so without this an interrupted run vanished from the list
        # rather than showing as interrupted.
        self._write_state()

        try:
            with open(self.log_path, "a", encoding="utf-8") as log:
                log.write(f"\n### webapp job {self.id}\n"
                          f"### {' '.join(self.argv)}\n")
                log.flush()
                # start_new_session puts the pipeline in its own process
                # group; cancel() then signals the whole group so GATK's
                # Java child processes go too.
                # DNA jobs run in cancer_pipeline, RNA jobs in
                # cancer_rna. Picked from the job rather than from the
                # script path so an operator running a custom script still
                # gets the environment their assay needs.
                env_name = (RNA_ENV_NAME if self.assay == "rna"
                            else PIPELINE_ENV_NAME)
                pipeline_env = find_conda_env(env_name)
                if pipeline_env is None:
                    log.write(f"### conda environment '{env_name}' not "
                              f"found; running with the current "
                              f"interpreter's environment\n")
                    log.flush()
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
                # `with` on the pipe, not just the loop. Reading to EOF
                # does not close the read end -- CPython releases it only
                # when the Popen object is collected, and this manager
                # keeps every Job for the life of the server so that
                # finished runs stay listable. Without this, each
                # completed run held a file descriptor until the process
                # exited, and a long-lived server working through a
                # worksheet ran out of them.
                with self._proc.stdout as stream:
                    for line in stream:
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
                    # A clean pipeline with a report still to come is
                    # "reporting", not "finished" -- the page keeps polling,
                    # keeps offering Cancel, and withholds the report
                    # buttons until there is something behind them.
                    if self.pcgr_form:
                        self.status = "reporting"
                        self.step = self.total_steps or self.step
                        self.step_name = "Clinical report (PCGR)"
                    else:
                        self.status = "finished"
                        self.step = self.total_steps or self.step
                else:
                    self.status = "failed"
                    self.error = f"pipeline exited {self.returncode}"

            # Persist the handover: a webapp restarted during the report
            # would otherwise find no state file at all for this run.
            self._write_state()

            # The clinical report runs only after a clean pipeline, and its
            # failure never fails the run: the call set is already on disk
            # and is the thing that took the hours.
            if self.status == "reporting":
                self._run_pcgr()
                with self._lock:
                    # Cancelling during the report is a cancelled run, not a
                    # finished one; anything else is now genuinely done.
                    if self.status == "reporting":
                        self.status = "finished"
                        self.step = self.total_steps or self.step
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

        manifest = latest_manifest(output_dir, self.assay)
        if manifest is None:
            self.pcgr_status = "skipped: no run manifest found"
            self._log_pcgr(f"[PCGR] Skipped: no {self.assay.upper()} run "
                           f"manifest in {output_dir}.\n")
            return

        # The two branches hand PCGR different molecular inputs -- a
        # somatic VCF, or a fusion table -- so they build different
        # commands. Everything after this point is shared: same
        # environment, same streaming, same status reporting.
        builder = (build_rna_pcgr_argv if self.assay == "rna"
                   else build_pcgr_argv)
        argv, note = builder(
            python_exe, os.path.dirname(self.script) or ".",
            manifest, self.pcgr_form, output_dir)
        if argv is None:
            self.pcgr_status = f"skipped: {note}"
            self._log_pcgr(f"[PCGR] Skipped: {note}.\n")
            return

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
                # Same reason as the pipeline loop above: close the pipe
                # here rather than leaving it to the collector.
                with proc.stdout as stream:
                    for line in stream:
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
            if self.status not in ("queued", "running", "reporting"):
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

    def submit(self, script, argv, patient, meta, pcgr_form=None,
               assay="dna"):
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
                  pcgr_form=pcgr_form, assay=assay)
        # Persisted NOW, while it is still only queued. load_existing()
        # rebuilds the run list from these files and skips a directory
        # without one, so a job that had not started when the server
        # restarted used to vanish from the list without trace -- after an
        # --update, which install.md says to follow with a restart, every
        # queued worksheet row but the running one disappeared. With the
        # record written up front it comes back marked "interrupted", which
        # is the truth, and can be resubmitted.
        job._write_state()
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
            return any(j.status in ("running", "queued", "reporting")
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
            job = Job(entry, run_dir, "", [], patient, saved.get("meta", {}),
                      assay=snap.get("assay")
                      or saved.get("meta", {}).get("assay") or "dna")
            job.status = snap.get("status", "finished")
            job.returncode = snap.get("returncode")
            job.step = snap.get("step", 0)
            job.total_steps = snap.get("total_steps", 0)
            job.started = snap.get("started")
            job.finished = snap.get("finished")
            job.error = snap.get("error")
            job.pcgr_status = snap.get("pcgr_status")
            job.step_name = snap.get("step_name", "")
            # An interrupted webapp leaves jobs stuck as "running"; they are
            # not, because the process died with the server. "reporting"
            # goes the same way -- PCGR died with it too, and the run has no
            # report to show for the status it was left in.
            if job.status in ("running", "queued", "reporting"):
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
PRESERVED_DIRS = ("mutect2", "annotated", "pcgr", "coverage", "metrics",
                  "stats", "logs", "fastp_reports", "reference")


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
# The RNA branch runs in its own environment: STAR and Arriba are not in
# cancer_pipeline, and a job launched in the wrong one fails at the first
# command with "STAR: not found" after the queue has already claimed the
# machine. Job.assay decides which is activated.
RNA_ENV_NAME = "cancer_rna"
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


def manifest_vcf(manifest):
    """
    The call set a DNA manifest points at -- the first one that EXISTS.

    Order of preference is COSMIC-annotated, then filtered: the annotated
    file carries the known-mutation identifiers and is the better input.
    But preference is not availability. Naming the preferred file and then
    testing only that one means a run whose COSMIC step was skipped -- no
    COSMIC VCF configured, or the step failed -- reports "no VCF on disk"
    while its filtered calls sit right there, which is a report lost to a
    file that was never going to be written.
    """
    for key in ("cosmic_vcf", "filtered_vcf"):
        path = (manifest or {}).get(key)
        if path and os.path.exists(path):
            return path
    return None


def latest_manifest(output_dir, assay="dna"):
    """
    The newest run manifest in output_dir, parsed, or None.

    The two branches write differently named manifests -- the DNA engine
    writes pipeline_manifest_*.json, the RNA engine rna_manifest_*.json --
    so the pattern follows the assay. Before this was assay-aware, an RNA
    run's report step found no pipeline manifest and skipped itself with
    "no run manifest found", which was true and completely misleading.
    """
    pattern = ("rna_manifest_*.json" if assay == "rna"
               else "pipeline_manifest_*.json")
    paths = sorted(glob.glob(os.path.join(output_dir, pattern)))
    if not paths:
        return None
    try:
        with open(paths[-1], encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def find_coverage_reports(output_dir):
    """
    The per-sample coverage reports for a run, newest naming first.

    Globbed rather than read from the manifest: they are written before the
    manifest is finalised, and the point of them is to be readable WHILE
    the run is still going.
    """
    if not output_dir:
        return []
    pattern = os.path.join(output_dir, "coverage", "*.coverage.html")
    return sorted(glob.glob(pattern))


def find_rna_reports(output_dir, kind):
    """
    The per-sample RNA reports for a run: 'fusion' or 'qc'.

    Globbed for the same reason the coverage reports are -- they are
    written before the run manifest is finalised, and the point of the
    library QC report in particular is to be readable while the run is
    still going. An operator who can see at step 5 that the library failed
    does not have to wait for step 6 to tell them the empty table means
    nothing.
    """
    if not output_dir:
        return []
    pattern = {
        "fusion": os.path.join(output_dir, "fusion_report", "*.fusions.html"),
        "qc": os.path.join(output_dir, "rna_qc", "*.rna_qc.html"),
    }.get(kind)
    return sorted(glob.glob(pattern)) if pattern else []


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
    vcf = manifest_vcf(manifest)
    if not vcf:
        return None, ("no variant calls on disk to report on (looked for "
                      "the COSMIC-annotated and filtered VCFs named in the "
                      "run manifest)")

    # PCGR rejects an input VCF with no VEP cache. The form refuses this
    # combination at submit time; this catches a job restored from an older
    # record, or one submitted straight to the API, and turns an hour-late
    # failure into a skip that names its cause.
    if not form.get("vep_dir"):
        return None, ("no VEP cache configured, and PCGR will not annotate a "
                      "VCF without one")

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
    # The TMB denominator, and the one place it was still going wrong.
    #
    # This app produces the clinical report itself, in PCGR's own conda
    # environment, AFTER the pipeline has exited -- so it does not inherit
    # the footprint the pipeline measured from the target BED. Without
    # this, a panel run through the web form reached PCGR with no target
    # size at all and PCGR divided by its assumed 34 Mb: on a 2 Mb panel,
    # a TMB roughly 17x too low, printed in a report as a number with no
    # indication that anything was assumed.
    #
    # The form's own value wins when somebody typed one -- that is a
    # considered figure, often callable bases from a validation, and a raw
    # BED footprint is a cruder measure. Otherwise the pipeline's
    # measurement is used, and only then PCGR's default.
    #
    # Note what is NOT passed: --panel. pcgr_report.py would apply the
    # profile's estimate toggles for any flag absent from this command
    # line, and an unticked checkbox is absent -- so a profile could turn
    # back on an estimate the operator had just turned off.
    target_size = form.get("pcgr_target_size_mb")
    measured = (manifest.get("panel_profile") or {}).get(
        "target_size_mb_measured")
    if not target_size and measured:
        target_size = measured
    if target_size:
        argv += ["--effective-target-size-mb", str(target_size)]
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


def build_rna_argv(python_exe, script, form):
    """
    Turn the submitted RNA form into a fusion_calling.py argv.

    Separate from build_pipeline_argv() rather than a branch inside it,
    for the same reason the engines are separate: almost nothing is shared.
    The DNA builder's vocabulary -- intervals, padding, allele fractions,
    known sites, panel of normals, PCGR -- has no meaning on this branch,
    and a single builder threading two option sets through one function is
    how an RNA run ends up with a --min-allele-fraction nobody reads.

    As on the DNA side, only options the form exposes are emitted; anything
    absent is left to the script's own defaults, so this wrapper never
    quietly changes the engine's behaviour.
    """
    argv = [python_exe, script,
            "--output-dir", form["output_dir"],
            "--reference", form["reference"],
            "--gtf", form["gtf"]]

    # The panel is passed through as well as being applied to the form, for
    # the same two reasons as on the DNA side: the run manifest then records
    # the assay as data, and the profile carries settings the form has no
    # field for (UMI layout, adapters, caller passthroughs).
    if form.get("panel"):
        argv += ["--panel", form["panel"]]
    for extra in str(form.get("panel_file", "")).split():
        argv += ["--panel-file", extra]

    if form.get("star_index"):
        argv += ["--star-index", form["star_index"]]
    if form.get("reference_dir"):
        argv += ["--reference-dir", form["reference_dir"]]
    if form.get("read_length"):
        argv += ["--read-length", str(form["read_length"])]
    if form.get("arriba_resources"):
        argv += ["--arriba-resources", form["arriba_resources"]]

    if form.get("manifest"):
        argv += ["--manifest", form["manifest"]]
        if form.get("sample"):
            argv += ["--sample", form["sample"]]
    else:
        argv += ["--input-dir", form["input_dir"]]
        if form.get("auto_discover"):
            argv += ["--auto-discover"]
            # Narrows discovery to one library -- a worksheet row, or a
            # lane-split sample run by name so every lane is merged.
            if form.get("sample"):
                argv += ["--sample", form["sample"]]
        else:
            argv += ["--sample", form["sample"], "--r1", form["r1"]]
            # Single-end is legitimate here (Ion Torrent), so an absent R2
            # is a library type rather than a missing field.
            if form.get("r2"):
                argv += ["--r2", form["r2"]]

    if form.get("fusion_min_confidence"):
        argv += ["--fusion-min-confidence", form["fusion_min_confidence"]]
    if form.get("min_fusion_reads"):
        argv += ["--min-fusion-reads", str(form["min_fusion_reads"])]
    if form.get("strandedness"):
        argv += ["--strandedness", form["strandedness"]]
    if form.get("rna_min_reads_millions"):
        argv += ["--rna-min-reads-millions",
                 str(form["rna_min_reads_millions"])]
    if form.get("rna_min_unique_mapped_pct"):
        argv += ["--rna-min-unique-mapped-pct",
                 str(form["rna_min_unique_mapped_pct"])]
    if form.get("min_read_length"):
        argv += ["--min-read-length", str(form["min_read_length"])]
    if form.get("star_extra_args"):
        argv += ["--star-extra-args", str(form["star_extra_args"])]
    if form.get("arriba_extra_args"):
        argv += ["--arriba-extra-args", str(form["arriba_extra_args"])]

    if form.get("threads"):
        argv += ["--threads", str(form["threads"])]
    if form.get("dry_run"):
        argv += ["--dry-run"]
    if form.get("skip_steps"):
        argv += ["--skip-steps"] + str(form["skip_steps"]).split()

    # No PCGR options, and not for the DNA branch's reason. PCGR is a
    # somatic SNV/indel reporter: it takes a VCF and has nothing whatever
    # to say about a rearrangement. The RNA branch's interpretation layer
    # is fusion_report.py, which the engine runs itself as step 6.
    return argv


def build_rna_pcgr_argv(python_exe, script_dir, manifest, form, output_dir):
    """
    Build the pcgr_report.py command for a finished RNA run.

    PCGR 2.x accepts RNA fusions as a molecular input in their own right,
    and requires a VCF only if you give it one -- so an RNA run gets a real
    clinical interpretation (actionable fusions placed against the tumour
    site) rather than only the pipeline's own ranked table.

    THE COMBINED REPORT is the case worth having. When the RNA run was
    linked to a DNA run on the form, the DNA library's filtered VCF is
    passed alongside the fusions, and PCGR produces ONE report covering
    both -- which is how a specimen is reported clinically, rather than as
    two documents somebody has to reconcile by eye.

    Returns (argv, note); argv is None when no report can be produced, and
    note then says why in one line for the run log.
    """
    samples = manifest.get("samples") or {}
    if not samples:
        return None, "no samples recorded in the RNA manifest"

    # One PCGR report per run, for THE sample the run is about. A batch of
    # RNA samples is a batch of specimens, and merging their fusions into
    # one report would attribute one specimen's finding to another -- and
    # so would picking "the first" of several, which is what this did.
    # The form now refuses a run that would discover several unnamed
    # samples, so a named sample or a single one is all that reaches here.
    if form.get("sample") in samples:
        sample = form["sample"]
    elif len(samples) == 1:
        sample = next(iter(samples))
    else:
        return None, (f"{len(samples)} samples in one run and none named, "
                      f"so no report can be attributed; use the worksheet "
                      f"for several specimens")
    record = samples[sample]
    fusion_tsv = (record.get("fusion_report") or {}).get("pcgr_tsv")
    if not fusion_tsv or not os.path.exists(fusion_tsv):
        return None, ("no fusions cleared the confidence bar, so there is "
                      "nothing for PCGR to interpret (it rejects an empty "
                      "fusion file)")

    if not form.get("pcgr_refdata_dir"):
        return None, "no PCGR reference bundle configured for this run"

    argv = [
        python_exe, os.path.join(script_dir, "pcgr_report.py"),
        "--input-rna-fusion", fusion_tsv,
        "--output-dir", os.path.join(output_dir, "pcgr"),
        "--sample-id", sample,
        "--pcgr-refdata-dir", form["pcgr_refdata_dir"],
    ]

    # The DNA half of the same specimen, when the operator linked it.
    #
    # RESOLVED HERE, NOT AT SUBMIT TIME. On a hybrid panel both runs are
    # queued together and the DNA run has not produced a VCF yet when the
    # RNA run is created -- so a pairing captured at submit time would
    # always be empty, and the combined report would silently become two
    # separate ones. What is stored is the DNA run's output DIRECTORY; its
    # manifest is read now, when the DNA run has actually finished.
    paired_vcf = form.get("paired_vcf")
    paired_normal = form.get("paired_normal")
    paired_dir = form.get("paired_output_dir")
    if not paired_vcf and paired_dir:
        dna_manifest = latest_manifest(paired_dir, "dna")
        if dna_manifest:
            paired_vcf = manifest_vcf(dna_manifest)
            paired_normal = dna_manifest.get("normal") or ""
            if not paired_vcf:
                print("[PCGR] the linked DNA run finished but left no "
                      "variant calls on disk; reporting the fusions alone.")
        else:
            # The DNA half failed, was cancelled, or is still going. Report
            # the fusions on their own rather than waiting or failing: half
            # a result now beats no result.
            print(f"[PCGR] the linked DNA run in {paired_dir} has no "
                  f"manifest yet; reporting the fusions alone.")

    if paired_vcf and os.path.exists(paired_vcf):
        argv += ["--input-vcf", paired_vcf, "--lift-tags"]
        if form.get("vep_dir"):
            # Required only because a VCF is now present; a fusion-only
            # report needs no VEP cache, there being nothing to annotate.
            argv += ["--vep-dir", form["vep_dir"]]
        else:
            return None, ("the linked DNA run's VCF needs a VEP cache and "
                          "none is configured; report the two runs "
                          "separately")
        if not paired_normal:
            argv.append("--tumour-only")

    if form.get("pcgr_tumour_site"):
        argv += ["--tumour-site", str(form["pcgr_tumour_site"])]
    if form.get("min_fusion_reads"):
        argv += ["--fusion-min-split-reads", str(form["min_fusion_reads"])]
    return argv, None


def role_args(form):
    """
    The tumour/normal names for a run whose FASTQs are found, not given.

    With --auto-discover or --manifest the engine used to be passed no
    names at all and chose the tumour alphabetically -- which puts PT01-N
    before PT01-T, and so analysed the normal as the tumour. It now refuses
    to guess, so the names this app resolved on the form travel with the
    command.
    """
    out = []
    if form.get("tumour_sample"):
        out += ["--tumour-sample", form["tumour_sample"]]
    if form.get("normal_sample"):
        out += ["--normal-sample", form["normal_sample"]]
    elif form.get("tumour_only"):
        out += ["--tumour-only"]
    return out


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

    # Only meaningful when --reference names a genome instead of a FASTA,
    # but harmless otherwise, and passing it always means a run that falls
    # back to the name still shares one download and one index.
    if form.get("reference_dir"):
        argv += ["--reference-dir", form["reference_dir"]]

    if form.get("manifest"):
        argv += ["--manifest", form["manifest"]]
        argv += role_args(form)
    else:
        argv += ["--input-dir", form["input_dir"]]
        if form.get("auto_discover"):
            argv += ["--auto-discover"]
            argv += role_args(form)
        else:
            argv += ["--tumour-sample", form["tumour_sample"],
                     "--tumour-r1", form["tumour_r1"],
                     "--tumour-r2", form["tumour_r2"]]
            if form.get("normal_sample"):
                argv += ["--normal-sample", form["normal_sample"],
                         "--normal-r1", form["normal_r1"],
                         "--normal-r2", form["normal_r2"]]

    # The panel is passed through as well as being applied to the form.
    # Two reasons, and neither is redundancy:
    #   * the run manifest then records the assay as data -- which kit,
    #     which chemistry, which caveats -- so a reader a year later can
    #     see what the numbers belong to without this app's run record;
    #   * a profile carries settings the form has no field for (UMI
    #     layout, adapters, BQSR on a panel too small to fit it, caller
    #     passthroughs), and the pipeline applies those itself.
    # It cannot conflict with the form: the pipeline only fills in what its
    # command line did not state, and everything the form exposes is
    # stated below.
    if form.get("panel"):
        argv += ["--panel", form["panel"]]
    for extra in str(form.get("panel_file", "")).split():
        argv += ["--panel-file", extra]

    for flag, key in (("--cosmic", "cosmic"),
                      ("--dbsnp", "dbsnp"),
                      ("--germline-resource", "germline_resource"),
                      ("--panel-of-normals", "panel_of_normals"),
                      ("--contamination-resource", "contamination_resource"),
                      ("--msi-models", "msi_models"),
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

    # The coverage check. Optional by design: with no BED the pipeline says
    # so and carries on, because a panel BED is not always to hand and a
    # run without one is still a run.
    if form.get("coverage_bed"):
        argv += ["--coverage-bed", form["coverage_bed"]]
        if form.get("coverage_min_depth"):
            argv += ["--coverage-min-depth",
                     str(form["coverage_min_depth"])]
    if form.get("min_depth"):
        argv += ["--min-depth", str(form["min_depth"])]
    if form.get("min_allele_fraction"):
        argv += ["--min-allele-fraction", str(form["min_allele_fraction"])]
    if form.get("mutect2_extra_args"):
        argv += ["--mutect2-extra-args", str(form["mutect2_extra_args"])]

    if form.get("threads"):
        argv += ["--threads", str(form["threads"])]
    if form.get("two_pass"):
        argv += ["--two-pass"]
    if form.get("resume"):
        argv += ["--resume"]
    if form.get("dry_run"):
        argv += ["--dry-run"]

    # No PCGR options are passed, deliberately. PCGR lives in its own conda
    # environment and cannot run inside the pipeline's, so this app always
    # produces the report itself afterwards -- see build_pcgr_argv(). Asking
    # the pipeline for it as well only ever hurt: on a normal install its
    # step 12 announced "[SKIP] PCGR skipped (pcgr ... unavailable)" after
    # doing the FORMAT->INFO lift for a report it could not run, and where
    # pcgr HAD been installed alongside the pipeline it ran a second, hours-
    # long PCGR into the same output_dir/pcgr with settings this builder
    # never emitted -- no tumour site, no TMB/MSI/signature estimates -- so
    # which report survived depended on which finished last. Leaving
    # --pcgr-refdata-dir off makes the pipeline say so plainly instead.

    if form.get("skip_steps"):
        argv += ["--skip-steps"] + form["skip_steps"].split()

    return argv
