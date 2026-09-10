#!/usr/bin/env python3
"""
maintenance.py
==============
Database checks and upgrades, started from the web interface.

WHAT THIS IS FOR
----------------
  check_db_updates.py answers "is anything out of date?" and stops there,
  on purpose. Answering it required someone to run the script, so the
  notice on the front page was only as fresh as the last time a human
  remembered -- and acting on it meant leaving the browser, activating an
  environment and re-running the installer with the right flags.

  This module closes both gaps: the app checks on startup and on demand,
  and can run the installer step that performs an upgrade.

WHY IT IS A SEPARATE, SERIAL RUNNER
-----------------------------------
  * ONE MAINTENANCE TASK AT A TIME, and never while a pipeline run is in
    progress. That rule is enforced by the caller (app.py), which knows
    about the run queue. Both constraints are about the same thing: an
    upgrade replaces the files a run reads, and re-downloading tens of
    gigabytes alongside a run competes for the same disk and network.
  * UPGRADES ARE NOT HOUSEKEEPING. A new VEP cache changes transcript
    sets, and a new COSMIC changes identifiers, so a report regenerated
    afterwards can legitimately disagree with one issued last week. The UI
    says so before it starts, and the log records exactly what ran.
  * THE INSTALLER REMAINS THE ONLY THING THAT WRITES TO THE DATA DIRECTORY.
    Nothing here downloads or deletes anything itself; it runs
    install_pipeline.py with a single --only step and reports what happened.
    There is one place that knows how each database is installed, and this
    is not it.
"""

import os
import subprocess
import sys
import threading
from collections import deque
from datetime import datetime, timezone

TAIL_LINES = 400

# What the front page can act on. Each entry maps a row from
# check_db_updates.py to the installer step that would upgrade it.
#
# Only steps that genuinely perform the upgrade are listed. The SnpEff row
# reports the bioconda PACKAGE version while the installer's snpeff step
# fetches the hg38 DATABASE -- two different things -- so that row carries
# an explicit note instead of a button that would do something other than
# what the row describes.
UPDATERS = {
    "Ensembl VEP cache": {
        "step": "pcgr-data",
        "label": "PCGR bundle and VEP cache",
        # Without a target release the installer would refetch the release
        # it is pinned to -- the very one the check just called outdated.
        "fields": [
            {"name": "vep_release", "flag": "--vep-release",
             "label": "VEP cache release to install", "from_latest": True},
            {"name": "pcgr_bundle", "flag": "--pcgr-bundle",
             "label": "PCGR bundle release (YYYYMMDD, blank to keep)"},
        ],
        "warning": "Re-downloads up to about 31 GB and can take hours. Two "
                   "cautions: a new VEP cache changes transcript sets, so "
                   "reports produced afterwards can differ from earlier ones "
                   "for the same variants; and the cache release must match "
                   "the Ensembl VEP that PCGR's own environment installs -- "
                   "a newer cache with an older VEP is refused by VEP "
                   "itself, so upgrade PCGR first or together with it.",
    },
    "PCGR": {
        "step": "pcgr-envs",
        "label": "PCGR and pcgrr conda environments",
        "fields": [
            {"name": "pcgr_version", "flag": "--pcgr-version",
             "label": "PCGR version to install", "from_latest": True},
        ],
        "warning": "Rebuilds both environments from upstream's pinned lock "
                   "files for that version (about 6 GB). PCGR and its "
                   "reference bundle are versioned together, and the bundle "
                   "release cannot be discovered automatically: check the "
                   "release notes for the bundle that pairs with this "
                   "version, then update the bundle and VEP cache too. A "
                   "mismatched pair is the classic failed first run.",
    },
    "MSIsensor2 models": {
        "step": "msi-models",
        "label": "MSIsensor2 hg38 models",
        "warning": "Re-fetches the models (about 250 MB). MSI scores from "
                   "before and after are not directly comparable.",
    },
    "GATK resource VCFs": {
        "step": "resources",
        "label": "GATK resource VCFs",
        "warning": "Re-downloads about 6.4 GB (dbSNP, gnomAD, the panel of "
                   "normals, known indels and the contamination resource). "
                   "A 'changed' row usually means the local copy is "
                   "truncated, which is worth fixing before the next run.",
    },
    "SnpEff": {
        "step": "snpeff",
        "label": "SnpEff hg38 database",
        "warning": "Re-downloads the hg38 DATABASE (about 450 MB). It does "
                   "NOT upgrade the SnpEff package itself: that is a conda "
                   "operation -- `conda update -n cancer_pipeline snpeff` -- "
                   "and the database should be refreshed after it.",
    },
    "COSMIC": {
        "step": "cosmic",
        "label": "COSMIC VCF",
        "needs_file": True,
        "warning": "COSMIC needs a registered account, so it cannot be "
                   "fetched here. Download the GRCh38 VCF yourself, then "
                   "give its path: the installer renames its contigs to UCSC "
                   "style and indexes it. A new COSMIC release changes "
                   "identifiers, so reports will not match older ones.",
    },
}


class Task:
    """One maintenance job: a short sequence of commands and their log."""

    def __init__(self, kind, label, stages, log_path):
        self.kind = kind                 # "check" | "update"
        self.label = label
        self.stages = stages             # [(stage label, argv)]
        self.log_path = log_path
        self.stage_label = stages[0][0] if stages else ""
        self.status = "running"          # running | finished | failed
        self.returncode = None
        self.error = None
        self.started = datetime.now(timezone.utc).isoformat()
        self.finished = None
        self._tail = deque(maxlen=TAIL_LINES)
        self._lock = threading.Lock()

    def snapshot(self):
        with self._lock:
            return {
                "kind": self.kind,
                "label": self.label,
                "stage": self.stage_label,
                "status": self.status,
                "returncode": self.returncode,
                "error": self.error,
                "started": self.started,
                "finished": self.finished,
                "log": os.path.basename(self.log_path),
            }

    def tail(self, limit=TAIL_LINES):
        with self._lock:
            return list(self._tail)[-limit:]


class Maintenance:
    """
    Runs one database task at a time, in the background.

    Failures are recorded, never raised: this is a side activity, and a
    failed upgrade must leave the app serving as before -- with the
    installer's own output on screen to say why it failed.
    """

    def __init__(self, log_dir, script_dir, data_dir, python_exe=None):
        self.log_dir = log_dir
        self.script_dir = script_dir
        self.data_dir = data_dir
        # The scripts this runs use only the standard library, so whichever
        # interpreter is serving the app can run them. install_pipeline.py
        # locates conda itself.
        self.python_exe = python_exe or sys.executable
        self._lock = threading.Lock()
        self._current = None
        self._last = None

    # -- state ------------------------------------------------------------
    def busy(self):
        with self._lock:
            return self._current is not None

    def task(self):
        with self._lock:
            return self._current or self._last

    def snapshot(self):
        task = self.task()
        return task.snapshot() if task else None

    def tail(self, limit=TAIL_LINES):
        task = self.task()
        return task.tail(limit) if task else []

    # -- the two things it can do ----------------------------------------
    def check_argv(self):
        # --print because the log is what the page shows while this runs;
        # without it the checker is silent on stdout (it is built for cron)
        # and the panel would sit empty until the status file appeared.
        return [self.python_exe,
                os.path.join(self.script_dir, "check_db_updates.py"),
                "--data-dir", self.data_dir, "--print"]

    def update_argv(self, step, cosmic=None, extra=None):
        argv = [self.python_exe,
                os.path.join(self.script_dir, "install_pipeline.py"),
                "--only", step, "--force",
                "--data-dir", self.data_dir,
                "--repo", self.script_dir,
                # Colour codes would be written into the log and shown as
                # escape sequences in the browser.
                "--no-colour"]
        for flag, value in (extra or []):
            argv += [flag, value]
        if cosmic:
            argv += ["--cosmic", cosmic]
        return argv

    def start_check(self, reason="requested"):
        """Refresh the status file. Cheap, read-only, network-bound."""
        return self._start("check", f"Database check ({reason})",
                           [("checking every source", self.check_argv())])

    def start_update(self, name, cosmic=None, extra=None):
        """
        Run the installer step that upgrades one database.

        The check is re-run afterwards in the same task, so the front page
        reflects what is now installed rather than what was found before
        the upgrade -- otherwise the panel keeps advertising an update that
        has just been applied.
        """
        spec = UPDATERS.get(name)
        if not spec:
            return None, f"{name} cannot be updated from here"
        if spec.get("needs_file") and not cosmic:
            return None, f"{name} needs the path to a file you downloaded"
        return self._start(
            "update", f"Update: {spec['label']}",
            [(f"installer --only {spec['step']}",
              self.update_argv(spec["step"], cosmic=cosmic, extra=extra)),
             ("re-checking what is installed", self.check_argv())])

    # -- execution --------------------------------------------------------
    def _start(self, kind, label, stages):
        with self._lock:
            if self._current is not None:
                return None, (f"{self._current.label} is already running; "
                              f"wait for it to finish")
            os.makedirs(self.log_dir, exist_ok=True)
            stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
            task = Task(kind, label, stages,
                        os.path.join(self.log_dir, f"{stamp}-{kind}.log"))
            self._current = task
        threading.Thread(target=self._run, args=(task,), daemon=True).start()
        return task, None

    def _run(self, task):
        try:
            with open(task.log_path, "a", encoding="utf-8") as log:
                log.write(f"### {task.label}\n")
                for stage_label, argv in task.stages:
                    with task._lock:
                        task.stage_label = stage_label
                    log.write(f"\n### {stage_label}\n### {' '.join(argv)}\n")
                    log.flush()
                    code = self._stream(argv, log, task)
                    task.returncode = code
                    if code != 0:
                        task.status = "failed"
                        task.error = (f"{stage_label} exited {code} -- see the "
                                      f"log for what the installer reported")
                        break
                else:
                    task.status = "finished"
        except Exception as exc:            # noqa: BLE001 - shown in the UI
            task.status = "failed"
            task.error = f"{type(exc).__name__}: {exc}"
        finally:
            task.finished = datetime.now(timezone.utc).isoformat()
            with self._lock:
                self._last = task
                self._current = None

    def _stream(self, argv, log, task):
        """Run one command, tee-ing its output to the log and the tail."""
        proc = subprocess.Popen(
            argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, cwd=self.script_dir,
            # Its own process group, so a future cancel can signal the whole
            # tree rather than orphaning conda's children.
            start_new_session=True)
        for line in proc.stdout:
            line = line.rstrip("\n")
            log.write(line + "\n")
            log.flush()
            with task._lock:
                task._tail.append(line)
        proc.wait()
        return proc.returncode
