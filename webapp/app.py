#!/usr/bin/env python3
"""
app.py
======
Flask front end for the cancer DNA pipeline: submit a run with patient
details, watch it progress, and produce a PDF report from the result.

RUNNING IT
----------
    pip install flask            # the only dependency this app adds
    python webapp/app.py         # then open http://127.0.0.1:5000

SECURITY POSTURE -- PLEASE READ
-------------------------------
  This app handles patient-identifying information and can start processes
  that read and write anywhere the user running it can reach. It is built
  for a single analyst on a workstation or a private lab machine, and its
  defaults reflect that:

    * Binds to 127.0.0.1 ONLY. Nothing outside the machine can reach it
      unless you deliberately pass --host. There is NO authentication, so
      binding to 0.0.0.0 publishes both the patient details and the ability
      to execute pipeline runs to everyone on the network. The app warns
      loudly if you do.
    * Debug mode is off. Flask's debugger allows arbitrary code execution
      through the browser, which is unacceptable anywhere near patient data.
    * Server-side paths, not uploads. FASTQs are tens of GB; they are
      referenced by path on the machine doing the analysis rather than
      pushed through the browser.
    * Path inputs are constrained to configured roots (see safe_path) so a
      typo -- or a hostile form post -- cannot walk into /etc.
    * Patient details are written to the run directory and never passed to
      the pipeline, so they cannot leak into logs or manifests.

  It is not a multi-user system, it has no audit trail beyond the run logs,
  and it is not a validated clinical application.
"""

import json
import glob
import os
import re
import sys

# The pipeline scripts sit one directory up; report.py and jobs.py sit here.
HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(HERE)
sys.path.insert(0, HERE)

try:
    from flask import (Flask, abort, jsonify, redirect, render_template,
                       request, send_file, url_for)
except ImportError:  # pragma: no cover - dependency guidance
    sys.exit("Flask is not installed. Install it with:  pip install flask")

import jobs                                              # noqa: E402
from jobs import JobManager, build_pipeline_argv          # noqa: E402
from report import (PATIENT_FIELDS, build_report,         # noqa: E402
                    find_pcgr_outputs, load_pipeline_manifest,
                    patient_warnings)

SCRIPT_VARIANT = os.path.join(PIPELINE_DIR, "comprehensive_variant_calling.py")

app = Flask(__name__)

# Configured at startup by main(); see safe_path().
app.config["ALLOWED_ROOTS"] = []
app.config["RUNS_DIR"] = os.path.join(PIPELINE_DIR, "webapp_runs")
manager = None


# ---------------------------------------------------------------------------
# Path safety
# ---------------------------------------------------------------------------
def safe_path(raw, must_exist=False):
    """
    Resolve a user-supplied path and confine it to the allowed roots.

    The form asks for server-side paths, which means the browser can post
    any string it likes. Without this, "../../etc" would be handed to a
    subprocess. Symlinks are resolved before the check so a link inside an
    allowed root cannot be used to escape it.

    Returns the resolved path, or raises ValueError.
    """
    if not raw:
        raise ValueError("empty path")
    resolved = os.path.realpath(os.path.expanduser(str(raw).strip()))
    roots = app.config["ALLOWED_ROOTS"]
    if roots:
        if not any(resolved == r or resolved.startswith(r + os.sep)
                   for r in roots):
            raise ValueError(
                f"path is outside the allowed directories: {resolved}")
    if must_exist and not os.path.exists(resolved):
        raise ValueError(f"path does not exist: {resolved}")
    return resolved


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------
@app.route("/")
def index():
    manager.load_existing()
    return render_template("index.html",
                           jobs=[j.snapshot() for j in manager.all()],
                           busy=manager.busy(),
                           db_notice=database_notice())


@app.route("/new")
def new_run():
    return render_template("new_run.html",
                           patient_fields=PATIENT_FIELDS,
                           roots=app.config["ALLOWED_ROOTS"],
                           default_runs=app.config["RUNS_DIR"],
                           defaults=dict(discover_defaults(),
                                         **ANALYSIS_TOGGLES))


# ---------------------------------------------------------------------------
# Resource discovery
# ---------------------------------------------------------------------------
# A run submitted with these fields blank still succeeds -- the pipeline skips
# each step cleanly -- and produces an annotated VCF and no clinical report.
# That is rarely what anyone wants, and it is not obvious from the form, so
# the standard install locations are pre-filled when they exist. Every value
# stays editable and clearable; this changes the default, not the rules.
_DATA = os.path.expanduser("~/data")

DEFAULT_RESOURCES = {
    "dbsnp": f"{_DATA}/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    "germline_resource": f"{_DATA}/resources/hg38/af-only-gnomad.hg38.vcf.gz",
    "panel_of_normals": f"{_DATA}/resources/hg38/1000g_pon.hg38.vcf.gz",
    "contamination_resource":
        f"{_DATA}/resources/hg38/small_exac_common_3.hg38.vcf.gz",
    "cosmic":
        f"{_DATA}/resources/hg38/Cosmic_GenomeScreensMutant_v103_GRCh38"
        f".chr.vcf.gz",
    "vep_dir": f"{_DATA}/vep_cache",
    # Tumour-only MSI. PCGR omits MSI on a panel, so this is the only route
    # to an MSI answer for a targeted assay.
    "msi_models": f"{_DATA}/msisensor2/models_hg38",
}

# --known-indels takes several files, so it is kept apart from the map
# above.
DEFAULT_KNOWN_INDELS = [
    f"{_DATA}/resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    f"{_DATA}/resources/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz",
]


def discover_defaults():
    """
    Standard resource paths that actually exist on this machine.

    Missing entries are omitted rather than offered, so the form never
    proposes a path that would fail validation. The PCGR bundle is found by
    globbing, because its directory is named for the bundle release date and
    changes with every refresh.
    """
    found = {k: v for k, v in DEFAULT_RESOURCES.items() if os.path.exists(v)}

    bundles = sorted(glob.glob(os.path.join(_DATA, "pcgr", "[0-9]" * 8)))
    if bundles:
        found["pcgr_refdata_dir"] = bundles[-1]

    # Clinical reports are the point of running this at all, so the estimate
    # boxes start ticked. MSI is left off: PCGR honours it only for WGS/WES
    # tumour-normal runs and ignoring that just puts a misleading blank in
    # the report.
    indels = [p for p in DEFAULT_KNOWN_INDELS if os.path.exists(p)]
    if indels:
        found["known_indels"] = " ".join(indels)
    return found


# Analysis toggles, not resources. They are checkboxes, and an unticked
# checkbox is simply absent from the submission -- indistinguishable from a
# field nobody filled in. Defaulting them the way paths are defaulted made
# them impossible to switch OFF: unticking sent nothing, and the blank was
# filled straight back in. They are kept separate, and only defaulted for a
# submission that never rendered the form.
ANALYSIS_TOGGLES = {
    # A clinical report is the point of running this, so TMB and signature
    # fitting start on.
    "pcgr_estimate_tmb": "on",
    "pcgr_estimate_signatures": "on",
    # MSI starts off deliberately. PCGR honours it only for WGS/WES
    # tumour-control runs; on a targeted or tumour-only query it warns and
    # skips, so defaulting it on would promise a section the report will
    # not contain.
    "pcgr_estimate_msi": "",
}


def apply_resource_defaults(form):
    """
    Fill in every installed resource the submission left blank.

    Pre-filling the form was not enough. Those values are only
    suggestions: a field cleared by accident, or a submission that never
    rendered the form, silently drops a whole step and the run still
    reports success. This applies them SERVER-SIDE at submit time, so
    the full resource flow is what actually runs rather than what was
    offered.

    Opt out with the "minimal run" checkbox, which is the supported way
    to run without a resource -- rather than clearing a field and hoping
    someone notices. Returns the keys that were filled in, so the run
    record can state what was added instead of leaving the reader to
    infer it from a command line.
    """
    if form.get("minimal_run"):
        return []
    applied = []
    for key, value in discover_defaults().items():
        if not form.get(key):
            form[key] = value
            applied.append(key)

    # Toggles are only defaulted for a submission that did not render the
    # form -- an API call or a script. The form posts a marker, and when it
    # is present the checkboxes are taken exactly as submitted, so unticking
    # one actually turns it off.
    if not form.get("form_rendered"):
        for key, value in ANALYSIS_TOGGLES.items():
            if value and not form.get(key):
                form[key] = value
                applied.append(key)
    return applied


@app.route("/submit", methods=["POST"])
def submit():
    form = {k: v.strip() for k, v in request.form.items()}
    patient = {key: form.get(key, "") for key, _ in PATIENT_FIELDS}

    # Before validation, so anything filled in here is path-checked exactly
    # like a value the user typed.
    auto_applied = apply_resource_defaults(form)

    errors = []
    # --- validate the paths before anything is started -----------------
    resolved = {}
    try:
        resolved["output_dir"] = safe_path(form.get("output_dir"))
    except ValueError as exc:
        errors.append(f"Output directory: {exc}")

    try:
        resolved["reference"] = safe_path(form.get("reference"),
                                          must_exist=True)
    except ValueError as exc:
        # 'hg38' is a genome name the pipeline downloads, not a path.
        if form.get("reference", "").strip() == "hg38":
            resolved["reference"] = "hg38"
        else:
            errors.append(f"Reference: {exc}")

    if form.get("manifest"):
        try:
            resolved["manifest"] = safe_path(form["manifest"], must_exist=True)
        except ValueError as exc:
            errors.append(f"QC manifest: {exc}")
    else:
        try:
            resolved["input_dir"] = safe_path(form.get("input_dir"),
                                              must_exist=True)
        except ValueError as exc:
            errors.append(f"Input directory: {exc}")
        if not form.get("auto_discover"):
            for field in ("tumour_sample", "tumour_r1", "tumour_r2"):
                if not form.get(field):
                    errors.append(f"{field.replace('_', ' ')} is required "
                                  f"unless auto-discover or a manifest is used")
            if form.get("normal_sample") and not (form.get("normal_r1")
                                                  and form.get("normal_r2")):
                errors.append("A normal sample needs both normal R1 and R2")

    for optional in ("cosmic", "dbsnp", "germline_resource",
                     "panel_of_normals", "contamination_resource",
                     "msi_models", "intervals", "pcgr_refdata_dir",
                     "vep_dir"):
        if form.get(optional):
            try:
                resolved[optional] = safe_path(form[optional], must_exist=True)
            except ValueError as exc:
                errors.append(f"{optional.replace('_', ' ')}: {exc}")

    if not patient.get("patient_id") and not patient.get("specimen_id"):
        errors.append("Give at least a patient/MRN or a specimen ID, so the "
                      "report can be attributed to something")

    # Resolve auto-discover here rather than letting the pipeline do it. See
    # discover_pairs() for why an ambiguous auto-discover is dangerous.
    batch_pairs = []
    if form.get("auto_discover") and not form.get("manifest") and \
            not errors and resolved.get("input_dir"):
        pairs, warns = discover_pairs(resolved["input_dir"],
                                      bool(form.get("recursive")))
        for w in warns:
            errors.append(f"input scan: {w}")
        if not pairs:
            errors.append("Auto-discover found no FASTQ pairs in "
                          f"{resolved['input_dir']}")
        elif len(pairs) == 1:
            pass                       # unambiguous: one tumour-only sample
        elif form.get("batch_mode"):
            batch_pairs = pairs        # one run per sample, queued in turn
        elif form.get("second_is_matched_normal"):
            if len(pairs) != 2:
                errors.append(
                    f"'Second sample is the matched normal' needs exactly two "
                    f"FASTQ pairs, but {len(pairs)} were found: "
                    f"{', '.join(p['sample'] for p in pairs)}.")
        else:
            errors.append(
                f"{len(pairs)} samples found in {resolved['input_dir']}: "
                f"{', '.join(p['sample'] for p in pairs)}. Choose 'Batch' to "
                f"run each as its own patient, or tick 'second sample is the "
                f"matched normal' if these are a tumour/normal pair from ONE "
                f"person. Running them unmarked would call somatic variants "
                f"on the difference between two people.")

    if errors:
        return render_template("new_run.html", errors=errors, form=form,
                               patient=patient,
                               patient_fields=PATIENT_FIELDS,
                               roots=app.config["ALLOWED_ROOTS"],
                               default_runs=app.config["RUNS_DIR"],
                               defaults=dict(discover_defaults(),
                                             **ANALYSIS_TOGGLES)), 400

    merged = dict(form)
    merged.update(resolved)
    merged["auto_applied_resources"] = " ".join(auto_applied)
    # The pipeline runs in the cancer_pipeline environment, not in whichever
    # interpreter is serving this app: launching it with sys.executable is
    # why a webapp started outside that env failed with
    # "Required tools not found: bwa-mem2, samtools, gatk, fastp".
    pipeline_env = jobs.find_conda_env(jobs.PIPELINE_ENV_NAME)
    pipeline_python = jobs.env_python(pipeline_env) or sys.executable
    argv = build_pipeline_argv(pipeline_python, SCRIPT_VARIANT, merged)

    # PCGR runs as a follow-up in its own conda environment -- it cannot run
    # inside the pipeline's, so passing --pcgr-refdata-dir alone would leave
    # step 11 reporting itself as skipped and no report on disk.
    pcgr_form = merged if merged.get("pcgr_refdata_dir") else None

    # --- batch: one queued job per sample, each its own output directory ---
    if batch_pairs:
        base_out = merged["output_dir"]
        first = None
        for pair in batch_pairs:
            per = dict(merged)
            per["output_dir"] = sample_output_dir(base_out, pair["sample"])
            per["tumour_sample"] = pair["sample"]
            per["tumour_r1"] = pair["r1"]
            per["tumour_r2"] = pair["r2"]
            # Explicit samples, so auto-discover must not also be passed --
            # it would override them and reintroduce the ambiguity.
            per.pop("auto_discover", None)
            per["batch_of"] = str(len(batch_pairs))
            per["batch_base"] = base_out

            per_patient = dict(patient)
            # Each run is a distinct specimen; keep the operator's patient
            # details but make the sample it belongs to unambiguous.
            per_patient["specimen_id"] = (
                f"{patient.get('specimen_id') or ''} {pair['sample']}".strip())

            per_argv = build_pipeline_argv(pipeline_python, SCRIPT_VARIANT,
                                           per)
            job = manager.submit(
                SCRIPT_VARIANT, per_argv, per_patient, meta=per,
                pcgr_form=per if per.get("pcgr_refdata_dir") else None)
            first = first or job
        return redirect(url_for("job_view", job_id=first.id))

    job = manager.submit(SCRIPT_VARIANT, argv, patient, meta=merged,
                         pcgr_form=pcgr_form)
    return redirect(url_for("job_view", job_id=job.id))


# ---------------------------------------------------------------------------
# Sample discovery and batch submission
# ---------------------------------------------------------------------------
# --auto-discover assigns the FIRST pair it finds as the tumour and the
# SECOND as that tumour's matched normal, then ignores the rest, and says
# nothing about any of it. Point it at a directory holding several patients
# and it calls somatic variants on the genetic difference between two
# unrelated people: their differing germline variants are reported as
# somatic, their shared real mutations are subtracted, and the output looks
# entirely ordinary. PCGR will then tier it and produce a confident report.
#
# The webapp therefore never submits an ambiguous auto-discover run. It
# resolves the samples itself and passes them explicitly.


def discover_pairs(input_dir, recursive=False):
    """
    FASTQ pairs in a directory, using the pipeline's own pairing rules.

    Imported from the pipeline rather than reimplemented so that the two
    cannot disagree about what constitutes a pair -- lane merging and R1/R2
    token handling are subtle enough that a second implementation would
    drift. Returns (pairs, warnings); pairs is empty on any failure.
    """
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_cvc_discovery", SCRIPT_VARIANT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.find_paired_fastqs(
            input_dir, "*_R1_*.fastq*", "*_R2_*.fastq*",
            recursive=recursive)
    except Exception as exc:                # noqa: BLE001 - surfaced to UI
        return [], ["could not scan %s: %s" % (input_dir, exc)]


def sample_output_dir(base, sample):
    """Per-sample output directory for a batch run."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", sample)
    return os.path.join(base, safe)


# ---------------------------------------------------------------------------
# Database update notice
# ---------------------------------------------------------------------------
# Read from a file that check_db_updates.py writes out of band. No network
# call happens while serving a page: a request must not hang because
# ftp.ensembl.org is slow, and a run must never be delayed or altered by an
# upgrade being available. The notice is information, and nothing in the
# pipeline consults it.


def database_notice():
    """
    What the last database check found, or None.

    Returns None when no check has ever run, so a machine that never set up
    the cron entry shows nothing rather than a permanent scolding. A stale
    result is surfaced as stale rather than quietly presented as current --
    "checked three months ago, all well" is not the same claim as "all well".
    """
    try:
        sys.path.insert(0, os.path.dirname(SCRIPT_VARIANT))
        from check_db_updates import load_status
    except Exception:
        return None
    status = load_status()
    if not status:
        return None
    return {
        "checked": (status.get("checked_utc") or "")[:16].replace("T", " "),
        "age_days": status.get("age_days"),
        "stale": status.get("stale"),
        "needs_attention": status.get("needs_attention", 0),
        # Only the rows worth showing. "current" is the boring majority and
        # listing it buries the one line that matters.
        #
        # NOT named "items": Jinja resolves dict.items to the built-in method
        # before it looks for a key of that name, so `db_notice.items` in a
        # template silently yields a bound method instead of this list.
        "flagged": [r for r in status.get("results", [])
                    if r.get("state") in ("update", "changed", "manual")],
    }


def run_warnings(meta):
    """
    What this run gave up by leaving fields blank.

    The form accepts an empty resource set and the pipeline then skips those
    steps by design, each with one line in a log nobody reads. The result is
    a VCF that looks finished and is quietly much weaker than it should be --
    which is exactly what happens on a first run, when the person driving the
    form does not yet know which boxes matter. These are warnings, not
    errors: a bare run is legitimate, it just should not be silent.
    """
    notes = []
    tumour_only = not meta.get("normal_sample")

    if not meta.get("germline_resource") and tumour_only:
        notes.append(
            "No germline resource (gnomAD). In tumour-only mode this is the "
            "main defence against germline variants being reported as "
            "somatic, and without it Mutect2's germline filter barely fires.")
    if not meta.get("dbsnp") and not meta.get("known_indels"):
        notes.append(
            "No dbSNP or known indels, so BQSR was skipped: base quality "
            "scores were used as the instrument reported them.")
    if not meta.get("intervals"):
        notes.append(
            "No target intervals. On a capture panel or exome Mutect2 walks "
            "the whole genome and calls on off-target reads, which can "
            "outnumber the real calls by two orders of magnitude.")
    if not meta.get("contamination_resource"):
        notes.append(
            "No contamination resource, so filtering assumed a contamination "
            "fraction of zero rather than measuring one.")
    if not meta.get("panel_of_normals"):
        notes.append("No panel of normals: recurrent artefacts were not "
                     "subtracted.")
    if not meta.get("cosmic"):
        notes.append("No COSMIC VCF: calls carry no known-mutation "
                     "identifiers.")
    if not meta.get("pcgr_refdata_dir"):
        notes.append(
            "No PCGR reference bundle, so no clinical report was produced -- "
            "only an annotated VCF.")

    # Say so rather than leaving a blank section to be interpreted.
    if not meta.get("pcgr_estimate_tmb"):
        notes.append(
            "Tumour mutational burden was not computed. PCGR omits the "
            "section entirely rather than reporting zero.")
    if not meta.get("pcgr_estimate_signatures"):
        notes.append("Mutational signatures were not fitted.")
    if not meta.get("msi_models"):
        notes.append(
            "No MSI scoring: --msi-models was not supplied. PCGR does not "
            "cover this on a panel, so the report has no MSI answer at all.")
    if meta.get("pcgr_estimate_msi") and (
            meta.get("pcgr_assay") != "WGS" and
            meta.get("pcgr_assay") != "WES" or not meta.get("normal_sample")):
        notes.append(
            "MSI was requested but PCGR restricts it to WGS/WES "
            "tumour-control runs, so it was skipped. A panel needs a "
            "dedicated caller such as MSIsensor2.")
    if meta.get("pcgr_estimate_tmb") and not meta.get("normal_sample"):
        notes.append(
            "TMB was computed from a tumour-only call set. PCGR warns that "
            "this is unreliable -- residual germline variants inflate it -- "
            "so treat the figure as internal QC, not a reportable result.")
    return notes


def render_job(job, notice=None, notice_kind="ok"):
    """
    Render the job page.

    Shared by the view and by actions that report a result, so an
    action can show its outcome without a session -- this app has no
    secret key and does not need one for a single-analyst tool.
    """
    entries, total = jobs.survey_intermediates(
        job.meta.get("output_dir"))
    return render_template(
        "job.html", job=job.snapshot(),
        patient=job.patient,
        patient_fields=PATIENT_FIELDS,
        char_warnings=patient_warnings(job.patient),
        run_notes=run_warnings(job.meta),
        cleanup={"entries": entries, "total": total},
        human_bytes=jobs.human_bytes,
        notice=notice, notice_kind=notice_kind,
        batch_of=job.meta.get("batch_of"),
        auto_applied=[k for k in
                      (job.meta.get("auto_applied_resources") or
                       "").split() if k])


@app.route("/job/<job_id>")
def job_view(job_id):
    return render_job(manager.get(job_id) or abort(404))


@app.route("/api/job/<job_id>")
def job_api(job_id):
    """Polled by the job page for live status and log tail."""
    job = manager.get(job_id) or abort(404)
    return jsonify({"job": job.snapshot(), "log": job.tail()})


@app.route("/job/<job_id>/cancel", methods=["POST"])
def job_cancel(job_id):
    job = manager.get(job_id) or abort(404)
    job.cancel()
    return redirect(url_for("job_view", job_id=job_id))


@app.route("/job/<job_id>/report", methods=["POST"])
def job_report(job_id):
    """Generate (or regenerate) the PDF for a finished run."""
    job = manager.get(job_id) or abort(404)
    output_dir = job.meta.get("output_dir")
    if not output_dir:
        abort(400, "this run has no recorded output directory")

    pdf_path = os.path.join(job.run_dir, f"report_{job_id}.pdf")
    try:
        result = build_report(job.patient, job.snapshot(), output_dir,
                              pdf_path)
    except Exception as exc:  # noqa: BLE001 - surfaced to the user
        return render_template("job.html", job=job.snapshot(),
                               patient=job.patient,
                               patient_fields=PATIENT_FIELDS,
                               char_warnings=patient_warnings(job.patient),
                               report_error=f"{type(exc).__name__}: {exc}"), 500

    with open(os.path.join(job.run_dir, "report_meta.json"), "w",
              encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, default=str)
    return redirect(url_for("job_view", job_id=job_id))


@app.route("/job/<job_id>/report.pdf")
def job_report_download(job_id):
    job = manager.get(job_id) or abort(404)
    pdf_path = os.path.join(job.run_dir, f"report_{job_id}.pdf")
    if not os.path.exists(pdf_path):
        abort(404, "no report has been generated for this run yet")
    label = job.patient.get("patient_id") or job.patient.get("specimen_id") \
        or job_id
    safe_label = "".join(c if c.isalnum() or c in "-_" else "_"
                         for c in str(label))
    return send_file(pdf_path, mimetype="application/pdf",
                     as_attachment=True,
                     download_name=f"report_{safe_label}.pdf")


@app.route("/job/<job_id>/cleanup", methods=["POST"])
def job_cleanup(job_id):
    """
    Delete a finished run's intermediates on request.

    Refused while the job is running or queued: those files are its
    working set. Never automatic -- rebuilding them costs about an hour,
    so the decision belongs to whoever knows the run is done with.
    """
    job = manager.get(job_id) or abort(404)
    snap = job.snapshot()
    if snap.get("status") in ("running", "queued"):
        return render_job(job, "Cannot clean up while the run is %s."
                          % snap["status"], "error"), 409

    keep_bam = bool(request.form.get("keep_bam"))
    freed, removed, errs = jobs.cleanup_intermediates(
        job.meta.get("output_dir"), keep_bam=keep_bam)

    if errs:
        return render_job(job, "Cleanup problems: %s" % "; ".join(errs),
                          "error")
    if not removed:
        return render_job(job, "Nothing to clean up -- the intermediates "
                          "are already gone.")
    return render_job(job, "Freed %s by removing %s.%s"
                      % (jobs.human_bytes(freed), ", ".join(removed),
                         " Recalibrated BAM kept." if keep_bam else ""))


@app.route("/job/<job_id>/results")
def job_results(job_id):
    """Show what the run produced, and where PCGR's own report lives."""
    job = manager.get(job_id) or abort(404)
    output_dir = job.meta.get("output_dir") or ""
    manifest = load_pipeline_manifest(output_dir) if output_dir else {}
    pcgr = find_pcgr_outputs(output_dir) if output_dir else {"html": [],
                                                             "tsv": [],
                                                             "dir": ""}
    pdf_path = os.path.join(job.run_dir, f"report_{job_id}.pdf")
    return render_template("results.html", job=job.snapshot(),
                           manifest=manifest, pcgr=pcgr,
                           has_pdf=os.path.exists(pdf_path),
                           patient=job.patient,
                           patient_fields=PATIENT_FIELDS)


@app.route("/job/<job_id>/pcgr")
def job_pcgr(job_id):
    """
    Serve PCGR's own HTML report.

    PCGR writes a self-contained page; it is served from the run's own
    directory only, never from an arbitrary path supplied by the request.
    """
    job = manager.get(job_id) or abort(404)
    output_dir = job.meta.get("output_dir") or abort(404)
    pcgr = find_pcgr_outputs(output_dir)
    if not pcgr["html"]:
        abort(404, "no PCGR report found for this run")
    return send_file(pcgr["html"][0])


@app.route("/job/<job_id>/log")
def job_log(job_id):
    job = manager.get(job_id) or abort(404)
    if not os.path.exists(job.log_path):
        abort(404, "no log for this run")
    return send_file(job.log_path, mimetype="text/plain")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Web interface for the cancer DNA pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--host", default="127.0.0.1",
                        help="Interface to bind. The default keeps the app "
                             "on this machine. There is NO authentication, "
                             "so any other value exposes patient details and "
                             "run submission to the network.")
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--runs-dir", default=None,
                        help="Where run directories, patient records and "
                             "generated PDFs are kept.")
    parser.add_argument("--allow-root", action="append", default=[],
                        metavar="DIR",
                        help="Restrict every path the form accepts to this "
                             "directory (repeatable). Strongly recommended: "
                             "without it, any path readable by this user can "
                             "be submitted.")
    parser.add_argument("--debug", action="store_true",
                        help="Flask debug mode. NEVER use this with real "
                             "patient data: the debugger allows arbitrary "
                             "code execution from the browser.")
    args = parser.parse_args()

    global manager
    runs_dir = os.path.abspath(
        os.path.expanduser(args.runs_dir or app.config["RUNS_DIR"]))

    # Fail with an explanation rather than a traceback. Getting this wrong
    # is the most likely first-run mistake -- a plausible-looking path like
    # /data usually needs root, or does not exist at all.
    try:
        os.makedirs(runs_dir, exist_ok=True)
    except PermissionError:
        sys.exit(
            f"Cannot create the run directory: {runs_dir}\n"
            f"  You do not have permission to create it.\n"
            f"  Choose a directory you own, for example:\n"
            f"      --runs-dir ~/cancer_runs\n"
            f"  or omit --runs-dir to use the default "
            f"({app.config['RUNS_DIR']}).")
    except OSError as exc:
        sys.exit(f"Cannot create the run directory: {runs_dir}\n  {exc}")
    if not os.access(runs_dir, os.W_OK):
        sys.exit(f"The run directory is not writable: {runs_dir}\n"
                 f"  Runs, patient records and generated reports are written "
                 f"here, so the app cannot start without it.")

    # An --allow-root that does not exist would silently reject every path
    # the form accepts, which looks like a mysterious validation bug.
    roots, bad = [], []
    for raw in args.allow_root:
        resolved = os.path.realpath(os.path.expanduser(raw))
        (roots if os.path.isdir(resolved) else bad).append(resolved)
    if bad:
        sys.exit("These --allow-root directories do not exist:\n" +
                 "".join(f"      {b}\n" for b in bad) +
                 "  Every submitted path is confined to these directories, so "
                 "a wrong one\n  would reject all input. Point it at where "
                 "your FASTQs and references live,\n  e.g. --allow-root "
                 "~/Coding++  (repeat the flag for several locations).")

    app.config["RUNS_DIR"] = runs_dir
    app.config["ALLOWED_ROOTS"] = roots
    manager = JobManager(runs_dir)
    manager.load_existing()

    print(f"[INFO] Pipeline scripts : {PIPELINE_DIR}")
    print(f"[INFO] Run directory    : {runs_dir}")
    if app.config["ALLOWED_ROOTS"]:
        print(f"[INFO] Allowed paths    : "
              f"{', '.join(app.config['ALLOWED_ROOTS'])}")
    else:
        print("[WARN] No --allow-root given: the form will accept any path "
              "this user can read. Pass --allow-root to constrain it.")
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"[WARN] Binding to {args.host}: this app has NO "
              f"authentication, so patient details and the ability to start "
              f"pipeline runs are exposed to anyone who can reach this port.")
    if args.debug:
        print("[WARN] Debug mode is ON. The Flask debugger permits arbitrary "
              "code execution from the browser -- never use it with real "
              "patient data.")
    print(f"[INFO] Open http://{args.host}:{args.port}")

    app.run(host=args.host, port=args.port, debug=args.debug,
            threaded=True)


if __name__ == "__main__":
    main()
