#!/usr/bin/env python3
# Created by Brainstorm, 2026.
"""
app.py
======
Flask front end for the Cancer Genomics Pipelines: submit a run with patient
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
    * Database upgrades are a deliberate, guarded action. The app can run
      install_pipeline.py to replace an installed database -- tens of GB,
      and a change to what every later report says -- so it refuses while
      a run is queued or in progress, requires a typed confirmation, and
      runs only the single installer step named. There is no route that
      deletes anything, and none that runs an arbitrary command.

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

import artefacts                                        # noqa: E402
import jobs                                              # noqa: E402
import maintenance                                        # noqa: E402
from jobs import JobManager, build_pipeline_argv          # noqa: E402
from report import (PATIENT_FIELDS, build_report,         # noqa: E402
                    find_pcgr_outputs, load_pipeline_manifest,
                    patient_warnings)
from rna_report import build_rna_report                   # noqa: E402

SCRIPT_VARIANT = os.path.join(PIPELINE_DIR, "comprehensive_variant_calling.py")

# Panel profiles live beside the pipeline scripts. Imported optionally: the
# app is fully usable without them (every setting a profile carries is also
# an ordinary form field), so a bundle missing the module loses the panel
# dropdown rather than the ability to start a run.
sys.path.insert(0, PIPELINE_DIR)
try:
    import panel_profiles                              # noqa: E402
except ImportError:                                     # pragma: no cover
    panel_profiles = None

app = Flask(__name__)

# Configured at startup by main(); see safe_path().
app.config["ALLOWED_ROOTS"] = []
app.config["RUNS_DIR"] = os.path.join(PIPELINE_DIR, "webapp_runs")
# Where the installer put the references and databases. Same default as
# _DATA below, which is what the form's resource paths are built from;
# main() can point both at another location with --data-dir.
app.config["DATA_DIR"] = os.path.expanduser("~/data")
manager = None
# Database checks and upgrades. Built in main() because it needs the run
# directory for its logs; the module-level None keeps import-time use (and
# the test client) honest about that.
db_tasks = None


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
                           db_notice=database_notice(),
                           db_task=db_tasks.snapshot() if db_tasks else None)


@app.route("/new")
def choose_pathway():
    """
    The first question: which pipeline, or both.

    WHY THIS PAGE EXISTS
      Before it, the interface had two entry points in the navigation bar
      and expected the operator to already know which one their assay
      needed. That is a reasonable expectation of a bioinformatician and an
      unreasonable one of a pathologist who has been handed a kit -- and it
      handled the commonest modern case worst of all, because a hybrid
      panel (TSO500, Oncomine Comprehensive, Archer) is ONE specimen whose
      two libraries have to be submitted as two runs and then linked, from
      memory, or the combined report never happens.

      So the choice is asked once, in the operator's own vocabulary --
      what am I looking for? -- rather than in the pipeline's.
    """
    return render_template("choose_pathway.html",
                           hybrid_pairs=hybrid_panel_pairs())


def hybrid_panel_pairs():
    """
    Profiles that are one half of a two-library kit, with their partner.

    Read off the profiles' own pairs_with field, so a site that adds its
    own hybrid kit gets it offered here without touching the app.
    """
    if panel_profiles is None:
        return []
    registry = panel_registry()
    pairs = []
    for profile in registry.values():
        partner_id = panel_profiles.partner_id(profile)
        if not partner_id or panel_profiles.assay_type(profile) != "dna":
            continue                      # list each kit once, DNA side
        partner = registry.get(partner_id)
        if not partner:
            continue
        pairs.append({
            "dna": profile["id"], "rna": partner["id"],
            "name": profile["name"].split("(")[0].strip(),
            "manufacturer": profile["manufacturer"],
        })
    return sorted(pairs, key=lambda p: (p["manufacturer"], p["name"]))


@app.route("/new/dna")
def new_run():
    return render_template("new_run.html",
                           patient_fields=PATIENT_FIELDS,
                           roots=app.config["ALLOWED_ROOTS"],
                           default_runs=app.config["RUNS_DIR"],
                           defaults=dict(discover_defaults(),
                                         **ANALYSIS_TOGGLES),
                           panels=panel_options(),
                           panel_settings=panel_settings_json(),
                           tumour_sites=TUMOUR_SITES)


# ---------------------------------------------------------------------------
# Panel profiles
# ---------------------------------------------------------------------------
# The form has a field for every setting a profile carries, which is
# exactly the problem: filling in nine of them correctly, per kit, from
# memory, is what a laboratory running four different assays gets wrong.
# The dropdown names the assay instead, and the SERVER applies it -- see
# apply_panel_to_form() for why the browser doing it as well is not enough.


def panel_registry():
    """Every known profile, or {} when the module is not importable."""
    if panel_profiles is None:
        return {}
    try:
        return panel_profiles.available_panels(
            warn=lambda msg: app.logger.warning("panel profile: %s", msg))
    except Exception:                     # noqa: BLE001 - never break the form
        return {}


def panel_options():
    """
    The dropdown's contents: (id, label, manufacturer) sorted for humans.

    Generic shapes first -- they are where somebody with an unlisted kit
    should start, and burying them under the vendor names makes the page
    look as though only listed kits are supported.
    """
    options = []
    for profile in panel_registry().values():
        options.append({
            "id": profile["id"],
            "name": profile["name"],
            "manufacturer": profile["manufacturer"],
            "chemistry": profile["chemistry"],
            "summary": profile.get("summary", ""),
            "notes": profile.get("notes", []),
            "requires": profile.get("requires", []),
            "nominal_target_size_mb": profile.get("nominal_target_size_mb"),
        })
    options.sort(key=lambda o: (o["manufacturer"] != "any",
                                o["manufacturer"].lower(), o["id"]))
    return options


def panel_settings_json():
    """
    Every profile's settings as JSON, for the form's live preview.

    The browser uses this only to SHOW what a choice implies while the
    operator is still deciding. It is not how the settings are applied --
    that happens server-side in apply_panel_to_form(), because a value the
    browser filled in is a value a stale page, a disabled script or a
    scripted POST can silently omit.
    """
    return json.dumps({p["id"]: {"settings": p.get("settings", {}),
                                 "notes": p.get("notes", []),
                                 "requires": p.get("requires", []),
                                 "summary": p.get("summary", "")}
                       for p in panel_registry().values()})


# Profile keys the form does not expose. They are still applied to the run,
# via --panel on the pipeline's own command line, but there is no box to
# fill in with them.
_PANEL_KEYS_NOT_ON_FORM = frozenset({
    "adapter_r1", "adapter_r2", "min_read_length", "cut_right", "correction",
    "umi_loc", "umi_len", "umi_skip", "foldback_min_match",
    "mutect2_extra_args", "skip_bqsr",
})


def apply_panel_to_form(form):
    """
    Fill the submitted form's blanks from the chosen profile, server-side.

    This mirrors apply_resource_defaults(), and for the same reason: a
    value that only ever existed in the browser is a value that can go
    missing without anybody noticing. A page left open while the profile
    changed, a script posting to /submit, a field cleared by accident --
    each produces a submission that looks deliberate and is not.

    A field the operator actually filled in is never touched, so the
    precedence is the same one the command line uses: what you typed beats
    the profile beats the default.

    Returns the list of field names that were filled in, for the run
    record, so it states what was added instead of leaving a later reader
    to work it out from a command line.
    """
    panel_id = form.get("panel")
    if not panel_id or panel_profiles is None:
        return []
    try:
        profile = panel_profiles.resolve_panel(panel_id, panel_registry())
    except Exception:                     # noqa: BLE001 - validated in submit
        return []

    # One BED, two fields: the run needs both a calling restriction and a
    # coverage statement, and they are the same file in all but unusual
    # laboratories. Asking for it twice is how one of them ends up empty.
    applied = []
    decided = []
    bed = form.get("panel_bed") or form.get("intervals") or \
        form.get("coverage_bed")
    if bed:
        for field in ("intervals", "coverage_bed"):
            if not form.get(field):
                form[field] = bed
                applied.append(field)

    for key, value in (profile.get("settings") or {}).items():
        if value is None or key in _PANEL_KEYS_NOT_ON_FORM:
            continue
        if key in ("intervals", "coverage_bed"):
            continue                      # handled above, from the BED field
        if key == "skip_steps":
            # Merged, never replaced: a profile that skips dedup and an
            # operator who skipped qc must get a run that skips both.
            current = [s for s in str(form.get("skip_steps", "")).split() if s]
            added = [s for s in value if s not in current]
            if added:
                form["skip_steps"] = " ".join(current + added)
                applied.append("skip_steps")
            continue
        if isinstance(value, bool):
            # Checkbox semantics: an unticked box is simply absent from the
            # submission, so a profile can only ever turn one ON. The
            # decision is RECORDED either way, so apply_resource_defaults()
            # does not helpfully turn back on something this profile
            # deliberately left off.
            #
            # And only for a submission that did NOT render the form. When
            # the form was rendered, its script has already ticked the
            # boxes this profile wants, the operator has seen them, and
            # what came back is their decision -- re-ticking one here would
            # make it impossible to untick.
            decided.append(key)
            if value and not form.get(key) and not form.get("form_rendered"):
                form[key] = "on"
                applied.append(key)
            continue
        if not form.get(key):
            form[key] = str(value)
            applied.append(key)
    form["_panel_toggles"] = " ".join(decided)
    return applied


# ---------------------------------------------------------------------------
# The RNA branch
# ---------------------------------------------------------------------------
# A separate form, a separate submit handler and a separate argv builder,
# mirroring the separation in the scripts themselves. Almost nothing is
# shared: the DNA form's vocabulary -- target intervals, padding, allele
# fractions, known sites, panel of normals, PCGR -- has no meaning on an RNA
# run, and one form carrying both would be a form most of whose fields are
# inapplicable whichever assay you picked.
#
# What IS shared is deliberate: the patient details, the path safety, the
# serial job queue and the panel registry. A DNA and an RNA library from one
# specimen are two runs that must be readable as one result, which is what
# the paired_run_id field below records.

SCRIPT_FUSION = os.path.join(PIPELINE_DIR, "fusion_calling.py")


def rna_panel_options():
    """Only the RNA profiles, for the RNA form's dropdown."""
    if panel_profiles is None:
        return []
    return [o for o in panel_options()
            if o["chemistry"] in panel_profiles.RNA_CHEMISTRIES]


def default_rna_resources(data=None):
    """
    The RNA reference data the installer puts on disk, if it is there.

    Same idea as discover_defaults() on the DNA side: the form is filled in
    from what is actually installed rather than from a guess, so an operator
    does not have to remember a path that only ever has one right value.
    """
    data = data or app.config["DATA_DIR"]
    found = {}

    gencode = os.path.join(data, "references", "gencode")
    if os.path.isdir(gencode):
        # Newest release present. A laboratory that has installed two is
        # mid-upgrade, and the newer one is what a new run should use --
        # but the field is editable, because the STAR index decides which
        # annotation is actually correct for a run.
        gtfs = sorted(glob.glob(os.path.join(gencode, "*.gtf")))
        if gtfs:
            found["gtf"] = gtfs[-1]

    indexes = sorted(glob.glob(os.path.join(data, "references",
                                            "star_hg38_*")))
    indexes = [d for d in indexes
               if os.path.exists(os.path.join(d, "SAindex"))]
    if indexes:
        found["star_index"] = indexes[-1]
        # Read length is in the directory name because that is what makes
        # two indexes different; recover it so the form can show what this
        # index is for.
        tail = os.path.basename(indexes[-1]).rsplit("_", 1)[-1]
        if tail.isdigit():
            found["read_length"] = tail

    reference = os.path.join(data, "references", "hg38",
                             "Homo_sapiens_assembly38.fasta")
    if os.path.exists(reference):
        found["reference"] = reference
    return found


@app.route("/new/rna")
def new_rna_run():
    return render_template("new_rna_run.html",
                           patient_fields=PATIENT_FIELDS,
                           roots=app.config["ALLOWED_ROOTS"],
                           default_runs=app.config["RUNS_DIR"],
                           defaults=default_rna_resources(),
                           panels=rna_panel_options(),
                           panel_settings=panel_settings_json(),
                           dna_runs=finished_dna_runs())


def finished_dna_runs():
    """
    DNA runs this specimen's RNA run could be paired with.

    Offered rather than typed, because the pairing is the point: a DNA and
    an RNA library from one specimen are one assay reported together, and a
    free-text field for the partner run is a free-text field somebody
    mistypes.
    """
    rows = []
    try:
        manager.load_existing()
        for job in manager.all():
            snap = job.snapshot()
            if snap.get("assay") == "rna":
                continue
            label = (job.patient.get("specimen_id")
                     or job.patient.get("patient_id") or snap["id"])
            rows.append({"id": snap["id"], "label": label,
                         "status": snap["status"]})
    except Exception:                     # noqa: BLE001 - never break the form
        return []
    return list(reversed(rows))[:40]


@app.route("/submit/rna", methods=["POST"])
def submit_rna(validate_only=False):
    """
    Validate and queue an RNA fusion run.

    `validate_only` stops before queueing and returns None when the
    submission is good -- see submit() for why the hybrid route needs it.

    Deliberately parallel to submit() rather than folded into it: the two
    validate different things, and the one field they must NOT share is the
    assay, because a DNA form posted to the RNA handler by accident would
    otherwise start an RNA run with DNA settings and no complaint.
    """
    form = {k: v.strip() for k, v in request.form.items()}
    patient = {key: form.get(key, "") for key, _ in PATIENT_FIELDS}

    panel_applied = apply_panel_to_form(form)
    # Fill in the installed reference data for anything left blank, the
    # same way the DNA form fills in its databases.
    for key, value in default_rna_resources().items():
        if not form.get(key):
            form[key] = value
            panel_applied.append(key)

    errors = []
    resolved = {}

    if form.get("panel"):
        if panel_profiles is None:
            errors.append(
                "A panel was chosen but panel_profiles.py is not importable.")
        elif form["panel"] not in panel_registry():
            errors.append(f"Unknown panel profile {form['panel']!r}.")
        elif panel_profiles.assay_type(
                panel_registry()[form["panel"]]) != "rna":
            # The guard that matters. A DNA profile here would configure
            # duplicate marking and allele-fraction floors for a pipeline
            # that has neither, and the run would start regardless.
            errors.append(
                f"{form['panel']!r} is a DNA profile. This form starts an "
                f"RNA fusion run; pick an RNA panel, or use the DNA form.")

    try:
        resolved["output_dir"] = safe_path(form.get("output_dir"))
    except ValueError as exc:
        errors.append(f"Output directory: {exc}")

    for field, label, required in (("reference", "Reference genome", True),
                                   ("gtf", "Annotation GTF", True),
                                   ("star_index", "STAR index", True),
                                   ("arriba_resources",
                                    "Arriba reference files", False)):
        if not form.get(field):
            if required:
                errors.append(
                    f"{label} is required. The installer puts it on disk "
                    f"with 'install_pipeline.py --with-rna'.")
            continue
        try:
            resolved[field] = safe_path(form[field], must_exist=True)
        except ValueError as exc:
            errors.append(f"{label}: {exc}")

    if form.get("manifest"):
        try:
            resolved["manifest"] = safe_path(form["manifest"],
                                             must_exist=True)
        except ValueError as exc:
            errors.append(f"QC manifest: {exc}")
    else:
        try:
            resolved["input_dir"] = safe_path(form.get("input_dir"),
                                              must_exist=True)
        except ValueError as exc:
            errors.append(f"Input directory: {exc}")
        if not form.get("auto_discover"):
            if not form.get("sample"):
                errors.append("Sample name is required unless auto-discover "
                              "or a manifest is used.")
            if not form.get("r1"):
                errors.append("Read 1 FASTQ is required.")
            for field in ("r1", "r2"):
                if form.get(field):
                    try:
                        resolved[field] = safe_path(form[field],
                                                    must_exist=True)
                    except ValueError as exc:
                        errors.append(f"{field.upper()}: {exc}")

    if not patient.get("patient_id") and not patient.get("specimen_id"):
        errors.append("Give at least a patient/MRN or a specimen ID, so the "
                      "report can be attributed to something.")

    if errors:
        return render_template("new_rna_run.html", errors=errors, form=form,
                               patient=patient,
                               patient_fields=PATIENT_FIELDS,
                               roots=app.config["ALLOWED_ROOTS"],
                               default_runs=app.config["RUNS_DIR"],
                               defaults=default_rna_resources(),
                               panels=rna_panel_options(),
                               panel_settings=panel_settings_json(),
                               dna_runs=finished_dna_runs()), 400

    if validate_only:
        return None                       # everything checked out

    merged = dict(form)
    merged.update(resolved)
    merged["assay"] = "rna"
    merged["panel_applied_fields"] = " ".join(panel_applied)
    if form.get("panel") and panel_profiles is not None:
        profile = panel_registry().get(form["panel"])
        if profile:
            merged["panel_label"] = f"{profile['name']} ({profile['id']})"
            merged["panel_chemistry"] = profile["chemistry"]
            merged["panel_notes"] = profile.get("notes", [])

    # PCGR 2.x takes RNA fusions as a molecular input in their own right --
    # verified against the installed version's own validator -- so an RNA
    # run gets a real clinical interpretation, not only the pipeline's
    # ranked table. The bundle is filled in from what is installed, exactly
    # as the DNA form does it.
    for key, value in discover_defaults().items():
        if key in ("pcgr_refdata_dir", "vep_dir") and not merged.get(key):
            merged[key] = value

    # When the operator linked a DNA run, find that run's variant calls so
    # PCGR can report both libraries of the specimen as ONE document.
    if merged.get("paired_run_id"):
        merged.update(paired_dna_inputs(merged["paired_run_id"]))

    rna_env = jobs.find_conda_env(jobs.RNA_ENV_NAME)
    rna_python = jobs.env_python(rna_env) or sys.executable
    argv = jobs.build_rna_argv(rna_python, SCRIPT_FUSION, merged)

    job = manager.submit(SCRIPT_FUSION, argv, patient, meta=merged,
                         pcgr_form=(merged if merged.get("pcgr_refdata_dir")
                                    else None),
                         assay="rna")
    return redirect(url_for("job_view", job_id=job.id))


def paired_dna_inputs(run_id):
    """
    The linked DNA run's variant calls, for a combined PCGR report.

    Returns {} whenever the pairing cannot be honoured -- the run is gone,
    it never finished, it produced no VCF. A combined report is a bonus;
    failing to build one must never cost the RNA report that would
    otherwise have been produced on its own.
    """
    try:
        dna = manager.get(run_id)
        if dna is None or dna.assay == "rna":
            return {}
        output_dir = dna.meta.get("output_dir")
        manifest = (jobs.latest_manifest(output_dir, "dna")
                    if output_dir else None)
        if not manifest:
            return {}
        vcf = manifest.get("cosmic_vcf") or manifest.get("filtered_vcf")
        if not vcf or not os.path.exists(vcf):
            return {}
        return {
            "paired_vcf": vcf,
            # PCGR has to know whether the DNA library had a matched
            # normal: without one it must apply germline filtering, and
            # getting that wrong fills the report with inherited variants.
            "paired_normal": manifest.get("normal") or "",
            "paired_sample": manifest.get("tumour") or "",
            "paired_output_dir": output_dir,
        }
    except Exception:                     # noqa: BLE001 - never block a run
        return {}


@app.route("/new/hybrid")
def new_hybrid_run():
    """The form for a kit whose specimen yields both a DNA and an RNA library."""
    return render_template("new_hybrid_run.html",
                           patient_fields=PATIENT_FIELDS,
                           roots=app.config["ALLOWED_ROOTS"],
                           default_runs=app.config["RUNS_DIR"],
                           defaults=dict(discover_defaults(),
                                         **default_rna_resources()),
                           dna_panels=[o for o in panel_options()
                                       if not o["chemistry"].startswith("rna-")],
                           rna_panels=rna_panel_options(),
                           hybrid_pairs=hybrid_panel_pairs(),
                           tumour_sites=TUMOUR_SITES)


@app.route("/submit/hybrid", methods=["POST"])
def submit_hybrid(validate_only=False):
    """
    Queue both halves of a hybrid panel as one action.

    `validate_only` stops after both halves have been checked and returns
    None, matching submit() and submit_rna(). Without it the worksheet --
    which calls all three handlers uniformly -- crashed on hybrid rows with
    a TypeError, and did so during QUEUEING, after earlier rows had already
    started. Every handler reachable from _submit_via() must accept it.

    WHAT THIS SAVES THE OPERATOR
      Two runs, submitted in the right order, with the RNA half linked back
      to the DNA half so that PCGR produces ONE report for the specimen.
      Done by hand that is two forms, two output directories and a linking
      step that is invisible if you forget it -- and forgetting it costs
      nothing visible at the time and produces two disconnected reports at
      the end.

    ORDER MATTERS AND IS NOT COSMETIC. The DNA run is submitted first so
    that it reaches the front of the serial queue first; by the time the
    RNA run's report step runs, the DNA half has finished and its VCF
    exists. The pairing itself is resolved late -- see build_rna_pcgr_argv()
    -- so even if that ordering is disturbed, the worst case is two
    separate reports rather than a failure.
    """
    form = {k: v.strip() for k, v in request.form.items()}
    patient = {key: form.get(key, "") for key, _ in PATIENT_FIELDS}
    errors = []

    if not patient.get("patient_id") and not patient.get("specimen_id"):
        errors.append("Give at least a patient/MRN or a specimen ID, so both "
                      "runs can be attributed to the same specimen.")

    # --- shared output root, one subdirectory per library ---------------
    base = None
    try:
        base = safe_path(form.get("output_dir"))
    except ValueError as exc:
        errors.append(f"Output directory: {exc}")

    # --- the DNA half ---------------------------------------------------
    dna_form = dict(form)
    dna_form["panel"] = form.get("dna_panel", "")
    dna_form["panel_bed"] = form.get("panel_bed", "")
    dna_form["input_dir"] = form.get("dna_input_dir", "")
    dna_form["tumour_sample"] = form.get("dna_sample", "")
    dna_form["tumour_r1"] = form.get("dna_r1", "")
    dna_form["tumour_r2"] = form.get("dna_r2", "")
    dna_form["auto_discover"] = form.get("dna_auto_discover", "")
    if base:
        dna_form["output_dir"] = os.path.join(base, "dna")

    # --- the RNA half ---------------------------------------------------
    rna_form = dict(form)
    rna_form["panel"] = form.get("rna_panel", "")
    rna_form["input_dir"] = form.get("rna_input_dir", "")
    rna_form["sample"] = form.get("rna_sample", "")
    rna_form["r1"] = form.get("rna_r1", "")
    rna_form["r2"] = form.get("rna_r2", "")
    rna_form["auto_discover"] = form.get("rna_auto_discover", "")
    if base:
        rna_form["output_dir"] = os.path.join(base, "rna")

    for label, half, panel_key, assay in (("DNA", dna_form, "panel", "dna"),
                                          ("RNA", rna_form, "panel", "rna")):
        chosen = half.get(panel_key)
        if chosen and panel_profiles is not None:
            profile = panel_registry().get(chosen)
            if not profile:
                errors.append(f"{label} half: unknown panel {chosen!r}.")
            elif panel_profiles.assay_type(profile) != assay:
                errors.append(
                    f"{label} half was given {chosen!r}, which is a "
                    f"{panel_profiles.assay_type(profile).upper()} profile. "
                    f"Each half needs a profile for its own library type.")

    if errors:
        return render_template(
            "new_hybrid_run.html", errors=errors, form=form, patient=patient,
            patient_fields=PATIENT_FIELDS,
            roots=app.config["ALLOWED_ROOTS"],
            default_runs=app.config["RUNS_DIR"],
            defaults=dict(discover_defaults(), **default_rna_resources()),
            dna_panels=[o for o in panel_options()
                        if not o["chemistry"].startswith("rna-")],
            rna_panels=rna_panel_options(),
            hybrid_pairs=hybrid_panel_pairs(),
            tumour_sites=TUMOUR_SITES), 400

    # The RNA engine needs a real FASTA: Arriba reads the sequence either
    # side of every breakpoint out of it, so it cannot resolve the genome
    # NAME the DNA form accepts. Supply the installed genome's path.
    if not rna_form.get("reference") or \
            not os.path.isabs(rna_form.get("reference", "")):
        installed = default_rna_resources().get("reference")
        if installed:
            rna_form["reference"] = installed
        else:
            errors.append(
                "The RNA half needs the genome FASTA itself, not the name "
                "'hg38' -- the fusion caller reads sequence from it. None "
                "was found in the configured data directory; install the "
                "reference, or start an RNA run on its own and give the "
                "path.")

    # VALIDATE BOTH HALVES BEFORE QUEUEING EITHER.
    #
    # Each half is checked by the SAME handler the single-assay forms use,
    # in validate-only mode. Two reasons, and the second is the one that
    # matters: re-implementing the checks here is how the paths would drift
    # apart, and submitting the DNA half before discovering the RNA half is
    # invalid leaves a run executing behind an error page that says nothing
    # was started.
    if not errors:
        for half, form_values in (("DNA", dna_form), ("RNA", rna_form)):
            handler = submit if half == "DNA" else submit_rna
            problem = _submit_via(handler, form_values, patient,
                                  validate_only=True)
            if problem is not None:
                # The handler rendered its own error page for that half.
                # Return it: it names the fields, which a summary would not.
                return problem

    if errors:
        return render_template(
            "new_hybrid_run.html", errors=errors, form=form, patient=patient,
            patient_fields=PATIENT_FIELDS,
            roots=app.config["ALLOWED_ROOTS"],
            default_runs=app.config["RUNS_DIR"],
            defaults=dict(discover_defaults(), **default_rna_resources()),
            dna_panels=[o for o in panel_options()
                        if not o["chemistry"].startswith("rna-")],
            rna_panels=rna_panel_options(),
            hybrid_pairs=hybrid_panel_pairs(),
            tumour_sites=TUMOUR_SITES), 400

    if validate_only:
        return None                       # both halves checked out

    # Both halves are good. DNA first, so it reaches the front of the
    # serial queue first and its VCF exists by the time the RNA half's
    # report step runs.
    dna_response = _submit_via(submit, dna_form, patient)
    if not hasattr(dna_response, "location"):
        return dna_response

    dna_id = str(dna_response.location).rsplit("/", 1)[-1]
    dna_job = manager.get(dna_id)
    rna_form["paired_run_id"] = dna_id
    if dna_job is not None:
        # The DIRECTORY, not the VCF: the VCF does not exist yet, and will
        # not until the DNA half finishes. build_rna_pcgr_argv() reads the
        # manifest out of it at report time.
        rna_form["paired_output_dir"] = dna_job.meta.get("output_dir", "")

    rna_response = _submit_via(submit_rna, rna_form, patient)
    if not hasattr(rna_response, "location"):
        return rna_response

    rna_id = str(rna_response.location).rsplit("/", 1)[-1]
    # Record the pairing on BOTH runs, so either page can reach the other.
    rna_job = manager.get(rna_id)
    if dna_job is not None and rna_job is not None:
        dna_job.meta["paired_run_id"] = rna_id
        dna_job.meta["hybrid"] = "1"
        rna_job.meta["hybrid"] = "1"
        dna_job._write_state()
        rna_job._write_state()
    return redirect(url_for("job_view", job_id=dna_id))


def _submit_via(handler, form_values, patient, validate_only=False):
    """
    Run one of the single-assay submit handlers against a synthetic form.

    The handlers read request.form, so the hybrid path pushes a request
    context carrying the half it wants submitted. Uglier than calling a
    shared function, and safer than the alternative: it means the hybrid
    route cannot validate a run differently from the form that normally
    submits it.
    """
    payload = {k: v for k, v in form_values.items() if v not in (None, "")}
    payload.update({k: v for k, v in patient.items() if v})
    with app.test_request_context("/", method="POST", data=payload):
        return handler(validate_only=validate_only)


# ---------------------------------------------------------------------------
# The multi-sample worksheet
# ---------------------------------------------------------------------------
# WHY THIS EXISTS
#   Batch mode already queued one run per sample. What it could not do was
#   give each sample its OWN patient. Every run in a batch inherited a
#   single patient record with the sample name appended to its specimen id,
#   which is wrong for the case batches actually are: a sequencing run
#   carrying several different people.
#
#   The worksheet is the lab's own artefact -- one row per specimen, filled
#   in before the run -- rendered as a form. Each row becomes its own job,
#   with its own patient details, its own assay, and its own report.
#
# SEQUENCING IS ALREADY SEQUENTIAL
#   Nothing here needs to orchestrate "finish one sample, then start the
#   next". JobManager is a serial queue by design (BWA-MEM2 wants ~32 GB
#   and STAR ~30 GB; two concurrent runs thrash rather than go faster), and
#   the clinical report runs as part of each job before the next begins. So
#   queueing N rows gives exactly the behaviour asked for: analysis and
#   reporting complete for one sample before the next starts.
#
# NOTHING IS QUEUED UNTIL EVERY ROW VALIDATES
#   The same rule the hybrid form follows, and for a stronger reason here:
#   a worksheet of twelve samples that queued the first four and then
#   rejected the fifth would leave the operator with a half-started batch
#   and no clear way back.

# Per-row fields that are NOT patient details. Everything in PATIENT_FIELDS
# is also per row; see worksheet_row_fields().
WORKSHEET_RUN_FIELDS = ("include", "sample", "assay", "panel",
                        "r1", "r2", "rna_r1", "rna_r2")

# The columns shown inline in the grid. The rest of the patient fields are
# per row too, behind a per-row expander -- a worksheet with eighteen
# visible columns is one nobody fills in correctly.
WORKSHEET_INLINE_PATIENT = ("patient_id", "specimen_id", "diagnosis")


def worksheet_row_fields():
    """Every field a worksheet row carries, run fields then patient ones."""
    return list(WORKSHEET_RUN_FIELDS) + [k for k, _ in PATIENT_FIELDS]


def parse_worksheet_rows(form):
    """
    Pull the rows out of a flat form submission.

    Fields are named row-<n>-<field>. A row is kept when it is ticked AND
    names a sample; anything else is a blank line in the worksheet, which
    is normal -- operators leave spare rows.
    """
    indices = set()
    for key in form:
        if key.startswith("row-"):
            parts = key.split("-", 2)
            if len(parts) == 3 and parts[1].isdigit():
                indices.add(int(parts[1]))

    rows = []
    for index in sorted(indices):
        row = {field: (form.get(f"row-{index}-{field}") or "").strip()
               for field in worksheet_row_fields()}
        row["_index"] = index
        if not row.get("sample"):
            continue                      # an empty line, not an error
        if not row.get("include"):
            continue                      # deliberately excluded
        rows.append(row)
    return rows


def worksheet_row_to_form(row, shared):
    """
    Turn one worksheet row into the form dict a submit handler expects.

    The shared settings (output root, reference, resources, threads) are
    the base; the row overrides what is specific to it. The per-sample
    output directory is derived from the sample name, so two rows cannot
    write into each other's results.
    """
    values = dict(shared)
    assay = (row.get("assay") or "dna").lower()

    base = shared.get("output_dir") or ""
    values["output_dir"] = sample_output_dir(base, row["sample"])

    # Never inherit the shared form's auto-discovery: the row names its
    # files explicitly, and auto-discover would override them and
    # reintroduce exactly the ambiguity the worksheet exists to remove.
    for key in ("auto_discover", "batch_mode", "second_is_matched_normal",
                "manifest"):
        values.pop(key, None)

    if row.get("panel"):
        values["panel"] = row["panel"]

    # The submit handlers require an input directory even when the FASTQs
    # are named explicitly -- it is how they resolve a bare filename. A
    # worksheet row names full paths, so the directory is implied by them;
    # deriving it here is what lets a row whose files live anywhere be
    # submitted by the same handler as a hand-filled form.
    def _dir_of(*paths):
        for candidate in paths:
            if candidate:
                return os.path.dirname(os.path.abspath(candidate))
        return ""

    if assay == "rna":
        values["sample"] = row["sample"]
        values["r1"] = row.get("r1", "")
        values["r2"] = row.get("r2", "")
        values["input_dir"] = _dir_of(row.get("r1"), row.get("r2"))
    else:
        values["tumour_sample"] = row["sample"]
        values["tumour_r1"] = row.get("r1", "")
        values["tumour_r2"] = row.get("r2", "")
        values["input_dir"] = _dir_of(row.get("r1"), row.get("r2"))

    # A hybrid row carries two libraries: R1/R2 are the DNA pair, and the
    # RNA pair has its own columns.
    if assay == "hybrid":
        values["dna_sample"] = row["sample"]
        values["dna_r1"] = row.get("r1", "")
        values["dna_r2"] = row.get("r2", "")
        values["rna_sample"] = row["sample"]
        values["rna_r1"] = row.get("rna_r1", "")
        values["rna_r2"] = row.get("rna_r2", "")
        # Each half resolves its own files, which may sit in different
        # folders -- a hybrid kit's DNA and RNA libraries frequently do.
        values["dna_input_dir"] = _dir_of(row.get("r1"), row.get("r2"))
        values["rna_input_dir"] = _dir_of(row.get("rna_r1"),
                                          row.get("rna_r2"))
        values["input_dir"] = values["dna_input_dir"]
    return values


def worksheet_row_patient(row, shared):
    """
    The patient record for one row.

    Row values win; anything the row leaves blank falls back to the shared
    header (referring clinician and specimen type are usually constant
    across a worksheet, and typing them twelve times invites typos).
    """
    patient = {}
    for key, _label in PATIENT_FIELDS:
        patient[key] = row.get(key) or shared.get(key, "") or ""
    return patient


@app.route("/new/worksheet")
def new_worksheet():
    """The worksheet form, optionally prefilled by scanning a directory."""
    return render_template(
        "worksheet.html",
        patient_fields=PATIENT_FIELDS,
        inline_patient=WORKSHEET_INLINE_PATIENT,
        roots=app.config["ALLOWED_ROOTS"],
        default_runs=app.config["RUNS_DIR"],
        defaults=dict(discover_defaults(), **default_rna_resources()),
        dna_panels=[o for o in panel_options()
                    if not o["chemistry"].startswith("rna-")],
        rna_panels=rna_panel_options(),
        tumour_sites=TUMOUR_SITES,
        rows=[], scanned=None, errors=None, form=None)


@app.route("/worksheet/scan", methods=["POST"])
def worksheet_scan():
    """
    Prefill the worksheet from a directory of FASTQs.

    The single biggest saving this page offers: point at the run folder and
    get one row per detected sample, with the file paths already correct.
    The operator then fills in who each sample belongs to -- which is the
    part only they can know.

    Uses the pipeline's own pairing rules via discover_pairs(), so the
    worksheet cannot disagree with what the engine would find.
    """
    form = {k: v.strip() for k, v in request.form.items()}
    errors = []
    rows = []
    scanned = None

    try:
        input_dir = safe_path(form.get("scan_dir"), must_exist=True)
        pairs, warnings = discover_pairs(input_dir,
                                         bool(form.get("scan_recursive")))
        for warning in warnings:
            errors.append(f"input scan: {warning}")
        scanned = input_dir
        default_assay = form.get("scan_assay") or "dna"
        for pair in pairs:
            rows.append({
                "include": "on",
                "sample": pair["sample"],
                "assay": default_assay,
                "r1": pair.get("r1", ""),
                "r2": pair.get("r2", ""),
            })
        if not pairs:
            errors.append(f"No FASTQ pairs found in {input_dir}.")
    except ValueError as exc:
        errors.append(f"Directory to scan: {exc}")

    return render_template(
        "worksheet.html",
        patient_fields=PATIENT_FIELDS,
        inline_patient=WORKSHEET_INLINE_PATIENT,
        roots=app.config["ALLOWED_ROOTS"],
        default_runs=app.config["RUNS_DIR"],
        defaults=dict(discover_defaults(), **default_rna_resources()),
        dna_panels=[o for o in panel_options()
                    if not o["chemistry"].startswith("rna-")],
        rna_panels=rna_panel_options(),
        tumour_sites=TUMOUR_SITES,
        rows=rows, scanned=scanned,
        errors=errors or None, form=form)


@app.route("/submit/worksheet", methods=["POST"])
def submit_worksheet():
    """
    Validate every row, then queue them in worksheet order.

    Two properties matter and both are deliberate:

    NOTHING STARTS UNTIL EVERYTHING VALIDATES. Each row is checked by the
    SAME handler that would submit it on its own, in validate-only mode.
    Re-implementing the checks here is how the worksheet path would drift
    away from the single-sample path, and the worksheet is the one nobody
    exercises by hand.

    ROWS RUN IN ORDER, ONE AT A TIME. That is not orchestrated here -- the
    job queue is serial and each job produces its report before the next
    begins. Queueing in worksheet order is all that is needed for "finish
    one sample, then the next".
    """
    form = {k: v.strip() for k, v in request.form.items()}
    rows = parse_worksheet_rows(form)

    # Shared header: everything not per row. The row fields are stripped so
    # a stray top-level value cannot leak into every run.
    shared = {k: v for k, v in form.items() if not k.startswith("row-")}
    shared.pop("scan_dir", None)
    shared.pop("scan_recursive", None)
    shared.pop("scan_assay", None)

    errors = []
    if not rows:
        errors.append(
            "No rows to run. Tick the samples you want and give each one a "
            "name, or scan a directory to fill the worksheet in.")

    seen = {}
    for row in rows:
        name = row["sample"]
        if name in seen:
            errors.append(
                f"Sample name {name!r} appears more than once (rows "
                f"{seen[name] + 1} and {row['_index'] + 1}). Each row writes "
                f"into a directory named after its sample, so duplicates "
                f"would overwrite each other's results.")
        seen[name] = row["_index"]

        if row.get("assay") not in ("dna", "rna", "hybrid"):
            errors.append(f"Row {row['_index'] + 1} ({name}): choose an "
                          f"assay.")
        if not row.get("patient_id") and not row.get("specimen_id") \
                and not shared.get("patient_id") \
                and not shared.get("specimen_id"):
            errors.append(
                f"Row {row['_index'] + 1} ({name}): give a patient/MRN or a "
                f"specimen ID, so the report can be attributed to "
                f"something.")

    # --- validate every row before starting any of them ----------------
    if not errors:
        for row in rows:
            assay = row["assay"]
            values = worksheet_row_to_form(row, shared)
            patient = worksheet_row_patient(row, shared)
            handler = {"dna": submit, "rna": submit_rna,
                       "hybrid": submit_hybrid}[assay]
            problem = _submit_via(handler, values, patient,
                                  validate_only=True)
            if problem is not None:
                # The handler rendered its own error page naming the
                # fields. Wrap it so the operator knows WHICH row failed --
                # on a twelve-row worksheet that is the whole message.
                errors.append(
                    f"Row {row['_index'] + 1} ({row['sample']}, "
                    f"{assay.upper()}) did not validate. Its details are "
                    f"below; nothing has been queued.")
                return render_worksheet(rows, form, errors,
                                        detail=problem)

    if errors:
        return render_worksheet(rows, form, errors), 400

    # --- queue, in worksheet order -------------------------------------
    queued = []
    first = None
    for row in rows:
        assay = row["assay"]
        values = worksheet_row_to_form(row, shared)
        patient = worksheet_row_patient(row, shared)
        values["worksheet_row"] = str(row["_index"] + 1)
        values["worksheet_of"] = str(len(rows))

        handler = {"dna": submit, "rna": submit_rna,
                   "hybrid": submit_hybrid}[assay]
        try:
            response = _submit_via(handler, values, patient)
        except Exception as exc:          # noqa: BLE001 - surfaced below
            # A crash mid-batch is the worst outcome here: some rows are
            # already running and the operator is looking at a stack
            # trace. Turn it into a page that says exactly how far it got.
            errors.append(
                f"Row {row['_index'] + 1} ({row['sample']}) raised "
                f"{type(exc).__name__}: {exc}. {len(queued)} earlier row(s) "
                f"ARE already running; the remaining rows were not "
                f"started.")
            return render_worksheet(rows, form, errors), 500
        if not hasattr(response, "location"):
            # Should not happen -- every row validated a moment ago -- but
            # a run that fails to queue must not be silently missing from a
            # batch the operator believes is complete.
            errors.append(
                f"Row {row['_index'] + 1} ({row['sample']}) validated but "
                f"could not be queued. {len(queued)} earlier row(s) ARE "
                f"running; the rest were not started.")
            return render_worksheet(rows, form, errors), 500
        job_id = str(response.location).rsplit("/", 1)[-1]
        queued.append(job_id)
        first = first or job_id

    return redirect(url_for("job_view", job_id=first))


def render_worksheet(rows, form, errors, detail=None, status=None):
    """Re-render the worksheet with what was typed and what went wrong."""
    return render_template(
        "worksheet.html",
        patient_fields=PATIENT_FIELDS,
        inline_patient=WORKSHEET_INLINE_PATIENT,
        roots=app.config["ALLOWED_ROOTS"],
        default_runs=app.config["RUNS_DIR"],
        defaults=dict(discover_defaults(), **default_rna_resources()),
        dna_panels=[o for o in panel_options()
                    if not o["chemistry"].startswith("rna-")],
        rna_panels=rna_panel_options(),
        tumour_sites=TUMOUR_SITES,
        rows=rows, scanned=None, errors=errors, form=form,
        detail=detail)


@app.route("/panels")
def panels_page():
    """
    Reference page: every profile, what it sets, and what it warns about.

    The dropdown can only show a line of summary. The caveats -- that this
    pipeline does not build UMI consensus reads, that duplicate marking is
    off for amplicon chemistry, that TMB under 1 Mb is not reportable --
    are the part somebody needs to have read before they use the result,
    so they get a page rather than a tooltip.
    """
    return render_template("panels.html",
                           panels=panel_options(),
                           details=json.loads(panel_settings_json()),
                           available=panel_profiles is not None)


# ---------------------------------------------------------------------------
# Tumour site
# ---------------------------------------------------------------------------
# PCGR tiers actionability against the tumour site: the same variant can be
# predictive in one tissue and merely oncogenic in another, so naming the site
# changes which evidence PCGR considers on-label. Code 0 is "Any", which is
# the honest default -- it keeps the site OPEN and reports every tier without
# claiming a tissue nobody specified.
#
# Copied from PCGR 2.3.2 (pcgr.pcgr_vars.tsites). It lives in the pcgr conda
# environment, which this app cannot import from, so it is duplicated here.
# Regenerate after a PCGR upgrade with:
#   conda run -n pcgr python -c \
#     "from pcgr import pcgr_vars; print(pcgr_vars.tsites)"
TUMOUR_SITES = [
    (0, "Any / not specified"),
    (1, "Adrenal Gland"),
    (2, "Ampulla of Vater"),
    (3, "Biliary Tract"),
    (4, "Bladder/Urinary Tract"),
    (5, "Bone"),
    (6, "Breast"),
    (7, "Cervix"),
    (8, "CNS/Brain"),
    (9, "Colon/Rectum"),
    (10, "Esophagus/Stomach"),
    (11, "Eye"),
    (12, "Head and Neck"),
    (13, "Kidney"),
    (14, "Liver"),
    (15, "Lung"),
    (16, "Lymphoid"),
    (17, "Myeloid"),
    (18, "Ovary/Fallopian Tube"),
    (19, "Pancreas"),
    (20, "Peripheral Nervous System"),
    (21, "Peritoneum"),
    (22, "Pleura"),
    (23, "Prostate"),
    (24, "Skin"),
    (25, "Soft Tissue"),
    (26, "Testis"),
    (27, "Thymus"),
    (28, "Thyroid"),
    (29, "Uterus"),
    (30, "Vulva/Vagina"),
]

TUMOUR_SITE_CODES = {code for code, _label in TUMOUR_SITES}


def tumour_site_label(code):
    """Human-readable site for a code, for the run record."""
    try:
        code = int(code)
    except (TypeError, ValueError):
        return None
    return dict(TUMOUR_SITES).get(code)


# ---------------------------------------------------------------------------
# Resource discovery
# ---------------------------------------------------------------------------
# A run submitted with these fields blank still succeeds -- the pipeline skips
# each step cleanly -- and produces an annotated VCF and no clinical report.
# That is rarely what anyone wants, and it is not obvious from the form, so
# the standard install locations are pre-filled when they exist. Every value
# stays editable and clearable; this changes the default, not the rules.
_DATA = os.path.expanduser("~/data")


def default_resources(data=None):
    """
    The installer's layout, under whichever data directory is in use.

    Derived rather than fixed at import, so --data-dir moves the suggested
    paths, the database check and the upgrade actions together instead of
    only the last two.
    """
    data = data or app.config["DATA_DIR"]
    return {
        "dbsnp":
            f"{data}/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        "germline_resource":
            f"{data}/resources/hg38/af-only-gnomad.hg38.vcf.gz",
        "panel_of_normals": f"{data}/resources/hg38/1000g_pon.hg38.vcf.gz",
        "contamination_resource":
            f"{data}/resources/hg38/small_exac_common_3.hg38.vcf.gz",
        "cosmic":
            f"{data}/resources/hg38/Cosmic_GenomeScreensMutant_v103_GRCh38"
            f".chr.vcf.gz",
        "vep_dir": f"{data}/vep_cache",
        # The genome install_pipeline.py downloads and indexes. Offering the
        # installed FASTA rather than the name "hg38" is what stops the
        # pipeline downloading and re-indexing a private copy for every new
        # output directory; --reference-dir covers the case where it must
        # download.
        "reference": f"{data}/references/hg38/Homo_sapiens_assembly38.fasta",
        "reference_dir": f"{data}/references",
        # Tumour-only MSI. PCGR omits MSI on a panel, so this is the only
        # route to an MSI answer for a targeted assay.
        "msi_models": f"{data}/msisensor2/models_hg38",
    }


def default_known_indels(data=None):
    """--known-indels takes several files, so it is kept apart from the map."""
    data = data or app.config["DATA_DIR"]
    return [
        f"{data}/resources/hg38/Mills_and_1000G_gold_standard.indels.hg38"
        f".vcf.gz",
        f"{data}/resources/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz",
    ]


def discover_defaults():
    """
    Standard resource paths that actually exist on this machine.

    Missing entries are omitted rather than offered, so the form never
    proposes a path that would fail validation. The PCGR bundle is found by
    globbing, because its directory is named for the bundle release date and
    changes with every refresh.
    """
    data = app.config["DATA_DIR"]
    found = {k: v for k, v in default_resources(data).items()
             if os.path.exists(v)}

    bundles = sorted(glob.glob(os.path.join(data, "pcgr", "[0-9]" * 8)))
    if bundles:
        found["pcgr_refdata_dir"] = bundles[-1]

    # Clinical reports are the point of running this at all, so the estimate
    # boxes start ticked. MSI is left off: PCGR honours it only for WGS/WES
    # tumour-normal runs and ignoring that just puts a misleading blank in
    # the report.
    indels = [p for p in default_known_indels(data) if os.path.exists(p)]
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
    #
    # A toggle the PANEL decided is left alone either way. A hotspot panel's
    # profile turns TMB off because 22 kb is not a denominator anybody can
    # divide by; letting the generic "clinical reports are the point" default
    # turn it straight back on would produce exactly the confidently wrong
    # number the profile exists to prevent.
    decided = set(str(form.get("_panel_toggles", "")).split())
    if not form.get("form_rendered"):
        for key, value in ANALYSIS_TOGGLES.items():
            if value and not form.get(key) and key not in decided:
                form[key] = value
                applied.append(key)
    return applied


@app.route("/submit", methods=["POST"])
def submit(validate_only=False):
    """
    Validate and queue a DNA run.

    `validate_only` runs every check and stops before anything is queued,
    returning None when the submission is good. The hybrid route uses it to
    check BOTH halves before starting either -- otherwise a bad RNA half
    leaves a DNA run already executing behind an error page that says
    nothing was started.
    """
    form = {k: v.strip() for k, v in request.form.items()}
    patient = {key: form.get(key, "") for key, _ in PATIENT_FIELDS}

    # The panel goes first: it decides the target BED and the analysis
    # toggles, and apply_resource_defaults() below must see those decisions
    # rather than overwrite them. Both run before validation, so anything
    # they fill in is path-checked exactly like a value the user typed.
    panel_applied = apply_panel_to_form(form)
    auto_applied = apply_resource_defaults(form)

    errors = []
    # A panel id that does not resolve is refused rather than ignored: a
    # run configured from a profile that silently did not exist is a run
    # configured from nothing, and it would look identical on the page.
    if form.get("panel"):
        if panel_profiles is None:
            errors.append(
                "A panel was chosen but panel_profiles.py is not importable "
                "from the pipeline directory, so none of its settings were "
                "applied. Clear the panel, or restore the file.")
        elif form["panel"] not in panel_registry():
            errors.append(
                f"Unknown panel profile {form['panel']!r}. Pick one from the "
                f"list, or add your own to "
                f"~/.config/cancer_pipeline/panels.")
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
        # Batch and "these two are a pair" both describe what to do with the
        # samples auto-discover finds, so neither means anything without it.
        # Ticked alone they used to fall through to the explicit-sample
        # checks below, which complained that "tumour sample is required" --
        # three errors naming neither box.
        needs_discovery = [
            label for key, label in
            (("batch_mode", "Batch"),
             ("second_is_matched_normal",
              "'The two samples are a tumour/normal pair from one person'"))
            if form.get(key)]
        if needs_discovery and not form.get("auto_discover"):
            errors.append(
                f"{' and '.join(needs_discovery)} decides what to do with the "
                f"samples found in the input directory, so it needs "
                f"'Auto-discover FASTQ pairs' ticked as well. Tick that, or "
                f"untick it and name the tumour (and normal) files "
                f"explicitly.")
        elif not form.get("auto_discover"):
            for field in ("tumour_sample", "tumour_r1", "tumour_r2"):
                if not form.get(field):
                    errors.append(f"{field.replace('_', ' ')} is required "
                                  f"unless auto-discover or a manifest is used")
            # The FASTQ fields are paths like any other, and were the one
            # set that skipped safe_path(): a value outside --allow-root was
            # forwarded to the pipeline, and a typo was only discovered when
            # the run failed minutes later. Both are caught here now.
            for field in ("tumour_r1", "tumour_r2",
                          "normal_r1", "normal_r2"):
                if form.get(field):
                    try:
                        resolved[field] = safe_path(form[field],
                                                    must_exist=True)
                    except ValueError as exc:
                        errors.append(f"{field.replace('_', ' ')}: {exc}")
            if form.get("normal_sample") and not (form.get("normal_r1")
                                                  and form.get("normal_r2")):
                errors.append("A normal sample needs both normal R1 and R2")

    # The target intervals and the coverage BED are usually copies of the
    # panel BED, put there by apply_panel_to_form(). Checking all three
    # separately reports one bad path as three errors naming two fields the
    # operator never filled in, so the copies take the panel BED's own
    # result instead of being re-checked. Note the loop order: panel_bed is
    # resolved before the fields that derive from it.
    derived_from_bed = {field for field in ("intervals", "coverage_bed")
                        if field in panel_applied
                        and form.get(field) == form.get("panel_bed")}

    for optional in ("cosmic", "dbsnp", "germline_resource",
                     "panel_of_normals", "contamination_resource",
                     "msi_models", "panel_bed", "intervals", "coverage_bed",
                     "pcgr_refdata_dir", "vep_dir", "reference_dir"):
        if optional in derived_from_bed:
            if "panel_bed" in resolved:
                resolved[optional] = resolved["panel_bed"]
            continue
        if form.get(optional):
            try:
                resolved[optional] = safe_path(form[optional], must_exist=True)
            except ValueError as exc:
                errors.append(f"{optional.replace('_', ' ')}: {exc}")

    # A bundle without a cache is a run that completes and then cannot
    # report. PCGR requires --vep_dir whenever it is given an input VCF, and
    # it says so only once it starts -- which is after the pipeline has spent
    # its hours. Catch it here, where it costs nothing.
    if form.get("pcgr_refdata_dir") and not form.get("vep_dir"):
        errors.append(
            "A PCGR reference bundle was given without a VEP cache. PCGR "
            "annotates with Ensembl VEP and refuses to start without one, so "
            "the run would finish and then fail to produce a report. Fill in "
            "the VEP cache (the installer puts it in ~/data/vep_cache), or "
            "clear the bundle to run without a clinical report.")

    # Validate the site code rather than forwarding whatever was posted:
    # PCGR rejects an out-of-range value, and it does so after the pipeline
    # has already run.
    if form.get("pcgr_tumour_site"):
        try:
            if int(form["pcgr_tumour_site"]) not in TUMOUR_SITE_CODES:
                raise ValueError
        except (TypeError, ValueError):
            errors.append(
                f"tumour site must be one of PCGR's codes 0-"
                f"{max(TUMOUR_SITE_CODES)}; 0 leaves it unspecified")

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
                                             **ANALYSIS_TOGGLES),
                               panels=panel_options(),
                               panel_settings=panel_settings_json(),
                               tumour_sites=TUMOUR_SITES), 400

    if validate_only:
        return None                       # everything checked out

    merged = dict(form)
    merged.update(resolved)
    merged["auto_applied_resources"] = " ".join(auto_applied)
    # Same reasoning as auto_applied_resources: a run record that shows a
    # dozen settings without saying where they came from leaves a later
    # reader to guess which were chosen and which were inherited.
    merged["panel_applied_fields"] = " ".join(panel_applied)
    if form.get("panel") and panel_profiles is not None:
        profile = panel_registry().get(form["panel"])
        if profile:
            merged["panel_label"] = f"{profile['name']} ({profile['id']})"
            merged["panel_chemistry"] = profile["chemistry"]
            merged["panel_notes"] = profile.get("notes", [])
    # Store the label too: a run record reading "Tumour site: 15" tells a
    # later reader nothing without the codebook.
    merged["pcgr_tumour_site_label"] = tumour_site_label(
        merged.get("pcgr_tumour_site"))
    # The pipeline runs in the cancer_pipeline environment, not in whichever
    # interpreter is serving this app: launching it with sys.executable is
    # why a webapp started outside that env failed with
    # "Required tools not found: bwa-mem2, samtools, gatk, fastp".
    pipeline_env = jobs.find_conda_env(jobs.PIPELINE_ENV_NAME)
    pipeline_python = jobs.env_python(pipeline_env) or sys.executable
    argv = build_pipeline_argv(pipeline_python, SCRIPT_VARIANT, merged)

    # PCGR runs as a follow-up in its own conda environment -- it cannot run
    # inside the pipeline's, so passing --pcgr-refdata-dir alone would leave
    # step 12 reporting itself as skipped and no report on disk.
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


def load_db_status():
    """
    The raw status check_db_updates.py last wrote, or None.

    Shared by the front-page notice and the databases page so both speak
    from the same file rather than checking anything themselves -- no page
    render ever touches the network.
    """
    try:
        sys.path.insert(0, os.path.dirname(SCRIPT_VARIANT))
        from check_db_updates import load_status
    except Exception:                      # noqa: BLE001 - optional sibling
        return None
    try:
        return load_status()
    except Exception:                      # noqa: BLE001 - never fatal
        return None


def database_notice():
    """
    What the last database check found, or None.

    Returns None when no check has ever run, so a machine that never set up
    the cron entry shows nothing rather than a permanent scolding. A stale
    result is surfaced as stale rather than quietly presented as current --
    "checked three months ago, all well" is not the same claim as "all well".
    """
    status = load_db_status()
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


def rna_run_warnings(meta):
    """
    What an RNA run gave up, and what its reader must not assume.

    The DNA equivalent is about resources that were left blank. Here the
    dominant risk is different and worse: an RNA run that produces nothing
    looks exactly like an RNA run that found nothing, and every artefact
    downstream -- the table, the PDF, the summary line -- is identical in
    the two cases.
    """
    notes = []

    if meta.get("panel_label"):
        for note in (meta.get("panel_notes") or []):
            notes.append(f"{meta['panel_label']}: {note}")
    else:
        notes.append(
            "No panel profile was named, so the library-size and "
            "mapping-rate floors that decide whether a negative is "
            "interpretable were left at their generic defaults. Those "
            "differ by an order of magnitude between an anchored-PCR panel "
            "and a whole transcriptome.")

    if not meta.get("rna_quality"):
        notes.append(
            "No RNA quality figure (DV200/RIN) was recorded. It is an "
            "instrument measurement of the extracted RNA, it cannot be "
            "recovered from the sequencing data, and it is the single best "
            "predictor of whether fusion detection could have worked at "
            "all. A negative result without it is hard to defend.")

    if not meta.get("arriba_resources"):
        notes.append(
            "No explicit Arriba reference directory was given. The files "
            "ship inside the conda environment and are normally found "
            "automatically -- but if they were not, the run had no "
            "blacklist, and recurrent read-through artefacts will appear in "
            "the table as high-confidence fusions. The run log says which "
            "files were located.")

    if not meta.get("paired_run_id"):
        notes.append(
            "This RNA run is not linked to a DNA run. Fusions and SNVs from "
            "one specimen are one result reported together; linking them "
            "makes that explicit in the record rather than leaving it to "
            "whoever reads the two reports.")

    if meta.get("skip_steps"):
        skipped = str(meta["skip_steps"]).split()
        if "rnaqc" in skipped:
            notes.append(
                "The RNA library QC step was skipped. Nothing in this run "
                "now says whether an empty fusion table is a negative "
                "result or an unusable library.")
        if "fusions" in skipped:
            notes.append("Fusion calling itself was skipped.")

    notes.append(
        "Fusion calls from RNA are evidence of a transcript, not proof of a "
        "genomic rearrangement. Callers disagree substantially on real "
        "data; orthogonal confirmation is expected before clinical "
        "reporting.")
    return notes


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
    # An RNA run shares almost none of the DNA warnings below -- it has no
    # germline resource, no BQSR, no panel of normals, no PCGR -- so it gets
    # its own list rather than a DNA list with most items suppressed.
    if meta.get("assay") == "rna":
        return rna_run_warnings(meta)

    notes = []
    tumour_only = not meta.get("normal_sample")

    # What assay these numbers belong to, and what that assay's profile
    # wanted the reader to know. The caveats are carried into the run
    # record on purpose: whoever reads this page months later will not have
    # the panels page open beside it.
    if meta.get("panel_label"):
        for note in (meta.get("panel_notes") or []):
            notes.append(f"{meta['panel_label']}: {note}")
        # The one combination that destroys a run silently. Amplicon reads
        # share primer coordinates, so duplicate marking flags nearly the
        # whole library and the caller then works from what little
        # survived. The profile skips the step; this catches an operator
        # who cleared the box afterwards.
        if meta.get("panel_chemistry") == "amplicon" and \
                "dedup" not in str(meta.get("skip_steps", "")).split():
            notes.append(
                "This is an amplicon panel and duplicate marking was NOT "
                "skipped. Amplicon reads all start and end at primer "
                "coordinates, so MarkDuplicates flags almost the entire "
                "library and the variant caller sees a small fraction of "
                "the real depth. Add 'dedup' to the skipped steps and "
                "re-run.")
    elif meta.get("intervals") or meta.get("coverage_bed"):
        notes.append(
            "No panel profile was named, so the chemistry-specific settings "
            "-- interval padding, the VAF and depth floors, whether "
            "duplicate marking is meaningful for this library -- were left "
            "at their general-purpose defaults. Choosing the assay on the "
            "new-run form sets them together.")

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
            "outnumber the real calls by two orders of magnitude. The panel "
            "BED field fills this and the coverage BED together.")
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
    if not meta.get("coverage_bed"):
        notes.append(
            "No target BED, so the run makes no coverage statement. Nothing "
            "distinguishes a region that was sequenced and is wild type from "
            "one the sequencing never reached: both are simply absent from "
            "the VCF and absent from the report.")
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


def coverage_links(job):
    """[(index, sample name)] for each coverage report this run has written."""
    return [(i, os.path.basename(p)[:-len(".coverage.html")])
            for i, p in enumerate(
                jobs.find_coverage_reports(job.meta.get("output_dir")))]


def render_job(job, notice=None, notice_kind="ok", report_error=None,
               status=200):
    """
    Render the job page.

    Shared by the view and by actions that report a result, so an
    action can show its outcome without a session -- this app has no
    secret key and does not need one for a single-analyst tool.

    EVERY path that shows this page goes through here. The PDF-failure
    path used to call render_template() directly with five of the
    template's variables, which rendered a page missing its output links,
    its coverage reports, its warnings and its progress -- all silently,
    because an undefined name in Jinja is simply falsy. The one moment an
    operator most needs the rest of the page is the moment a report failed
    to build.
    """
    entries, total = jobs.survey_intermediates(
        job.meta.get("output_dir"))
    items = run_artefacts(job)
    return render_template(
        "job.html", job=job.snapshot(),
        patient=job.patient,
        patient_fields=PATIENT_FIELDS,
        char_warnings=patient_warnings(job.patient),
        run_notes=run_warnings(job.meta),
        cleanup={"entries": entries, "total": total},
        human_bytes=jobs.human_bytes,
        notice=notice, notice_kind=notice_kind,
        report_error=report_error,
        batch_of=job.meta.get("batch_of"),
        auto_applied=[k for k in
                      (job.meta.get("auto_applied_resources") or
                       "").split() if k],
        coverage=coverage_links(job),
        coverage_bed=job.meta.get("coverage_bed"),
        # The RNA branch's two reports, listed the same way and for the
        # same reason: both are readable while the run is still going.
        rna_fusion_reports=[
            (i, os.path.basename(p)[:-len(".fusions.html")])
            for i, p in enumerate(jobs.find_rna_reports(
                job.meta.get("output_dir"), "fusion"))],
        rna_qc_reports=[
            (i, os.path.basename(p)[:-len(".rna_qc.html")])
            for i, p in enumerate(jobs.find_rna_reports(
                job.meta.get("output_dir"), "qc"))],
        paired_run_id=job.meta.get("paired_run_id"),
        # Every output of the run, grouped and described in plain English.
        # See webapp/artefacts.py: the group ORDER is a clinical argument,
        # not a filing convention.
        artefact_groups=artefacts.grouped(items),
        artefact_absent=artefacts.absent_notes(job.snapshot(), items,
                                               job.meta),
        human_size=jobs.human_bytes,
        # Which assay this run was configured for, and which fields the
        # profile filled in. Shown for the same reason the auto-applied
        # resources are: a setting that arrived from somewhere the reader
        # cannot see is what makes two runs disagree inexplicably.
        panel_label=job.meta.get("panel_label"),
        panel_chemistry=job.meta.get("panel_chemistry"),
        panel_applied=[k for k in
                       (job.meta.get("panel_applied_fields") or
                        "").split() if k])


@app.route("/job/<job_id>")
def job_view(job_id):
    return render_job(manager.get(job_id) or abort(404))


@app.route("/api/job/<job_id>")
def job_api(job_id):
    """Polled by the job page for live status and log tail."""
    job = manager.get(job_id) or abort(404)
    # The coverage reports are written mid-run, before PCGR; the page shows
    # the link the moment the file lands rather than at the end of the run.
    return jsonify({"job": job.snapshot(), "log": job.tail(),
                    "coverage": [{"index": i, "sample": name}
                                 for i, name in coverage_links(job)]})


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
    # The buttons are hidden during the report phase, but the URL is not, and
    # a PDF built now would state "PCGR did not run for this analysis" about
    # a PCGR that is running as it is being written.
    if job.status == "reporting":
        abort(409, "PCGR is still running for this job; the PDF would record "
                   "that no clinical report exists. Wait for the run to "
                   "finish.")

    pdf_path = os.path.join(job.run_dir, f"report_{job_id}.pdf")
    try:
        if job.assay == "rna":
            # A different document, not a flag on the same one: an RNA run
            # has no variant table, no target coverage and no PCGR, and
            # half a report saying "not applicable" teaches a reader to
            # skim. See webapp/rna_report.py.
            result = build_rna_report(job.snapshot(), job.patient,
                                      output_dir, pdf_path)
        else:
            result = build_report(job.patient, job.snapshot(), output_dir,
                                  pdf_path)
    except Exception as exc:  # noqa: BLE001 - surfaced to the user
        return render_job(job,
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
    """
    Show what the run produced, and where the interpretation lives.

    ASSAY-AWARE, because it was not and the result was actively
    misleading: on an RNA run it read a DNA manifest that does not exist
    (so every field was blank), then reported "no coverage report -- this
    run was submitted without a target BED", which on an RNA run is not a
    missing input but a concept that does not apply. A page that invents a
    shortcoming teaches an operator to distrust the ones that are real.
    """
    job = manager.get(job_id) or abort(404)
    output_dir = job.meta.get("output_dir") or ""
    assay = job.assay
    manifest = {}
    if output_dir:
        manifest = (jobs.latest_manifest(output_dir, "rna") if assay == "rna"
                    else load_pipeline_manifest(output_dir)) or {}
    pcgr = find_pcgr_outputs(output_dir) if output_dir else {"html": [],
                                                             "tsv": [],
                                                             "dir": ""}
    pdf_path = os.path.join(job.run_dir, f"report_{job_id}.pdf")

    # Links are built from the run's artefact list so they carry the same
    # STABLE ids the job page uses. Addressing these by position in a glob
    # was the bug: the set grows while the run is going, and a link then
    # serves a different file than the one it is labelled with.
    def links(paths, suffix):
        return [(artefacts.artefact_id(p, job.run_dir, output_dir),
                 os.path.basename(p)[:-len(suffix)]) for p in paths]

    # For an RNA run the question "was the sequencing good enough?" is
    # answered by the library-adequacy report, not by target coverage.
    quality = []
    if assay == "rna":
        quality = links(jobs.find_rna_reports(output_dir, "qc"),
                        ".rna_qc.html")
    fusions = links(jobs.find_rna_reports(output_dir, "fusion"),
                    ".fusions.html")

    return render_template("results.html", job=job.snapshot(),
                           assay=assay,
                           manifest=manifest, pcgr=pcgr,
                           coverage=links(
                               jobs.find_coverage_reports(output_dir),
                               ".coverage.html"),
                           coverage_bed=job.meta.get("coverage_bed"),
                           rna_quality=quality, rna_fusions=fusions,
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




def run_artefacts(job):
    """Every file this run produced, described and ordered for reading."""
    return artefacts.collect(job.snapshot(), job.run_dir,
                             job.meta.get("output_dir"))


@app.route("/job/<job_id>/file/<artefact>")
def job_file(job_id, artefact):
    """
    Serve one of a run's output files.

    Addressed by INDEX into the discovered list, never by a path from the
    request -- the same rule the PCGR, coverage and RNA routes follow. The
    index is re-resolved and re-checked for containment on every request,
    so a job record edited on disk cannot turn this into a file-read
    primitive.

    HTML and PDF open in the browser; everything else downloads, because a
    browser rendering 40 MB of VCF inline helps nobody.
    """
    job = manager.get(job_id) or abort(404)
    items = run_artefacts(job)
    path, item = artefacts.resolve(items, artefact, job.run_dir,
                                   job.meta.get("output_dir"))
    if not path or not os.path.exists(path):
        abort(404, "no such output file for this run")
    if item["action"] == artefacts.VIEW:
        return send_file(path)
    return send_file(path, as_attachment=True,
                     download_name=os.path.basename(path))




# ---------------------------------------------------------------------------
# Reference databases: check, and upgrade
# ---------------------------------------------------------------------------
# The check is safe -- it reads remote metadata and writes one small status
# file. An upgrade is not: it replaces the files runs read, downloads tens of
# gigabytes, and shifts annotations for every report produced afterwards. So
# an upgrade is refused while any run is queued or in progress, and asks for
# a typed confirmation rather than being one click away from a stray tap.


def db_rows():
    """
    Each checked database, with what can be done about it here.

    Built from the last check's results so the page shows the same facts as
    the front-page notice, plus the action for the rows that have one.
    """
    notice = load_db_status()
    rows = []
    for r in (notice or {}).get("results", []):
        spec = maintenance.UPDATERS.get(r["name"])
        rows.append({
            "name": r["name"],
            "state": r.get("state", "unknown"),
            "installed": r.get("installed"),
            "latest": r.get("latest"),
            "note": r.get("note"),
            "url": r.get("url"),
            "action": spec["label"] if spec else None,
            "warning": spec["warning"] if spec else None,
            "needs_file": bool(spec and spec.get("needs_file")),
            # Prefilled with the version the check found, because the
            # installer would otherwise reinstall the release it is pinned
            # to -- the one just reported as out of date.
            "fields": [dict(f, value=(r.get("latest") or ""
                                      if f.get("from_latest") else ""))
                       for f in (spec or {}).get("fields", [])],
        })
    return rows


@app.route("/databases")
def databases():
    return render_template(
        "databases.html",
        rows=db_rows(),
        notice=request.args.get("notice"),
        db_notice=database_notice(),
        db_task=db_tasks.snapshot() if db_tasks else None,
        run_busy=manager.busy() if manager else False,
        roots=app.config["ALLOWED_ROOTS"])


@app.route("/databases/check", methods=["POST"])
def databases_check():
    if db_tasks is None:
        abort(503, "database tasks are not configured")
    _task, error = db_tasks.start_check(reason="requested")
    return redirect(url_for("databases", notice=error or "Checking now."))


@app.route("/databases/update", methods=["POST"])
def databases_update():
    if db_tasks is None:
        abort(503, "database tasks are not configured")
    name = request.form.get("name", "")

    # An upgrade mid-run would swap the databases underneath it, and would
    # compete for the same disk and network the run needs.
    if manager.busy():
        return redirect(url_for(
            "databases",
            notice="A run is queued or in progress. Upgrading now would "
                   "replace the databases it is reading; wait for it to "
                   "finish."))

    # Typed, not a checkbox: this is tens of gigabytes and hours, and it
    # changes what every later report says.
    if request.form.get("confirm", "").strip().upper() != "UPDATE":
        return redirect(url_for(
            "databases",
            notice="Type UPDATE in the confirmation box to start an upgrade."))

    cosmic = None
    if request.form.get("cosmic"):
        try:
            cosmic = safe_path(request.form["cosmic"], must_exist=True)
        except ValueError as exc:
            return redirect(url_for("databases", notice=f"COSMIC VCF: {exc}"))

    # Version targets go on a command line and into a directory name, so
    # they are constrained here rather than trusted from the form.
    extra = []
    for field in maintenance.UPDATERS.get(name, {}).get("fields", []):
        value = (request.form.get(field["name"]) or "").strip()
        if not value:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,31}", value):
            return redirect(url_for(
                "databases",
                notice=f"{field['label']}: {value!r} is not a version "
                       f"(letters, digits, dot, dash, underscore)."))
        extra.append((field["flag"], value))

    _task, error = db_tasks.start_update(name, cosmic=cosmic, extra=extra)
    return redirect(url_for(
        "databases",
        notice=error or f"Started: {name}. It keeps running if you close "
                        f"this page."))


@app.route("/api/databases")
def databases_api():
    """Polled by the databases page while a task runs."""
    return jsonify({
        "task": db_tasks.snapshot() if db_tasks else None,
        "log": db_tasks.tail() if db_tasks else [],
        "run_busy": manager.busy() if manager else False,
    })


@app.route("/databases/log")
def databases_log():
    task = db_tasks.task() if db_tasks else None
    if not task or not os.path.exists(task.log_path):
        abort(404, "no database task has run yet")
    return send_file(task.log_path, mimetype="text/plain")


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
        description="Web interface for the Cancer Genomics Pipelines.",
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
    parser.add_argument("--data-dir", default=app.config["DATA_DIR"],
                        help="Where the installer put the references and "
                             "databases. Used for the database check, for "
                             "the upgrade actions, and as the source of the "
                             "form's suggested resource paths.")
    parser.add_argument("--no-db-check", dest="db_check", action="store_false",
                        help="Do not check for database updates at startup. "
                             "The front page then shows the last check's "
                             "result, however old it is.")
    parser.add_argument("--debug", action="store_true",
                        help="Flask debug mode. NEVER use this with real "
                             "patient data: the debugger allows arbitrary "
                             "code execution from the browser.")
    args = parser.parse_args()

    global manager
    data_dir = os.path.abspath(os.path.expanduser(args.data_dir))
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
    app.config["DATA_DIR"] = data_dir
    manager = JobManager(runs_dir)
    manager.load_existing()

    global db_tasks
    db_tasks = maintenance.Maintenance(
        log_dir=os.path.join(runs_dir, "_database_tasks"),
        script_dir=PIPELINE_DIR, data_dir=data_dir)

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
    print(f"[INFO] Data directory   : {data_dir}")
    print(f"[INFO] Open http://{args.host}:{args.port}")

    # Check the databases as the app comes up, in the background: the notice
    # on the front page is then about today rather than about whenever
    # somebody last remembered to run the checker. It only reads remote
    # metadata and writes one small status file, so it cannot disturb a run
    # -- and because it is a thread, a slow or unreachable source delays
    # nothing. In debug mode Flask's reloader runs main() twice; the child
    # process is the one that serves, so only it checks.
    serving = not args.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true"
    if args.db_check and serving:
        db_tasks.start_check(reason="startup")
        print("[INFO] Checking reference databases in the background; the "
              "front page shows the result when it lands.")
    elif not args.db_check:
        print("[INFO] Startup database check disabled (--no-db-check). The "
              "front page shows whatever the last check found.")

    app.run(host=args.host, port=args.port, debug=args.debug,
            threaded=True)


if __name__ == "__main__":
    main()
