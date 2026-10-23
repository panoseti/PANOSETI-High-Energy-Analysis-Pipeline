"""
CORSIKA -> sim_telarray simulation pipeline for PANOSETI, run as SLURM job arrays (or one run at
a time locally).

Configured by a YAML like configs/simulation.yaml: paths to CORSIKA, sim_telarray, the PANOSETI
sim_telarray config, the CORSIKA card template and the output/scratch directories, plus the
production's run numbers, CORSIKA keyword overrides, sim_telarray settings and SLURM resources.

Stages:
    full    CORSIKA for each run, then sim_telarray on its output (production.mode files), or
            CORSIKA piped straight into sim_telarray (production.mode pipe).
    simtel  sim_telarray only, on existing CORSIKA files (paths.corsika_input), e.g. to rerun a
            files-mode production with a changed sim_telarray config under a new tag.

On its way to sim_telarray, CORSIKA's output always goes through trace_photons.py, which traces
the photons through the lens (optical_model.py, settings in sim_telarray.optics) to the camera;
sim_telarray then skips its own optics (bypass_optics = 1 in its config, checked in its log) and
points at the zenith, so the real pointing is written next to the events.

prepare checks the config's paths and writes each run's CORSIKA card. submit (after prepare, for
the full stage) writes a SLURM array script for the stage, one array task per run or per CORSIKA
file, and submits it. Each submission's jobs read their own copy of the config, so editing the
config or submitting again (e.g. the simtel stage with a new tag) doesn't change queued jobs.

Output layout:
    <output_dir>/inputs/<name>_run<NNNNNN>.inp                 CORSIKA cards
    <output_dir>/corsika/<name>_run<NNNNNN>.corsika.zst        CORSIKA IACT files (files mode)
    <output_dir>/simtel/<tag>/<name>_run<NNNNNN>.simtel.zst    sim_telarray events
    <output_dir>/simtel/<tag>/<name>_run<NNNNNN>.hdata.zst     sim_telarray histograms
    <output_dir>/simtel/<tag>/<name>_run<NNNNNN>.pointing.json real pointing and optics settings
    <output_dir>/logs/<name>_run<NNNNNN>.corsika.log
    <output_dir>/logs/<tag>/<name>_run<NNNNNN>.simtel.log
    <output_dir>/logs/<tag>/<name>_run<NNNNNN>.trace.log       trace_photons.py
    <output_dir>/slurm/<job>.sbatch, <job>.yaml (config copy), <job>.files (simtel inputs)
    <output_dir>/slurm/logs/<job>_<SLURM job id>_<task>.out
with <job> = <stage>_<tag>_<YYYYmmdd-HHMMSS>.

Each run works in its own scratch directory (paths.work_dir) and its outputs are moved to
output_dir when it succeeds; a failed run keeps its scratch directory for debugging. Runs whose
outputs already exist are skipped unless --overwrite is given.

Usage:
    python simulation_tools/run_simulation.py prepare --config configs/simulation.yaml
    python simulation_tools/run_simulation.py submit --config configs/simulation.yaml
    python simulation_tools/run_simulation.py submit --config configs/simulation.yaml --stage simtel --dry-run
    python simulation_tools/run_simulation.py run --config configs/simulation.yaml --run 1
    python simulation_tools/run_simulation.py simtel --config configs/simulation.yaml --file some.corsika.zst
"""
import argparse
import glob
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import zlib
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

SECTIONS = {
    "paths": {"corsika_run_dir", "corsika_binary", "sim_telarray_dir", "sim_telarray_config",
              "corsika_template", "output_dir", "corsika_input", "work_dir", "processed_dir"},
    "production": {"name", "mode", "first_run", "n_runs", "seed"},
    "corsika": None, # keyword -> value
    "sim_telarray": {"telescope_theta", "telescope_phi", "wobble_offset", "wobble_direction", "tag", "defines",
                     "parameters", "optics"},
    "slurm": {"partition", "account", "time", "mem", "max_concurrent", "mail_user", "mail_type", "python",
              "extra_options", "setup"},
}
MODES = {"files", "pipe"}
STAGES = {"full", "simtel"}
# Keywords the pipeline sets itself, so they can't be overridden in the config's corsika section.
MANAGED_KEYWORDS = {"RUNNR", "SEED", "TELFIL", "DIRECT", "EXIT"}
MAX_RUNNR = 999999
CORSIKA_SUFFIX = ".corsika.zst"
PIPE_SCRIPT = "simtel_pipe.sh"
CORSIKA_DONE = "END OF RUN"
SIMTEL_DONE = "Sim_telarray finished"
# sim_telarray must take the traced photons' positions as camera positions (bypass_optics = 1 in
# its config); it says so in its log.
SIMTEL_BYPASS = "Bypassing optics"
TRACE_SCRIPT = Path(__file__).resolve().parent / "trace_photons.py"
TRACE_DONE = "trace_photons done"
OPTICS_KEYS = {"thick_lens", "focal_offset"}
# sim_telarray log lines repeated in the pipeline's output when a run finishes.
SIMTEL_SUMMARY = re.compile(r"^(Run\(s\) completed|Sim_telarray finished|\s*Tel\. triggered:)")
SHOWER_DONE = re.compile(r"END OF SHOWER NO\s+(\d+)")
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _check_keys(section, config, required):
    if not isinstance(config, dict):
        raise ValueError(f"config section {section!r}: expected a mapping, got {config!r}")
    missing = required - set(config)
    unknown = set(config) - required
    if missing or unknown:
        raise ValueError(f"config section {section!r}: missing {sorted(missing)}, unknown {sorted(unknown)}")


def _resolve(path):
    if path is None:
        return None
    path = Path(os.path.expandvars(str(path))).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_config(path):
    """
    Load and check a simulation config.

    Parameters:
        path: config YAML, see configs/simulation.yaml

    Returns:
        the config's sections as nested dicts, with paths.* as Path (relative paths resolved
        against the project root), except paths.corsika_input, a glob string (default
        <output_dir>/corsika/*.corsika.zst); paths.processed_dir and sim_telarray.optics None if null
        or missing (running sim_telarray needs optics); corsika keys upper-cased; and the sim_telarray defines/parameters and slurm
        extra_options/setup defaulting to empty
    """
    with open(path) as f:
        config = yaml.safe_load(f)
    if isinstance(config, dict) and isinstance(config.get("paths"), dict):
        config["paths"].setdefault("processed_dir", None)
    if isinstance(config, dict) and isinstance(config.get("sim_telarray"), dict):
        config["sim_telarray"].setdefault("optics", None)
    _check_keys("top level", config, set(SECTIONS))
    for section, required in SECTIONS.items():
        if required is None:
            config[section] = config[section] or {}
            if not isinstance(config[section], dict):
                raise ValueError(f"config section {section!r}: expected a mapping, got {config[section]!r}")
        else:
            _check_keys(section, config[section], required)

    paths = config["paths"]
    for key in paths:
        if key not in {"corsika_binary", "corsika_input"}:
            paths[key] = _resolve(paths[key])
    if paths["corsika_input"] is None:
        paths["corsika_input"] = str(paths["output_dir"] / "corsika" / f"*{CORSIKA_SUFFIX}")
    else:
        paths["corsika_input"] = str(_resolve(paths["corsika_input"]))

    production = config["production"]
    if production["mode"] not in MODES:
        raise ValueError(f"production.mode: {production['mode']!r} not one of {sorted(MODES)}")
    # The name goes into CORSIKA's TELFIL line, which CORSIKA cuts at blanks.
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", str(production["name"])):
        raise ValueError(f"production.name: {production['name']!r} may only contain letters, digits and _.+-")
    last_run = production["first_run"] + production["n_runs"] - 1
    if production["first_run"] < 1 or production["n_runs"] < 1 or last_run > MAX_RUNNR:
        raise ValueError(f"production: runs {production['first_run']}..{last_run} not within 1..{MAX_RUNNR}")

    config["corsika"] = {key.upper(): value for key, value in config["corsika"].items()}
    managed = MANAGED_KEYWORDS & set(config["corsika"])
    if managed:
        raise ValueError(f"corsika: {sorted(managed)} are set by the pipeline")

    sim_telarray = config["sim_telarray"]
    sim_telarray["defines"] = sim_telarray["defines"] or {}
    sim_telarray["parameters"] = sim_telarray["parameters"] or {}
    if sim_telarray["optics"] is not None:
        _check_keys("sim_telarray.optics", sim_telarray["optics"], OPTICS_KEYS)
    if "bypass_optics" in sim_telarray["parameters"]:
        raise ValueError("sim_telarray.parameters: bypass_optics belongs in the sim_telarray config (always 1, see trace_photons.py)")
    slurm = config["slurm"]
    slurm["extra_options"] = slurm["extra_options"] or []
    slurm["setup"] = slurm["setup"] or []
    return config


def check_paths(config):
    """Raise if any of the programs or files the config points to is missing."""
    paths = config["paths"]
    required = {
        "CORSIKA binary": paths["corsika_run_dir"] / paths["corsika_binary"],
        "sim_telarray binary": paths["sim_telarray_dir"] / "bin" / "sim_telarray",
        "sim_telarray config": paths["sim_telarray_config"],
        "CORSIKA template": paths["corsika_template"],
    }
    missing = [f"{what}: {path}" for what, path in required.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing " + "; ".join(missing))


def run_name(config, run):
    return f"{config['production']['name']}_run{run:06d}"


def _output_name(corsika_file):
    """<name>_run000001.corsika.zst -> <name>_run000001 (also DATxxx.telescope -> DATxxx)."""
    return re.sub(r"(\.(corsika|telescope|eventio|dat))?(\.(zst|gz|bz2|xz|lz4))?$", "", Path(corsika_file).name)


def _corsika_seeds(config, run, n):
    rng = np.random.default_rng([config["production"]["seed"], run, 0])
    return [int(seed) for seed in rng.integers(1, 900_000_000, n)]


def _simtel_seed(config, name):
    """Seeded by the output name, so a simtel-stage rerun of a run reuses its seed."""
    rng = np.random.default_rng([config["production"]["seed"], zlib.crc32(name.encode()), 1])
    return int(rng.integers(1, 2**31 - 1))


def _trace_seed(config, name):
    """trace_photons.py's seed (its wavelengths), like _simtel_seed()."""
    rng = np.random.default_rng([config["production"]["seed"], zlib.crc32(name.encode()), 2])
    return int(rng.integers(1, 2**31 - 1))


def _card_value(value):
    if isinstance(value, (list, tuple)):
        return " ".join(_card_value(v) for v in value)
    if isinstance(value, bool):
        return "T" if value else "F"
    return str(value)


def _override_lines(keyword, value):
    if keyword == "TELESCOPE":
        return [f"TELESCOPE {_card_value(row)}" for row in value]
    return [f"{keyword} {_card_value(value)}"]


def render_card(config, run):
    """
    The CORSIKA input card for one run: the template with RUNNR, SEED (derived from
    production.seed and the run number) and TELFIL set, and the config's corsika overrides
    applied (keywords not in the template are added before EXIT).
    """
    lines = [line.strip() for line in config["paths"]["corsika_template"].read_text().splitlines() if line.strip()]
    keywords = [line.split()[0].upper() for line in lines]
    for keyword in ("RUNNR", "SEED", "TELFIL", "EXIT"):
        if keyword not in keywords:
            raise ValueError(f"CORSIKA template {config['paths']['corsika_template']} has no {keyword} line")

    if config["production"]["mode"] == "pipe":
        telfil = f"|./{PIPE_SCRIPT}"
    else:
        telfil = f"./{run_name(config, run)}{CORSIKA_SUFFIX}"
    seeds = iter(_corsika_seeds(config, run, keywords.count("SEED")))
    overrides = config["corsika"]
    card = []
    replaced = set()
    for line, keyword in zip(lines, keywords):
        if keyword == "EXIT":
            for missing in overrides.keys() - set(keywords):
                card += _override_lines(missing, overrides[missing])
            card.append("EXIT")
            break
        if keyword == "RUNNR":
            card.append(f"RUNNR {run}")
        elif keyword == "SEED":
            card.append(f"SEED {next(seeds)} 0 0")
        elif keyword == "TELFIL":
            card.append(f"TELFIL {telfil}")
        elif keyword in overrides:
            # Repeated keywords (TELESCOPE) are replaced as a block at their first line.
            if keyword not in replaced:
                card += _override_lines(keyword, overrides[keyword])
                replaced.add(keyword)
        else:
            card.append(line)
    return "\n".join(card) + "\n"


def write_card(config, run):
    path = config["paths"]["output_dir"] / "inputs" / f"{run_name(config, run)}.inp"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_card(config, run))
    return path


def _simtel_value(value):
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def _corsika_header(config, corsika_file):
    """
    A CORSIKA file's number of telescopes and primary direction, from extract_corsika_tel.

    Returns:
        (n_telescopes, theta, phi): the direction is the middle of the THETAP/PHIP ranges (the
        VIEWCONE axis) in sim_telarray's convention, as in _card_direction()
    """
    command = [str(config["paths"]["sim_telarray_dir"] / "bin" / "extract_corsika_tel"),
               "--header-only", "--only-telescopes", "1-999", str(corsika_file)]
    header = subprocess.run(command, capture_output=True, text=True).stdout
    values = dict(re.findall(r"^(CORSIKA_TELESCOPES|CORSIKA_THETA|CORSIKA_PHI)=(\S+)$", header, re.MULTILINE))
    if len(values) < 3:
        raise RuntimeError(f"could not read the number of telescopes and primary direction from {corsika_file}")
    return int(values["CORSIKA_TELESCOPES"]), float(values["CORSIKA_THETA"]), float(values["CORSIKA_PHI"])


def _card_direction(card_text):
    """
    The primary direction of a CORSIKA card in sim_telarray's convention, computed like
    extract_corsika_tel: zenith = middle of THETAP, azimuth (from north, clockwise, where the
    primaries come from) = 180 + ARRANG - middle of PHIP (CORSIKA's direction of motion,
    counter-clockwise from magnetic north).
    """
    keywords = {line.split()[0].upper(): line.split()[1:] for line in card_text.splitlines() if line.strip()}
    theta = sum(float(v) for v in keywords.get("THETAP", ["0", "0"])[:2]) / 2
    phi_low, phi_high = (float(v) for v in keywords.get("PHIP", ["0", "0"])[:2])
    if phi_high < phi_low:
        phi_high += 360.
    phi = (180. + float(keywords.get("ARRANG", ["0"])[0]) - (phi_low + phi_high) / 2) % 360.
    return theta, phi


def _separation(theta1, phi1, theta2, phi2):
    """Angle [deg] between two directions (zenith angle, azimuth) [deg]."""
    theta1, phi1, theta2, phi2 = np.radians([theta1, phi1, theta2, phi2])
    cos_separation = np.cos(theta1) * np.cos(theta2) + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    return np.degrees(np.arccos(np.clip(cos_separation, -1., 1.)))


def _wobble(theta, phi, offset, direction):
    """
    The direction (zenith angle, azimuth) [deg] offset [deg] from (theta, phi) towards direction
    [deg] (from north, clockwise: 0 = N, 90 = E, 180 = S). North and east are the ground's,
    carried from the zenith to (theta, phi) by the rotation between them, so exact at the zenith;
    on the meridian, north and south are towards the zenith or horizon.
    """
    if offset == 0.:
        return theta, phi
    theta, phi, offset, direction = np.radians([theta, phi, offset, direction])
    source = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]) # x north, y east, z up
    axis = np.array([-np.sin(phi), np.cos(phi), 0.]) # zenith -> source is a rotation by theta about it
    def rotate(v):
        return v * np.cos(theta) + np.cross(axis, v) * np.sin(theta) + axis * (axis @ v) * (1. - np.cos(theta))
    towards = np.cos(direction) * rotate(np.array([1., 0., 0.])) + np.sin(direction) * rotate(np.array([0., 1., 0.]))
    wobbled = source * np.cos(offset) + towards * np.sin(offset)
    return (float(np.degrees(np.arccos(np.clip(wobbled[2], -1., 1.)))),
            float(np.degrees(np.arctan2(wobbled[1], wobbled[0]))) % 360.)


def _pointing(config, theta, phi):
    """
    Telescope pointing: sim_telarray.wobble_offset from the primaries' direction (theta, phi)
    towards wobble_direction (see _wobble()), with telescope_theta/phi overriding where not null.
    Warns if the pointing is not wobble_offset from the primaries' direction.
    """
    sim_telarray = config["sim_telarray"]
    wobble = _wobble(theta, phi, sim_telarray["wobble_offset"], sim_telarray["wobble_direction"])
    pointing = (wobble[0] if sim_telarray["telescope_theta"] is None else sim_telarray["telescope_theta"],
                wobble[1] if sim_telarray["telescope_phi"] is None else sim_telarray["telescope_phi"])
    offset = _separation(*pointing, theta, phi)
    if abs(offset - sim_telarray["wobble_offset"]) > 0.01:
        _say(f"WARNING: telescopes point at theta={pointing[0]:g}, phi={pointing[1]:g} deg, {offset:.2f} deg from the "
             f"CORSIKA primaries' direction theta={theta:g}, phi={phi:g} deg, not wobble_offset "
             f"{sim_telarray['wobble_offset']:g} deg")
    return pointing


def trace_command(config, input_file, name, pointing):
    """
    trace_photons.py command line: CORSIKA input (a file, or '-' for standard input) traced
    through the optics for telescopes pointing at pointing = (theta, phi) [deg] (see _pointing()),
    written to standard output for sim_telarray.
    """
    optics = config["sim_telarray"]["optics"]
    if optics is None:
        raise ValueError("config has no sim_telarray.optics")
    command = [sys.executable, str(TRACE_SCRIPT), str(input_file), "-", f"--theta={pointing[0]}", f"--phi={pointing[1]}",
               f"--focal-offset={optics['focal_offset']}", f"--seed={_trace_seed(config, name)}"]
    if optics["thick_lens"]:
        command.append("--thick-lens")
    return command


def simtel_command(config, name, output_file, hdata_file, n_telescopes):
    """
    sim_telarray command line for n_telescopes telescopes (sim_telarray otherwise reads only the
    first 4), reading trace_photons.py's output from standard input. No pointing: with
    bypass_optics sim_telarray points at the zenith whatever it's given.
    """
    paths = config["paths"]
    sim_telarray = config["sim_telarray"]
    command = [str(paths["sim_telarray_dir"] / "bin" / "sim_telarray"), f"-I{paths['sim_telarray_config'].parent}"]
    command += [f"-D{key}" if value is None else f"-D{key}={value}" for key, value in sim_telarray["defines"].items()]
    command += ["-c", str(paths["sim_telarray_config"])]
    parameters = {
        "maximum_telescopes": n_telescopes,
        "random_seed": _simtel_seed(config, name),
        **sim_telarray["parameters"],
        "output_file": output_file,
        "histogram_file": hdata_file,
        "show": "all",
    }
    for key, value in parameters.items():
        command += ["-C", f"{key}={_simtel_value(value)}"]
    command.append("-")
    return command


def _simtel_outputs(config, name):
    """
    Final (output_dir) paths of one sim_telarray run's events, histograms, pointing (the real one,
    which the events' header doesn't have, see write_pointing()) and logs (sim_telarray's and
    trace_photons.py's).
    """
    output_dir = config["paths"]["output_dir"]
    tag = config["sim_telarray"]["tag"]
    return {
        "events": output_dir / "simtel" / tag / f"{name}.simtel.zst",
        "histograms": output_dir / "simtel" / tag / f"{name}.hdata.zst",
        "pointing": output_dir / "simtel" / tag / f"{name}.pointing.json",
        "log": output_dir / "logs" / tag / f"{name}.simtel.log",
        "trace_log": output_dir / "logs" / tag / f"{name}.trace.log",
    }


def write_pointing(config, path, pointing):
    """
    The telescopes' pointing = (theta, phi) [deg] and the optics settings, as JSON: sim_telarray
    with bypass_optics records the zenith instead (read by heap.read_simtel.read_run_header()).
    """
    sim_telarray = config["sim_telarray"]
    path.write_text(json.dumps({
        "theta": pointing[0], "phi": pointing[1], "alt": 90. - pointing[0], "az": pointing[1],
        "wobble_offset": sim_telarray["wobble_offset"], "wobble_direction": sim_telarray["wobble_direction"],
        "optics": sim_telarray["optics"],
    }, indent=2) + "\n")


def _work_dir(config, name):
    base = config["paths"]["work_dir"]
    if base is None:
        base = Path(os.environ["TMPDIR"]) if os.environ.get("TMPDIR") else config["paths"]["output_dir"] / "work"
    # Retried, since another run finishing can remove an empty base (_remove_work) in between.
    for _ in range(10):
        base.mkdir(parents=True, exist_ok=True)
        try:
            return Path(tempfile.mkdtemp(prefix=f"{name}_", dir=base))
        except FileNotFoundError:
            pass
    raise RuntimeError(f"could not create a work directory in {base}")


def _remove_work(config, work):
    """
    Delete a finished run's scratch directory, then the work directory holding it if that's
    now empty (not when it's $TMPDIR, which belongs to the system or SLURM).
    """
    shutil.rmtree(work)
    if config["paths"]["work_dir"] is None and os.environ.get("TMPDIR"):
        return
    try:
        work.parent.rmdir()
    except OSError:
        pass # still holds other runs' directories, running or failed


def _move(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))


def _say(message):
    """Print with a timestamp, flushed so SLURM logs show progress as it happens."""
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {message}", flush=True)


def _duration(seconds):
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m"


def _run(command, log_path, cwd, stdin=None, verbose=False, on_line=None, close_stdin=False):
    """
    Run command in cwd, writing its output (stdout and stderr) to log_path, echoed to the
    terminal if verbose; on_line(line) is called for every output line. close_stdin: close our
    copy of stdin (another process's output pipe) once command has it, so that process sees a
    broken pipe if command exits early.

    Returns:
        exit code
    """
    if verbose:
        _say(f"$ {shlex.join(command)}")
    with open(log_path, "w") as log:
        process = subprocess.Popen(command, stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   cwd=cwd, text=True, errors="replace")
        if close_stdin:
            stdin.close()
        for line in process.stdout:
            log.write(line)
            if verbose:
                print(f"  | {line}", end="", flush=True)
            if on_line is not None:
                on_line(line)
        return process.wait()


def _run_traced(trace, command, trace_log, log_path, cwd, verbose=False):
    """
    trace | command in cwd: trace_photons.py (trace) feeding sim_telarray (command), with
    trace_photons.py's messages in trace_log and sim_telarray's output handled as by _run().

    Returns:
        trace_photons.py's exit code if not 0, else sim_telarray's
    """
    if verbose:
        _say(f"$ {shlex.join(trace)} 2> {trace_log.name} |")
    with open(trace_log, "w") as messages:
        tracer = subprocess.Popen(trace, stdout=subprocess.PIPE, stderr=messages, cwd=cwd)
        returncode = _run(command, log_path, cwd, stdin=tracer.stdout, verbose=verbose, close_stdin=True)
        return tracer.wait() or returncode


def _simtel_summary(log):
    for line in log.read_text(errors="replace").splitlines():
        if SIMTEL_SUMMARY.match(line):
            print(f"    {line.strip()}", flush=True)


def _finish_simtel(work, outputs, returncode=0):
    """
    Check a trace_photons.py | sim_telarray run in work finished, with sim_telarray taking the
    traced photons as camera positions (bypass_optics), then move its files to outputs.
    """
    log, trace_log = work / outputs["log"].name, work / outputs["trace_log"].name
    traced = trace_log.exists() and TRACE_DONE in trace_log.read_text(errors="replace")
    simtel_log = log.read_text(errors="replace") if log.exists() else ""
    if returncode != 0 or not traced or SIMTEL_DONE not in simtel_log or not (work / outputs["events"].name).exists():
        raise RuntimeError(f"trace_photons.py | sim_telarray failed (exit {returncode}), see {trace_log} and {log}; "
                           "work directory kept")
    if SIMTEL_BYPASS not in simtel_log:
        raise RuntimeError(f"sim_telarray traced the photons through its own optics again (no {SIMTEL_BYPASS!r} in {log}): "
                           "set bypass_optics = 1 in the sim_telarray config; work directory kept")
    _simtel_summary(log)
    print(f"    {trace_log.read_text(errors='replace').strip().splitlines()[-1]}", flush=True)
    for key in ("events", "histograms", "pointing", "log", "trace_log"):
        _move(work / outputs[key].name, outputs[key])


def run_simtel(config, corsika_file, overwrite=False, verbose=False):
    """
    Run sim_telarray on one CORSIKA file.

    Parameters:
        config: from load_config()
        corsika_file: CORSIKA IACT file; outputs are named after it with its extensions removed
        overwrite: rerun even if this tag's outputs already exist
        verbose: echo sim_telarray's command and output (always written to its log)

    Returns:
        {"events", "histograms", "log"}: output paths, see _simtel_outputs()
    """
    corsika_file = Path(corsika_file).resolve()
    name = _output_name(corsika_file)
    outputs = _simtel_outputs(config, name)
    if outputs["events"].exists() and not overwrite:
        _say(f"{outputs['events']} exists, skipping sim_telarray")
        return outputs
    n_telescopes, theta, phi = _corsika_header(config, corsika_file)
    pointing = _pointing(config, theta, phi)
    trace = trace_command(config, corsika_file, name, pointing)
    work = _work_dir(config, name)
    write_pointing(config, work / outputs["pointing"].name, pointing)
    command = simtel_command(config, name, work / outputs["events"].name, work / outputs["histograms"].name, n_telescopes)
    _say(f"trace_photons.py | sim_telarray: {corsika_file} ({n_telescopes} telescopes, pointing theta={pointing[0]:g}, "
         f"phi={pointing[1]:g} deg, tag {config['sim_telarray']['tag']})")
    _say(f"    work directory {work}")
    start = time.time()
    returncode = _run_traced(trace, command, work / outputs["trace_log"].name, work / outputs["log"].name, work, verbose)
    _say(f"trace_photons.py | sim_telarray finished in {_duration(time.time() - start)} (exit {returncode})")
    _finish_simtel(work, outputs, returncode)
    _remove_work(config, work)
    _say(f"sim_telarray outputs: {outputs['events']}, {outputs['histograms']}, {outputs['log']}")
    return outputs


def run_corsika(config, run, overwrite=False, verbose=False):
    """
    Run CORSIKA for one run, then sim_telarray: on its output file (files mode), or reading
    from CORSIKA through a pipe (pipe mode).

    Parameters:
        config: from load_config()
        run: CORSIKA run number
        overwrite: redo the run even if its outputs already exist
        verbose: echo CORSIKA's and sim_telarray's commands and output (always written to their
            logs; in pipe mode sim_telarray's output only goes to its log)
    """
    paths = config["paths"]
    name = run_name(config, run)
    pipe = config["production"]["mode"] == "pipe"
    corsika_file = paths["output_dir"] / "corsika" / f"{name}{CORSIKA_SUFFIX}"
    outputs = _simtel_outputs(config, name)
    if not overwrite and (outputs["events"].exists() if pipe else corsika_file.exists()):
        _say(f"{outputs['events'] if pipe else corsika_file} exists, skipping CORSIKA")
        if not pipe:
            run_simtel(config, corsika_file, overwrite, verbose)
        return

    card = write_card(config, run)
    work = _work_dir(config, name)
    # CORSIKA reads its data tables from the working directory. DAT* are earlier outputs.
    for item in paths["corsika_run_dir"].iterdir():
        if not item.name.startswith("DAT"):
            (work / item.name).symlink_to(item)
    if pipe:
        # CORSIKA starts this (TELFIL |./simtel_pipe.sh) and writes its output to its stdin.
        # TELFIL can't hold shell syntax, hence the script.
        n_telescopes = sum(line.split()[0].upper() == "TELESCOPE" for line in card.read_text().splitlines() if line.strip())
        pointing = _pointing(config, *_card_direction(card.read_text()))
        write_pointing(config, work / outputs["pointing"].name, pointing)
        trace = trace_command(config, "-", name, pointing)
        command = simtel_command(config, name, work / outputs["events"].name, work / outputs["histograms"].name, n_telescopes)
        script = work / PIPE_SCRIPT
        script.write_text("#!/bin/bash\nset -o pipefail\n"
                          f"{shlex.join(trace)} 2> {shlex.quote(str(work / outputs['trace_log'].name))} | "
                          f"{shlex.join(command)} > {shlex.quote(str(work / outputs['log'].name))} 2>&1\n")
        script.chmod(0o755)

    keywords = {line.split()[0].upper(): line.split()[1:] for line in card.read_text().splitlines() if line.strip()}
    n_showers = int(keywords["NSHOW"][0]) if "NSHOW" in keywords else None
    _say(f"run {run} ({name}), {config['production']['mode']} mode: CORSIKA {paths['corsika_binary']}, "
         f"{n_showers} showers, then sim_telarray (tag {config['sim_telarray']['tag']})")
    _say(f"    card {card}")
    _say(f"    work directory {work}")
    corsika_log = work / f"{name}.corsika.log"
    start = time.time()

    def progress(line):
        match = SHOWER_DONE.search(line)
        if match is None:
            return
        _say(f"CORSIKA shower {match.group(1)}/{n_showers or '?'} done ({_duration(time.time() - start)} elapsed)")

    with open(card) as stdin:
        returncode = _run([f"./{paths['corsika_binary']}"], corsika_log, work, stdin=stdin, verbose=verbose,
                          on_line=progress)
    _say(f"CORSIKA finished in {_duration(time.time() - start)} (exit {returncode})")
    if returncode != 0 or CORSIKA_DONE not in corsika_log.read_text(errors="replace"):
        raise RuntimeError(f"CORSIKA failed (exit {returncode}), see {corsika_log}; work directory kept")
    _move(corsika_log, paths["output_dir"] / "logs" / corsika_log.name)

    if pipe:
        _finish_simtel(work, outputs)
        _remove_work(config, work)
        _say(f"sim_telarray outputs: {outputs['events']}, {outputs['histograms']}, {outputs['log']}")
    else:
        _move(work / corsika_file.name, corsika_file)
        _remove_work(config, work)
        _say(f"CORSIKA output: {corsika_file}")
        run_simtel(config, corsika_file, overwrite, verbose)
    _say(f"run {run} done in {_duration(time.time() - start)}")


def prepare(config):
    """Check the config's paths and write every run's CORSIKA card to <output_dir>/inputs/."""
    check_paths(config)
    production = config["production"]
    for run in range(production["first_run"], production["first_run"] + production["n_runs"]):
        write_card(config, run)
    print(f"wrote {production['n_runs']} CORSIKA cards to {config['paths']['output_dir'] / 'inputs'}")


def write_sbatch(config, config_path, stage, verbose=False, overwrite=False):
    """
    Write a SLURM array script for one stage, one array task per run (full) or per CORSIKA file
    matching paths.corsika_input (simtel), as <output_dir>/slurm/<job>.sbatch with
    <job> = <stage>_<tag>_<YYYYmmdd-HHMMSS>. The jobs read a copy of the config, <job>.yaml, so
    later edits or submissions don't change them; the simtel stage's file list is frozen to
    <job>.files.

    Parameters:
        verbose: jobs echo CORSIKA's and sim_telarray's output to their SLURM logs
        overwrite: jobs redo outputs that already exist (each replaced when its run succeeds)

    Returns:
        (sbatch path, number of array tasks)
    """
    files = sorted(glob.glob(config["paths"]["corsika_input"])) if stage == "simtel" else None
    if files == []:
        raise FileNotFoundError(f"no CORSIKA files match {config['paths']['corsika_input']}")
    slurm = config["slurm"]
    slurm_dir = config["paths"]["output_dir"] / "slurm"
    (slurm_dir / "logs").mkdir(parents=True, exist_ok=True)
    job = f"{stage}_{config['sim_telarray']['tag']}_{datetime.now():%Y%m%d-%H%M%S}"
    # Relative paths resolve against the project root, not the config's location, so the copy
    # resolves the same as the original.
    config_copy = slurm_dir / f"{job}.yaml"
    shutil.copyfile(config_path, config_copy)
    python = slurm["python"] or sys.executable
    base = [python, str(Path(__file__).resolve())]
    flags = (" --verbose" if verbose else "") + (" --overwrite" if overwrite else "")
    if stage == "full":
        n_tasks = config["production"]["n_runs"]
        task = shlex.join(base + ["run", "--config", str(config_copy)]) + ' --index "$SLURM_ARRAY_TASK_ID"' + flags
    else:
        file_list = slurm_dir / f"{job}.files"
        file_list.write_text("\n".join(files) + "\n")
        n_tasks = len(files)
        task = (f'file="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {shlex.quote(str(file_list))})"\n'
                + shlex.join(base + ["simtel", "--config", str(config_copy)]) + ' --file "$file"' + flags)

    array = f"0-{n_tasks - 1}" + (f"%{slurm['max_concurrent']}" if slurm["max_concurrent"] else "")
    options = [
        f"--job-name={config['production']['name']}-{stage}",
        f"--array={array}",
        f"--output={slurm_dir / 'logs'}/{job}_%A_%a.out",
        "--ntasks=1",
        "--cpus-per-task=1",
        f"--time={slurm['time']}",
        f"--mem={slurm['mem']}",
    ]
    if slurm["partition"]:
        options.append(f"--partition={slurm['partition']}")
    if slurm["account"]:
        options.append(f"--account={slurm['account']}")
    if slurm["mail_user"]:
        options += [f"--mail-user={slurm['mail_user']}", f"--mail-type={slurm['mail_type']}"]
    options += slurm["extra_options"]
    lines = ["#!/bin/bash", *(f"#SBATCH {option}" for option in options), "", "set -e", *slurm["setup"], "", task]
    path = slurm_dir / f"{job}.sbatch"
    path.write_text("\n".join(lines) + "\n")
    return path, n_tasks


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    for command, help in [("prepare", "check paths and write the CORSIKA cards"),
                          ("submit", "prepare, write the SLURM array script and submit it"),
                          ("run", "full stage (CORSIKA + sim_telarray) for one run"),
                          ("simtel", "sim_telarray only, for one CORSIKA file")]:
        sub = commands.add_parser(command, help=help)
        sub.add_argument("--config", required=True, help="simulation config YAML, see configs/simulation.yaml")
        if command == "submit":
            sub.add_argument("--stage", choices=sorted(STAGES), default="full")
            sub.add_argument("--dry-run", action="store_true", help="write the sbatch script but don't submit it")
            sub.add_argument("--verbose", action="store_true", help="jobs echo CORSIKA's and sim_telarray's output to their SLURM logs")
            sub.add_argument("--overwrite", action="store_true", help="jobs redo outputs that already exist")
        if command in {"run", "simtel"}:
            which = sub.add_mutually_exclusive_group(required=True)
            if command == "run":
                which.add_argument("--run", type=int, help="CORSIKA run number")
                which.add_argument("--index", type=int, help="run first_run + index (SLURM array task id)")
            else:
                which.add_argument("--file", help="CORSIKA file")
                which.add_argument("--index", type=int, help="index into the sorted files matching paths.corsika_input")
            sub.add_argument("--overwrite", action="store_true", help="redo outputs that already exist")
            sub.add_argument("--verbose", action="store_true",
                             help="echo CORSIKA's and sim_telarray's commands and output (always written to their logs)")
    args = parser.parse_args()

    config = load_config(args.config)
    if args.command == "prepare":
        prepare(config)
    elif args.command == "submit":
        if args.stage == "full":
            prepare(config)
        else:
            check_paths(config)
        path, n_tasks = write_sbatch(config, args.config, args.stage, args.verbose, args.overwrite)
        print(f"wrote {path} ({n_tasks} array tasks)")
        if not args.dry_run:
            subprocess.run(["sbatch", str(path)], check=True)
    elif args.command == "run":
        check_paths(config)
        run = args.run if args.run is not None else config["production"]["first_run"] + args.index
        run_corsika(config, run, args.overwrite, args.verbose)
    else:
        check_paths(config)
        corsika_file = args.file
        if corsika_file is None:
            files = sorted(glob.glob(config["paths"]["corsika_input"]))
            if not 0 <= args.index < len(files):
                raise IndexError(f"--index {args.index}: {len(files)} files match {config['paths']['corsika_input']}")
            corsika_file = files[args.index]
        run_simtel(config, corsika_file, args.overwrite, args.verbose)


if __name__ == "__main__":
    main()
