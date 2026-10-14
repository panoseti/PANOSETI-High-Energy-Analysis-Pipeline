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

prepare checks the config's paths and writes each run's CORSIKA card. submit (after prepare, for
the full stage) writes a SLURM array script for the stage, one array task per run or per CORSIKA
file, and submits it. Each submission's jobs read their own copy of the config, so editing the
config or submitting again (e.g. the simtel stage with a new tag) doesn't change queued jobs.

Output layout:
    <output_dir>/inputs/<name>_run<NNNNNN>.inp                 CORSIKA cards
    <output_dir>/corsika/<name>_run<NNNNNN>.corsika.zst        CORSIKA IACT files (files mode)
    <output_dir>/simtel/<tag>/<name>_run<NNNNNN>.simtel.zst    sim_telarray events
    <output_dir>/simtel/<tag>/<name>_run<NNNNNN>.hdata.zst     sim_telarray histograms
    <output_dir>/logs/<name>_run<NNNNNN>.corsika.log
    <output_dir>/logs/<tag>/<name>_run<NNNNNN>.simtel.log
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
              "corsika_template", "output_dir", "corsika_input", "work_dir"},
    "production": {"name", "mode", "first_run", "n_runs", "seed"},
    "corsika": None, # keyword -> value
    "sim_telarray": {"telescope_theta", "telescope_phi", "tag", "defines", "parameters"},
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
        <output_dir>/corsika/*.corsika.zst); corsika keys upper-cased; and the sim_telarray
        defines/parameters and slurm extra_options/setup defaulting to empty
    """
    with open(path) as f:
        config = yaml.safe_load(f)
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


def _n_telescopes(config, corsika_file):
    """Number of telescopes in a CORSIKA file's header."""
    command = [str(config["paths"]["sim_telarray_dir"] / "bin" / "extract_corsika_tel"),
               "--header-only", "--only-telescopes", "1-999", str(corsika_file)]
    header = subprocess.run(command, capture_output=True, text=True).stdout
    match = re.search(r"^CORSIKA_TELESCOPES=(\d+)$", header, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"could not read the number of telescopes from {corsika_file}")
    return int(match.group(1))


def simtel_command(config, input_file, name, output_file, hdata_file, n_telescopes):
    """
    sim_telarray command line for one CORSIKA input (a file, or '-' for standard input) with
    n_telescopes telescopes (sim_telarray otherwise reads only the first 4).
    """
    paths = config["paths"]
    sim_telarray = config["sim_telarray"]
    command = [str(paths["sim_telarray_dir"] / "bin" / "sim_telarray"), f"-I{paths['sim_telarray_config'].parent}"]
    command += [f"-D{key}" if value is None else f"-D{key}={value}" for key, value in sim_telarray["defines"].items()]
    command += ["-c", str(paths["sim_telarray_config"])]
    parameters = {
        "maximum_telescopes": n_telescopes,
        "telescope_theta": sim_telarray["telescope_theta"],
        "telescope_phi": sim_telarray["telescope_phi"],
        "random_seed": _simtel_seed(config, name),
        **sim_telarray["parameters"],
        "output_file": output_file,
        "histogram_file": hdata_file,
        "show": "all",
    }
    for key, value in parameters.items():
        command += ["-C", f"{key}={_simtel_value(value)}"]
    command.append(str(input_file))
    return command


def _simtel_outputs(config, name):
    """Final (output_dir) paths of one sim_telarray run's events, histograms and log."""
    output_dir = config["paths"]["output_dir"]
    tag = config["sim_telarray"]["tag"]
    return {
        "events": output_dir / "simtel" / tag / f"{name}.simtel.zst",
        "histograms": output_dir / "simtel" / tag / f"{name}.hdata.zst",
        "log": output_dir / "logs" / tag / f"{name}.simtel.log",
    }


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


def _run(command, log_path, cwd, stdin=None, verbose=False, on_line=None):
    """
    Run command in cwd, writing its output (stdout and stderr) to log_path, echoed to the
    terminal if verbose; on_line(line) is called for every output line.

    Returns:
        exit code
    """
    if verbose:
        _say(f"$ {shlex.join(command)}")
    with open(log_path, "w") as log:
        process = subprocess.Popen(command, stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   cwd=cwd, text=True, errors="replace")
        for line in process.stdout:
            log.write(line)
            if verbose:
                print(f"  | {line}", end="", flush=True)
            if on_line is not None:
                on_line(line)
        return process.wait()


def _simtel_summary(log):
    for line in log.read_text(errors="replace").splitlines():
        if SIMTEL_SUMMARY.match(line):
            print(f"    {line.strip()}", flush=True)


def _finish_simtel(work, outputs, returncode=0):
    """Check a sim_telarray run in work finished, then move its files to outputs."""
    log = work / outputs["log"].name
    finished = log.exists() and SIMTEL_DONE in log.read_text(errors="replace")
    if returncode != 0 or not finished or not (work / outputs["events"].name).exists():
        raise RuntimeError(f"sim_telarray failed (exit {returncode}), see {log}; work directory kept")
    _simtel_summary(log)
    for key in ("events", "histograms", "log"):
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
    work = _work_dir(config, name)
    n_telescopes = _n_telescopes(config, corsika_file)
    command = simtel_command(config, corsika_file, name, work / outputs["events"].name, work / outputs["histograms"].name,
                             n_telescopes)
    _say(f"sim_telarray: {corsika_file} ({n_telescopes} telescopes, tag {config['sim_telarray']['tag']})")
    _say(f"    work directory {work}")
    start = time.time()
    returncode = _run(command, work / outputs["log"].name, work, verbose=verbose)
    _say(f"sim_telarray finished in {_duration(time.time() - start)} (exit {returncode})")
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
        command = simtel_command(config, "-", name, work / outputs["events"].name, work / outputs["histograms"].name,
                                 n_telescopes)
        script = work / PIPE_SCRIPT
        script.write_text(f"#!/bin/bash\nexec {shlex.join(command)} > {shlex.quote(str(work / outputs['log'].name))} 2>&1\n")
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


def write_sbatch(config, config_path, stage, verbose=False):
    """
    Write a SLURM array script for one stage, one array task per run (full) or per CORSIKA file
    matching paths.corsika_input (simtel), as <output_dir>/slurm/<job>.sbatch with
    <job> = <stage>_<tag>_<YYYYmmdd-HHMMSS>. The jobs read a copy of the config, <job>.yaml, so
    later edits or submissions don't change them; the simtel stage's file list is frozen to
    <job>.files.

    Parameters:
        verbose: jobs echo CORSIKA's and sim_telarray's output to their SLURM logs

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
    flags = " --verbose" if verbose else ""
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
        path, n_tasks = write_sbatch(config, args.config, args.stage, args.verbose)
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
