"""process_simulation

The pipeline of heap.process_dataset.process_dataset() for sim_telarray files: the same image processing
(process_image()) and the same outputs, so simulated images are analysed like data. The steps that only
apply to real data are left out: there are no packet loss, rate spikes, meridian flips or mount logs in a
simulation, the pedestals and pedvars are calculated from sim_telarray's pedestal frames (pedestal_events)
instead of the shower frames, and the gain is flat (1).

Also loads the processed outputs back into array events (load_events()) and plots them against the MC truth.
Camera coordinates are sim_telarray's, fixed to the alt-az mount: +x towards higher altitude, +y towards
lower azimuth (measured from the images' displacement towards the shower core), see camera_position().
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from heap import events as ev
from heap.make_pedestals import calculate_pedestal_and_pedvar
from heap.parameterize import CAMERA_CMAP, CAMERA_EXTENT, CAMERA_HALF_WIDTH, TAB10_COLORS, draw_params
from heap.process_dataset import process_image
from heap.read_simtel import read_events, read_run_header

MIN_PEDESTAL_FRAMES = 100 # as build_calibrations() needs per chunk of data

def process_simulation(
        simtel_files,
        out_dir,
        name: str,
        nsig: float = 5.0,
        fit_gaussian: bool = True,
        x: float = None,
        y: float = None,
        image_threshold: float = 4.0,
        border_threshold: float = 2.0,
        keep_brightest_island: bool = False,
        reprocess: bool = False,
):
    """
    Run the pipeline for every telescope of a simulation production.

    For each telescope, calculates one pedestal and pedvar from its pedestal frames of all simtel_files
    (make_pedestals.calculate_pedestal_and_pedvar()), with a flat gain, then runs process_image() over
    every shower image. Writes <out_dir>/T<telescope id>/<name>/<name>.npz per telescope, holding
    cleaned_images plus one array per params_df column, as process_dataset() does, next to that
    telescope's calibrations.npz (as build_calibrations() writes it), and <out_dir>/<name>_showers.npz with
    the MC truth (see read_simtel.read_events()). params_df has Run and EventID, the keys into the showers,
    in place of process_dataset()'s Timestamp.

    The settings (nsig, fit_gaussian, x, y and the cleaning settings) are recorded in
    <out_dir>/<name>_processing.json, as heap.events.process_night() does for the data: if out_dir
    already has this production, it is reused when it was processed with the same settings, and made
    with other (or unrecorded) settings it raises, unless reprocess.

    Parameters:
        simtel_files: .simtel(.zst) files of the production, simulated with sim_telarray's pedestal_events
            (at least MIN_PEDESTAL_FRAMES pedestal frames per telescope in total)
        out_dir: directory to write to (created if missing)
        name: name of the production, used for the file names
        nsig, fit_gaussian: passed through to make_pedestals.calculate_pedestal_and_pedvar()
        x, y, image_threshold, border_threshold, keep_brightest_island: passed through to process_image()
        reprocess: rerun even if out_dir already has this production

    Returns:
        {telescope_id: (cleaned_images, params_df, raw_images)}, one entry per telescope with images.
        raw_images is the pre-cleaning image for each row of cleaned_images/params_df.
        showers: the MC truth, see read_simtel.read_events()
        Both None if the existing output was reused.
    """

    out_dir = Path(out_dir)
    settings = {
        "nsig": nsig, "fit_gaussian": fit_gaussian, "x": x, "y": y,
        "image_threshold": image_threshold, "border_threshold": border_threshold,
        "keep_brightest_island": keep_brightest_island,
    }
    settings_path = out_dir / f"{name}_processing.json"
    if not reprocess and (out_dir / f"{name}_showers.npz").exists():
        previous = json.loads(settings_path.read_text()) if settings_path.exists() else None
        if previous == settings:
            print(f"{out_dir} already has {name}, processed with the same settings: reusing it")
            return None, None
        raise ValueError(
            f"{out_dir} was processed with {previous or 'unrecorded settings'}, not {settings}; "
            "reprocess, or use another out_dir"
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    telescopes, showers = read_events(simtel_files)
    np.savez(out_dir / f"{name}_showers.npz", **showers)

    results = {}

    for tel_id, tel in telescopes.items():
        if len(tel["images"]) == 0:
            print(f"  Skipping telescope {tel_id}: no images")
            continue
        pedestal_frames = tel["pedestal_frames"].reshape(-1, 1024)
        if len(pedestal_frames) < MIN_PEDESTAL_FRAMES:
            raise ValueError(f"telescope {tel_id}: {len(pedestal_frames)} pedestal frames, need at least {MIN_PEDESTAL_FRAMES} "
                             "(set sim_telarray's pedestal_events, in the simulation config's sim_telarray parameters)")

        pedestal, pedvar = calculate_pedestal_and_pedvar(pedestal_frames, nsig=nsig, fit_gaussian=fit_gaussian)
        gain = np.ones(1024, dtype=np.float32)
        print(f"  Telescope {tel_id}: pedestal and pedvar from {len(pedestal_frames)} pedestal frames, flat gain, {len(tel['images'])} images")

        raw_images = tel["images"].reshape(-1, 1024)
        tel_dir = out_dir / f"T{tel_id}" / name
        tel_dir.mkdir(parents=True, exist_ok=True)
        np.savez(
            tel_dir / "calibrations.npz",
            pedestals=np.broadcast_to(pedestal.astype(np.float32), raw_images.shape),
            pedestal_variances=np.broadcast_to(pedvar.astype(np.float32), raw_images.shape),
            gain=gain.reshape(32, 32),
            gain_source="flat",
            gain_caption="Flat gain (1.0) - simulation",
        )

        cleaned_images = []
        all_params = []

        for i, (event, run, event_id) in enumerate(zip(raw_images, tel["run"], tel["event_id"])):
            cleaned, params = process_image(
                event, pedestal, pedvar, gain,
                x=x, y=y,
                image_threshold=image_threshold,
                border_threshold=border_threshold,
                keep_brightest_island=keep_brightest_island,
            )

            cleaned_images.append(cleaned)

            params['Event'] = i
            params['Run'] = run
            params['EventID'] = event_id
            all_params.append(params)

        cleaned_images = np.stack(cleaned_images, axis=0).astype(np.float32)
        params_df = pd.DataFrame(all_params)

        np.savez(
            tel_dir / f"{name}.npz",
            cleaned_images=cleaned_images,
            **{col: params_df[col].values for col in params_df.columns},
        )

        results[tel_id] = (cleaned_images, params_df, raw_images)

    settings_path.write_text(json.dumps(settings, indent=2))
    return results, showers

def run_table(simtel_files):
    """
    read_simtel.read_run_header() of every file, one row per run (indexed by run), with the file's
    path and size (MB)
    """
    if len(simtel_files) == 0:
        raise ValueError("no .simtel(.zst) files (check the config's output_dir, production.name and sim_telarray.tag)")
    rows = [read_run_header(f) | {"File": Path(f), "MB": Path(f).stat().st_size / 1e6} for f in simtel_files]
    return pd.DataFrame(rows).set_index("run")

def product_table(out_dir, name: str, telescopes: list):
    """
    What process_simulation() wrote per telescope: Images, Cleaned (the fraction with any pixel left
    after cleaning), the camera's median pedestal and pedvar, and Settings (<name>_processing.json, the
    same for every telescope)

    Parameters:
        out_dir, name: as passed to process_simulation()
        telescopes: "T<telescope id>"
    """
    settings_path = Path(out_dir) / f"{name}_processing.json"
    settings = json.loads(settings_path.read_text()) if settings_path.exists() else None
    rows = []
    for tel in telescopes:
        tel_dir = Path(out_dir) / tel / name
        if not (tel_dir / f"{name}.npz").exists():
            rows.append({"Telescope": tel, "Images": 0, "Settings": settings})
            continue
        n_pix = np.load(tel_dir / f"{name}.npz")["N_pix"]
        calib = np.load(tel_dir / "calibrations.npz")
        rows.append({
            "Telescope": tel,
            "Images": len(n_pix),
            "Cleaned": np.mean(n_pix > 0),
            "Pedestal": np.median(calib["pedestals"][0]),
            "Pedvar": np.median(calib["pedestal_variances"][0]),
            "Settings": settings,
        })
    return pd.DataFrame(rows).set_index("Telescope")

def plot_calibrations(out_dir, name: str, telescopes: list):
    """
    Each telescope's pedestal and pedvar maps (one per production, from its pedestal frames), each
    row on one color scale

    Parameters:
        out_dir, name: as passed to process_simulation()
        telescopes: "T<telescope id>", one column each

    Returns:
        the Figure
    """
    telescopes = [tel for tel in telescopes if (Path(out_dir) / tel / name / "calibrations.npz").exists()]
    calibs = {tel: np.load(Path(out_dir) / tel / name / "calibrations.npz") for tel in telescopes}
    fig, axs = plt.subplots(2, len(telescopes), figsize=(3.6*len(telescopes), 6.5), squeeze=False, layout="constrained")
    for row, (key, label) in enumerate([("pedestals", "Pedestal (ADC)"), ("pedestal_variances", "Pedvar (ADC)")]):
        maps = {tel: calib[key][0].reshape(32, 32) for tel, calib in calibs.items()}
        vmin, vmax = np.percentile(np.stack(list(maps.values())), [1, 99])
        for ax, (tel, values) in zip(axs[row], maps.items()):
            im = ax.imshow(values, origin="lower", extent=CAMERA_EXTENT, cmap=CAMERA_CMAP, vmin=vmin, vmax=vmax)
            ax.set(title=f"{tel} {label.split()[0].lower()}", xlabel="X (degrees)", ylabel="Y (degrees)")
        fig.colorbar(im, ax=axs[row], label=label, shrink=0.8)
    return fig

def plot_image_grid(out_dir, name: str, telescope: str, indices, n_cols: int = 4):
    """
    heap.diagnostics.plot_image_grid() for one telescope of a simulation: cleaned images with their
    Hillas ellipse, titled with their run and event id

    Parameters:
        out_dir, name: as passed to process_simulation()
        telescope: "T<telescope id>"
        indices: rows of the telescope's cleaned_images to show

    Returns:
        the Figure
    """
    npz = np.load(Path(out_dir) / telescope / name / f"{name}.npz")
    params = pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"})
    images = npz["cleaned_images"]
    n_rows = int(np.ceil(len(indices) / n_cols))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(3.5*n_cols, 3.2*n_rows), squeeze=False)
    for ax, idx in zip(axs.flat, indices):
        row = params.iloc[idx]
        im = ax.imshow(images[idx], origin="lower", extent=CAMERA_EXTENT, cmap=CAMERA_CMAP, vmin=0)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.sca(ax)
        draw_params(fig, row, lw=1.5)
        ax.set_title(f"{idx}: N_pix={int(row.N_pix)}, size={row['size']:.0f}\nrun {row.Run:.0f}, event id {row.EventID:.0f}", fontsize=8)
    for ax in axs.flat[len(indices):]:
        ax.axis("off")
    fig.suptitle(telescope)
    fig.tight_layout()
    return fig

def camera_position(alt, az, pointing_alt, pointing_az):
    """
    (alt, az) [deg] -> sim_telarray camera coordinates [deg], gnomonic about the pointing: +x towards
    higher altitude, +y towards lower azimuth
    """
    alt, az, pointing_alt, pointing_az = (np.radians(np.asarray(a, dtype=float)) for a in (alt, az, pointing_alt, pointing_az))
    daz = az - pointing_az
    cos_c = np.sin(alt) * np.sin(pointing_alt) + np.cos(alt) * np.cos(pointing_alt) * np.cos(daz)
    x = (np.sin(alt) * np.cos(pointing_alt) - np.cos(alt) * np.sin(pointing_alt) * np.cos(daz)) / cos_c
    y = -np.cos(alt) * np.sin(daz) / cos_c
    return np.degrees(x), np.degrees(y)

def impact_distance(showers, pos):
    """
    Distance [m] of a telescope from each shower's axis

    Parameters:
        showers: rows with alt, az [deg] and core_x, core_y [m], see read_simtel.read_events()
        pos: the telescope's position [m], x north, y west, z up (CORSIKA's frame)
    """
    alt, az = np.radians(showers.alt.to_numpy()), np.radians(showers.az.to_numpy())
    axis = np.stack([np.cos(alt) * np.cos(az), -np.cos(alt) * np.sin(az), np.sin(alt)], axis=-1) # towards the primary
    r = np.asarray(pos)[None, :] - np.stack([showers.core_x, showers.core_y, np.zeros(len(showers))], axis=-1)
    return np.linalg.norm(r - (r * axis).sum(axis=1)[:, None] * axis, axis=1)

def load_events(out_dir, name: str, headers, telescopes: list = None):
    """
    Read process_simulation()'s outputs back as array events, in the layout of
    heap.events.build_array_events(), for heap.events.apply_cuts() and heap.reconstruction.

    Parameters:
        out_dir, name: as passed to process_simulation()
        headers: DataFrame of read_simtel.read_run_header() per run, indexed by run
        telescopes: "T<telescope id>" to load (default = every telescope in headers)

    Returns:
        images: DataFrame with one row per image: Telescope, ImageIndex (row in the telescope's
            cleaned_images), Event (unique per Run and EventID), Run, EventID, Date (= name, which
            reconstruct_directions() carries), the Hillas parameters, and impact, the telescope's
            distance [m] from the shower axis (impact_distance())
        showers: DataFrame of the MC truth (read_simtel.read_events()) with Event, Run, EventID and
            true_x, true_y, the primary's direction in camera coordinates (camera_position())
    """
    out_dir = Path(out_dir)
    header = headers.iloc[0]
    tel_pos = {f"T{tel_id}": pos for tel_id, pos in zip(header.tel_id, header.tel_pos)}

    showers = pd.DataFrame(dict(np.load(out_dir / f"{name}_showers.npz"))).rename(columns={"run": "Run", "event_id": "EventID"})
    pointing = headers.loc[showers.Run, ["alt", "az"]].to_numpy()
    showers["true_x"], showers["true_y"] = camera_position(showers.alt, showers.az, pointing[:, 0], pointing[:, 1])
    showers.insert(0, "Event", np.arange(len(showers)))

    frames = []
    for tel in telescopes or tel_pos:
        path = out_dir / tel / name / f"{name}.npz"
        if not path.exists():
            print(f"{tel}: no images")
            continue
        npz = np.load(path)
        p = pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"})
        p = p.rename(columns={"Event": "ImageIndex"}).merge(showers[["Run", "EventID", "Event"]], on=["Run", "EventID"])
        p.insert(0, "Telescope", tel)
        p["impact"] = impact_distance(showers.set_index("Event").loc[p.Event], tel_pos[tel])
        frames.append(p)
    images = pd.concat(frames, ignore_index=True)
    images["Date"] = name
    return images, showers

def simulated_counts(headers, energy_bins):
    """
    Events simulated per energy bin over every run: n_showers * n_use per run, on the power law of
    energy_range and spectral_index (see read_simtel.read_run_header())
    """
    header = headers.iloc[0]
    e_min, e_max = header.energy_range
    index = header.spectral_index # dN/dE ~ E^index
    clipped = np.clip(energy_bins, e_min, e_max) # none simulated outside energy_range
    frac = (clipped[1:]**(index + 1) - clipped[:-1]**(index + 1)) / (e_max**(index + 1) - e_min**(index + 1))
    return (headers.n_showers * headers.n_use).sum() * frac

def plot_energies(headers, samples: dict, energy_bins):
    """
    Energy histograms of the simulated events and of each sample, and each sample's effective area
    (pi * core_range[1]^2 times the fraction of simulated events in it)

    Parameters:
        headers: DataFrame of read_simtel.read_run_header() per run
        samples: {label: showers}, e.g. {"triggered": ..., "after cuts": ...}, rows of load_events()'s showers
        energy_bins: bin edges [TeV]

    Returns:
        the Figure
    """
    simulated = simulated_counts(headers, energy_bins)
    area_simulated = np.pi * headers.iloc[0].core_range[1]**2 # m^2
    centers = np.sqrt(energy_bins[1:] * energy_bins[:-1])

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    axs[0].stairs(simulated, energy_bins, color="grey", label=f"simulated N={simulated.sum():.0f}")
    for i, (label, sample) in enumerate(samples.items()):
        counts, _ = np.histogram(sample.energy, energy_bins)
        axs[0].stairs(counts, energy_bins, color=f"C{i}", label=f"{label} N={len(sample)}")
        with np.errstate(divide="ignore", invalid="ignore"):
            area = area_simulated * counts / simulated
            err = area_simulated * np.sqrt(counts) / simulated
        axs[1].errorbar(centers, area, yerr=err, fmt="o", color=f"C{i}", label=label)
    axs[0].set(xscale="log", yscale="log", xlabel="energy (TeV)", ylabel="events", title="Energy")
    axs[1].set(xscale="log", yscale="log", xlabel="energy (TeV)", ylabel="effective area (m$^2$)", title="Effective area")
    for ax in axs:
        ax.legend()
    return fig

def plot_cores(headers, samples: dict, colors: dict = None):
    """
    Core positions on the ground, east right and north up, one panel per sample, with the telescopes
    and the circle cores are simulated in

    Parameters:
        headers: DataFrame of read_simtel.read_run_header() per run
        samples: {label: showers}, see plot_energies()
        colors: optional {"T<telescope id>": color}

    Returns:
        the Figure
    """
    header = headers.iloc[0]
    colors = colors or {}
    fig, axs = plt.subplots(1, len(samples), figsize=(6.5*len(samples), 6), sharex=True, sharey=True, squeeze=False)
    for ax, (label, sample) in zip(axs[0], samples.items()):
        ax.scatter(-sample.core_y, sample.core_x, s=2, alpha=0.4, color="grey") # CORSIKA's x north, y west
        for tel_id, (x, y, _) in zip(header.tel_id, header.tel_pos):
            ax.scatter(-y, x, marker="^", s=80, color=colors.get(f"T{tel_id}", "k"), edgecolor="k", zorder=3)
            ax.annotate(f"T{tel_id}", (-y, x), xytext=(5, 5), textcoords="offset points")
        ax.add_patch(plt.Circle((0, 0), header.core_range[1], fill=False, ls="--", color="k"))
        ax.set(aspect="equal", xlabel="east (m)", ylabel="north (m)", title=f"Cores, {label} (N={len(sample)})")
    return fig

def plot_multiplicity(samples: dict, n_telescopes: int):
    """
    Images per event of each sample

    Parameters:
        samples: {label: images}, rows of load_events()'s images
        n_telescopes: number of telescopes, the largest multiplicity

    Returns:
        the Figure
    """
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for i, (label, sample) in enumerate(samples.items()):
        ax.hist(sample.groupby("Event").size(), bins=np.arange(0.5, n_telescopes + 1.5), histtype="step", color=f"C{i}", label=label)
    ax.set(xlabel="images per event", ylabel="events", title="Multiplicity")
    ax.legend()
    return fig

def plot_image_params(images, showers, title="", colors: dict = None):
    """
    Each image's size against the true energy, and its length/width against its distance (from the
    camera center, see process_simulation())

    Parameters:
        images, showers: from load_events() (images e.g. after heap.events.apply_cuts())
        colors: optional {"T<telescope id>": color}

    Returns:
        the Figure
    """
    colors = colors or {}
    images = images.merge(showers[["Event", "energy"]], on="Event")
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for name, tel in images.groupby("Telescope"):
        axs[0].scatter(tel.energy, tel["size"], s=4, alpha=0.5, color=colors.get(name), label=name)
        axs[1].scatter(tel.distance, tel.length / tel.width, s=4, alpha=0.5, color=colors.get(name), label=name)
    axs[0].set(xscale="log", yscale="log", xlabel="energy (TeV)", ylabel="size (ADC)", title="Size vs energy")
    axs[1].set(xlabel="distance (deg)", ylabel="length / width", title="Length/width vs distance")
    axs[1].legend(markerscale=4)
    fig.suptitle(title)
    return fig

def true_alpha(images, showers):
    """
    Each image's alpha [deg] to its shower's true direction (true_x, true_y, see load_events()): the
    angle between the image's major axis and the line from its centroid to the true direction. The
    same as the Hillas alpha without a wobble, which is from the camera center.
    """
    true = showers.set_index("Event").loc[images.Event, ["true_x", "true_y"]].to_numpy()
    dx, dy = true[:, 0] - images.x_c.to_numpy(), true[:, 1] - images.y_c.to_numpy()
    phi = np.radians(images.phi.to_numpy())
    cos_alpha = np.abs(np.cos(phi) * dx + np.sin(phi) * dy) / np.hypot(dx, dy)
    return np.degrees(np.arccos(np.clip(cos_alpha, 0, 1)))

def plot_alpha(samples: dict, showers, colors: dict = None, bins: int = 18):
    """
    Alpha to the true direction (true_alpha()), one panel per sample, one line per telescope and all
    telescopes pooled filled in grey

    Parameters:
        samples: {label: images}, rows of load_events()'s images, e.g. {"before cuts": ..., "after cuts": ...}
        showers: load_events()'s showers
        colors: optional {"T<telescope id>": color}
        bins: bins over 0-90 deg

    Returns:
        the Figure
    """
    colors = colors or {}
    fig, axs = plt.subplots(1, len(samples), figsize=(6.5*len(samples), 4.5), sharey=True, squeeze=False)
    for ax, (label, images) in zip(axs[0], samples.items()):
        alpha = true_alpha(images, showers)
        ax.hist(alpha, bins=bins, range=(0, 90), density=True, histtype="stepfilled", color="black", alpha=0.4,
                label=f"All telescopes N={len(images)}")
        for i, (name, tel) in enumerate(images.groupby("Telescope")):
            ax.hist(alpha[(images.Telescope == name).to_numpy()], bins=bins, range=(0, 90), density=True, histtype="step",
                    color=colors.get(name, TAB10_COLORS[i % len(TAB10_COLORS)]), label=f"{name} N={len(tel)}")
        ax.set(xlabel=r"$\alpha$ to the true direction (deg)", ylabel="normalized counts", title=f"Alpha, {label}")
    axs[0, -1].legend(fontsize=8)
    return fig

def plot_directions(directions, energy_bins, radius: float = 3, min_events: int = 10):
    """
    Reconstructed directions in the camera, their theta^2 to the true direction, and the 68%
    containment radius per energy bin

    Parameters:
        directions: heap.reconstruction.reconstruct_directions() merged with load_events()'s showers
            (energy, true_x, true_y), with theta2
        energy_bins: bin edges [TeV]
        radius: deg from the camera center for the direction histogram, and sqrt of the theta^2 range
        min_events: energy bins with fewer reconstructed events are left out

    Returns:
        the Figure
    """
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    edges = np.linspace(-radius, radius, 61)
    _, _, _, im = axs[0].hist2d(directions.Xoffset, directions.Yoffset, bins=[edges, edges], cmap=CAMERA_CMAP)
    fig.colorbar(im, ax=axs[0], label="events")
    axs[0].scatter(directions.true_x.mean(), directions.true_y.mean(), marker="*", facecolor="none", edgecolor="w", s=200, label="true (mean)")
    axs[0].set(aspect="equal", xlabel="X (degrees)", ylabel="Y (degrees)", title="Reconstructed directions")
    axs[0].legend()

    axs[1].hist(directions.theta2, bins=60, range=(0, radius**2), histtype="step")
    axs[1].axvline(directions.theta2.quantile(0.68), color="k", ls="--", label="68%")
    axs[1].set(xlabel=r"$\theta^2$ (deg$^2$)", ylabel="events", title=r"$\theta^2$ to the true direction")
    axs[1].legend()

    centers = np.sqrt(energy_bins[1:] * energy_bins[:-1])
    by_energy = directions.groupby(pd.cut(directions.energy, energy_bins), observed=False)
    r68 = np.sqrt(by_energy.theta2.quantile(0.68)).to_numpy()
    n = by_energy.size().to_numpy()
    axs[2].scatter(centers[n >= min_events], r68[n >= min_events])
    axs[2].set(xscale="log", xlabel="energy (TeV)", ylabel="68% containment (deg)",
               title=f"Angular resolution (bins with {min_events}+ events)")
    return fig

def plot_event(event, direction, shower, out_dir, name: str, colors: dict = None):
    """
    heap.events.plot_event() for a simulated event: each cleaned image, then every image's Hillas
    ellipse on one camera plane, all with the reconstructed and true directions.

    Parameters:
        event: one event's images, see heap.events.apply_cuts()
        direction: that event's row from heap.reconstruction.reconstruct_directions()
        shower: that event's row of load_events()'s showers
        out_dir, name: as passed to process_simulation(), for the cleaned images
        colors: optional {"T<telescope id>": color} for the combined Hillas panel

    Returns:
        the Figure
    """
    half = CAMERA_HALF_WIDTH
    colors = colors or {}
    source_xy = (shower.true_x, shower.true_y)

    fig, axs = plt.subplots(1, len(event) + 1, figsize=(4*(len(event) + 1), 4), squeeze=False)
    for ax, (_, tel) in zip(axs[0], event.iterrows()):
        img = ev._cleaned_images(Path(out_dir) / tel.Telescope / name / f"{name}.npz")[tel.ImageIndex]
        im = ax.imshow(img, origin="lower", extent=[-half, half, -half, half], cmap=CAMERA_CMAP, vmin=0)
        fig.colorbar(im, ax=ax, label="ADU", fraction=0.046, pad=0.04)
        ev.draw_hillas(ax, tel, "white")
        ev.draw_positions(ax, direction, "white", source_xy, "true")
        ax.set(xlim=[-half, half], ylim=[-half, half], xlabel="X (degrees)", ylabel="Y (degrees)")
        ax.set_title(f"{tel.Telescope}, impact {tel.impact:.0f} m, size {tel['size']:.0f}", fontsize=9)

    ax = axs[0, -1]
    for i, (_, tel) in enumerate(event.iterrows()):
        ev.draw_hillas(ax, tel, colors.get(tel.Telescope, TAB10_COLORS[i % len(TAB10_COLORS)]), fill_alpha=0.3)
    ev.draw_positions(ax, direction, "k", source_xy, "true")
    ax.set(xlim=[-half, half], ylim=[-half, half], aspect="equal", xlabel="X (degrees)", ylabel="Y (degrees)")
    ax.set_title("Hillas parameterization", fontsize=9)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.suptitle(f"Event {direction.Event} (run {shower.Run:.0f}, event id {shower.EventID:.0f}), {shower.energy:.1f} TeV", fontsize=10)
    fig.tight_layout()
    return fig
