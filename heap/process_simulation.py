"""process_simulation

The pipeline of heap.process_dataset.process_dataset() for sim_telarray files: the same image processing
(process_image()) and the same outputs, so simulated images are analysed like data. The steps that only
apply to real data are left out: there are no packet loss, rate spikes, meridian flips or mount logs in a
simulation, the pedestals and pedvars are calculated from sim_telarray's pedestal frames (pedestal_events)
instead of the shower frames, and the gain is flat (1).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from heap.make_pedestals import calculate_pedestal_and_pedvar
from heap.process_dataset import process_image
from heap.read_simtel import read_events

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

    Parameters:
        simtel_files: .simtel(.zst) files of the production, simulated with sim_telarray's pedestal_events
            (at least MIN_PEDESTAL_FRAMES pedestal frames per telescope in total)
        out_dir: directory to write to (created if missing)
        name: name of the production, used for the file names
        nsig, fit_gaussian: passed through to make_pedestals.calculate_pedestal_and_pedvar()
        x, y, image_threshold, border_threshold, keep_brightest_island: passed through to process_image()

    Returns:
        {telescope_id: (cleaned_images, params_df, raw_images)}, one entry per telescope with images.
        raw_images is the pre-cleaning image for each row of cleaned_images/params_df.
        showers: the MC truth, see read_simtel.read_events()
    """

    out_dir = Path(out_dir)
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

    return results, showers
