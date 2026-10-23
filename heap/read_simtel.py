"""read_simtel

Functions for reading sim_telarray (.simtel.zst) files into numpy arrays, in the same layout as the
pff and root data: one (N, 32, 32) image per telescope, indexed [row, col], with columns along the
camera's x and rows along its y (sim_telarray's pixel id is row * 32 + col). Images are ADC with
sim_telarray's artificial offset (fadc_pedestal, there so negative noise isn't clipped) removed, and
NSB and noise left in, like pff data. sim_telarray's pedestal frames (pedestal_events) are kept separate
from the shower images, for the pedestals and pedvars (see heap.process_simulation).

Telescopes are sim_telarray's telescope ids (1 ... n, the order of the CORSIKA TELESCOPE lines).
Events are matched across telescopes by event_id, not by timestamps: a simulated event's times are
when sim_telarray processed it.
"""

import numpy as np
from eventio import SimTelFile

SKY_PEDESTAL = 1 # calibration_type of sim_telarray's pedestal_events (lid open, with NSB)

def _as_list(simtel_files):
    """One path, or a list of paths -> list of paths"""
    if isinstance(simtel_files, (str, bytes)) or not hasattr(simtel_files, "__iter__"):
        return [simtel_files]
    return list(simtel_files)

def _tel_image(pixels):
    """(1024,) pixel values in sim_telarray's pixel order -> (32, 32) image indexed [row, col]"""
    return np.asarray(pixels).reshape(32, 32)

def read_events(simtel_files, true_images=False):
    """
    Read every telescope's images and the simulated showers from one or more .simtel(.zst) files

    Parameters:
        simtel_files: path to a .simtel(.zst) file, or a list of them (e.g. a production's runs)
        true_images: also return each image's true number of photoelectrons per pixel (before noise and electronics)

    Returns:
        telescopes: {telescope_id: dict} with
            - images: the telescope's triggered shower images, ADC minus fadc_pedestal, numpy array with shape (N, 32, 32)
            - run, event_id: each image's run and event, the keys into showers
            - true_pe: (if true_images) Cherenkov photoelectrons per pixel of each image, shape (N, 32, 32),
              NaN for images whose telescope got fewer than sim_telarray's store_photoelectrons
            - pedestal_frames: the telescope's sky pedestal events (sim_telarray's pedestal_events: NSB
              and noise only, no shower), ADC minus fadc_pedestal, shape (M, 32, 32), for the pedestals and pedvars
        showers: dict of numpy arrays, one entry per shower event (with at least one telescope triggered):
            - run, event_id
            - energy: primary energy [TeV]
            - alt, az: primary direction [deg], az from north, clockwise
            - core_x, core_y: core position on the ground [m], x north and y west (CORSIKA's frame)
            - h_first_int: height of the first interaction [m]
            - xmax: depth of the shower maximum [g/cm^2]
            - n_triggered: number of telescopes with an image
    """
    telescopes = {}
    showers = {key: [] for key in ["run", "event_id", "energy", "alt", "az", "core_x", "core_y", "h_first_int", "xmax", "n_triggered"]}

    for simtel_file in _as_list(simtel_files):
        with SimTelFile(str(simtel_file)) as f:
            run = int(f.header["run"])
            tel_ids = f.header["tel_id"] # photoelectrons are keyed by the telescope's index in this list
            for event in f:
                pedestal = event["type"] == "calibration"
                if pedestal and event["calibration_type"] != SKY_PEDESTAL:
                    continue # dark, single-p.e. and flat-field events
                if not pedestal:
                    shower = event["mc_shower"]
                    mc_event = event["mc_event"]
                    showers["run"].append(run)
                    showers["event_id"].append(event["event_id"])
                    showers["energy"].append(shower["energy"])
                    showers["alt"].append(np.degrees(shower["altitude"]))
                    showers["az"].append(np.degrees(shower["azimuth"]))
                    showers["core_x"].append(mc_event["xcore"])
                    showers["core_y"].append(mc_event["ycore"])
                    showers["h_first_int"].append(shower["h_first_int"])
                    showers["xmax"].append(shower["xmax"])
                    showers["n_triggered"].append(len(event["telescope_events"]))

                if true_images and not pedestal:
                    true_pe = {int(tel_ids[index]): pe["photoelectrons"] for index, pe in event["photoelectrons"].items()}

                for tel_id, tel_event in event["telescope_events"].items():
                    if "adc_sums" not in tel_event:
                        raise ValueError(f"{simtel_file}: event {event['event_id']} telescope {tel_id} has no ADC sums")
                    tel = telescopes.setdefault(tel_id, {"images": [], "run": [], "event_id": []} | ({"true_pe": []} if true_images else {}) | {"pedestal_frames": []})
                    offset = event["camera_monitorings"][tel_id]["pedestal"][0] # fadc_pedestal * fadc_sum_bins, per pixel
                    image = _tel_image(tel_event["adc_sums"][0] - offset).astype(np.float32)
                    if pedestal:
                        tel["pedestal_frames"].append(image)
                        continue
                    tel["images"].append(image)
                    tel["run"].append(run)
                    tel["event_id"].append(event["event_id"])
                    if true_images: # NaN if sim_telarray didn't store them (fewer than store_photoelectrons)
                        tel["true_pe"].append(_tel_image(true_pe[tel_id]).astype(np.float32) if tel_id in true_pe else np.full((32, 32), np.nan, dtype=np.float32))

    telescopes = {
        tel_id: {key: np.array(values) if len(values) else np.empty((0, 32, 32) if key in ("images", "true_pe", "pedestal_frames") else 0)
                 for key, values in tel.items()}
        for tel_id, tel in sorted(telescopes.items())
    }
    showers = {key: np.array(values) for key, values in showers.items()}
    return telescopes, showers

def read_tel_data(simtel_file, telescope_id):
    """
    Read the images of a single telescope from a .simtel(.zst) file and return a numpy array

    Parameters:
        simtel_file: .simtel(.zst) file written by sim_telarray
        telescope_id: sim_telarray telescope id

    Returns:
        data: ADC minus fadc_pedestal, numpy array with shape (N, 32, 32)
    """

    telescopes, _ = read_events(simtel_file)
    if telescope_id not in telescopes:
        return np.empty((0, 32, 32), dtype=np.float32)
    return telescopes[telescope_id]["images"]

def read_run_header(simtel_file):
    """
    Read the run's pointing, telescope positions and simulated showers from a .simtel(.zst) file

    Parameters:
        simtel_file: .simtel(.zst) file written by sim_telarray

    Returns dict with:
        - run: run number
        - alt, az: array pointing [deg], az from north, clockwise
        - tel_id: telescope ids, shape (n_tel,)
        - tel_pos: telescope positions [m], shape (n_tel, 3), x north, y west, z up (CORSIKA's frame)
        - n_showers, n_use: CORSIKA showers in the run and times each is reused, so n_showers * n_use events were simulated
        - energy_range: [TeV], spectral_index
        - core_range: [m], cores are simulated within core_range[1] of the array origin
        - viewcone: [deg], primary directions are simulated between viewcone[0] and viewcone[1] off the cone's axis
        - obsheight: observation level [m]
    """

    with SimTelFile(str(simtel_file)) as f:
        header = f.header
        mc_header = f.mc_run_headers[0]

    return {
        "run": int(header["run"]),
        "alt": float(np.degrees(header["direction"][1])),
        "az": float(np.degrees(header["direction"][0])),
        "tel_id": np.array(header["tel_id"], dtype=int),
        "tel_pos": np.array(header["tel_pos"], dtype=float),
        "n_showers": int(mc_header["n_showers"]),
        "n_use": int(mc_header["n_use"]),
        "energy_range": np.array(mc_header["E_range"], dtype=float),
        "spectral_index": float(mc_header["spectral_index"]),
        "core_range": np.array(mc_header["core_range"], dtype=float),
        "viewcone": np.array(mc_header["viewcone"], dtype=float),
        "obsheight": float(mc_header["obsheight"]),
    }
