# PANOSETI-High-Energy-Analysis-Pipeline
Tools to analyze images of air showers captured by [PANOSETI](https://panoseti.ucsd.edu/) telescopes.

## Data Pipeline (WIP)
### Dependencies
* [pypff](https://github.com/panoseti/pypff.git) 

### Installation
Clone this repo and install with:
```
pip install .
```

For development, you may want to install an editable version with:
```
pip install -e .
```
This will let you run jupyter notebooks after making changes to the package without needing to reinstall. You will still need to restart your kernel.

### This machine's paths
Where the data and simulations are differs from machine to machine, so it is set only in `configs/local.yaml`, which git ignores; the tracked configs and notebooks hold no machine paths. Copy the example and set your paths before running any notebook:
```
cp configs/local.example.yaml configs/local.yaml
```
* `raw_data_dir`: holds one `<YYYYMMDD>/` folder of runs per night.
* `processed_dir`: where processing writes. Each analysis config writes to `<processed_dir>/<its file name>/`, e.g. `<processed_dir>/crab/` for `configs/analysis/crab.yaml`.
* `analysis`: the analysis config the `analysis_tools/` notebooks and `simulation_tools/notebooks/data_mc.ipynb` load (`heap.config.local_analysis()`), e.g. `configs/analysis/crab.yaml`.
* `simulations`: the simulation configs the simulation notebooks load, by name (`heap.config.local_simulation()`), e.g. `gamma: configs/simulation/gamma_pointlike.yaml`: `gamma` for `process_sims.ipynb` and `explore_sims.ipynb`, `gamma` and `proton` for `gamma_hadron.ipynb` and `data_mc.ipynb`.

To run a notebook on another analysis config or production for a while, set its `CONFIG` (or `GAMMA_CONFIG`, `SIM_CONFIGS`, ...) to that file's path instead.

### Layout
* `heap/` - the analysis package. All pipeline logic lives here; notebooks and scripts only call into it.
* `configs/` - YAML configs. Paths in them are resolved relative to the project root unless absolute. At the top are this machine's `local.yaml` and the templates to copy (`*.example.*`); the folders hold the configs in use.
    * `analysis/` - analysis configs, e.g. `crab.yaml`.
    * `simulation/` - the simulation configs of each production, with their CORSIKA cards.
    * `nextday/` - the next-day pipeline's config.
    * `analysis/source_run_map.yaml` - each run's source and flip side, by run folder name, for runs the mount data can't identify (`heap.process_dataset.load_fallback_map()`). Kept by hand and shared by every analysis.
* `NextDay/` - the next-day quick-look pipeline.
* `analysis_tools/` - pff processing and analysis notebooks (see pff analysis) and `pcap_to_root_analysis.ipynb`.
* `example_notebooks/` - step-by-step examples of individual `heap` modules (pre-cleaning, pedestals, gain, image cleaning, coincidences, meridian flip correction).
* `test_data/` - the data the example notebooks and `pcap_to_root_analysis.ipynb` run on.
* `simulation_tools/` - the CORSIKA -> sim_telarray simulation pipeline (see simulation_tools below).
* `misc_tools/` - standalone notebooks (PSF, array layout, stellar spectra, transmission, ...).

Raw data is expected as one folder per night, holding one folder per run (acquisition session, e.g. either side of a meridian flip):
```
<raw_data_dir>/<YYYYMMDD>/[pff/]<run_id>.pffd/*.pff
```
Runs are read from `<YYYYMMDD>/pff/` if it exists (newer nights keep allsky, weather and obslogs alongside), otherwise directly from `<YYYYMMDD>/`.

### Next-day pipeline
Produces quick-look plots for one night: spike cuts, timing corrections, pairwise and triple coincidence rates, pedestals/pedvars/gains, Hillas parameter histograms and a sample of coincident events. Plots, a `manifest.json` and a static `index.html` gallery are written to `<output_dir>/<date>/`.
```
python NextDay/run_nextday_pipeline.py --date 20260116
python NextDay/run_nextday_pipeline.py --date 20260116 --config my_config.yaml
```
Configured by `configs/nextday/nextday.yaml` (data paths, module number -> telescope name, reference telescope and image cleaning thresholds).

### pff analysis
Four notebooks in `analysis_tools/` process and analyze one source over one or more nights of pff data. All four read the same analysis config, a YAML like `configs/analysis/crab.yaml` (`configs/local.yaml`'s `analysis`, or `CONFIG` at the top of each notebook), and share its output folder, `output_dir` below (`<processed_dir>/<config name>/`):

1. `process_pff.ipynb` - image processing per telescope, per night: pedestals, pedvars and gain maps, cleaning and Hillas parameters (`heap.events.process_night()`, then `heap.events.update_gain_maps()` for nights where a closer night's gain map has become available, `heap.process_dataset`, `heap.image_cleaning`, `heap.parameterize`). Shows diagnostics of the inputs and outputs: each run's source/flip side/pointing, spike cuts, gain maps and which fallback made them, pedestals and pedvars over time, rates, and timing corrections and coincidences.
2. `inspect_data_products.ipynb` - what processing writes to `output_dir` and how to load it for further analysis (`heap.diagnostics.load_products()`, `heap.events.load_camera_frame()`, `heap.events.build_array_events()`).
3. `explore_events.ipynb` - Hillas parameter distributions, image centroid histograms and reconstructed events (direction histogram and event displays) for a selection of telescopes, nights and cuts, overridable in the notebook without editing the config.
4. `pff_analysis.ipynb` - the source analysis. Processes any nights not already in `output_dir`, then:
    * applies each camera's pointing correction, then corrects timestamps against the timing reference and builds array events from coincident images (`heap.coincidences`, `heap.events`)
    * applies image cuts, then reconstructs each event's arrival direction from the intersection of image axes (`heap.reconstruction`), and its shower core from the same axes through the telescopes' positions (the config's `telescopes.<module>.position`), with the camera's axes on the ground from the config's `site` (`heap.events.sky_frame()`)
    * counts on/off regions (reflected or ring), computes Li & Ma significance and makes sky maps (`heap.significance`)

Processing writes, per night, telescope and source:
```
<output_dir>/<YYYYMMDD>/<telescope>/processing.json             settings used
<output_dir>/<YYYYMMDD>/<telescope>/<source>/<source>.npz        cleaned images + Hillas parameters, one row per frame
<output_dir>/<YYYYMMDD>/<telescope>/<source>/calibrations.npz    pedestals and pedvars per frame, gain map per night
```
Nights already in `output_dir` are not reprocessed, and processing them with different settings raises an error, so to process with other telescopes or cleaning settings, copy the config under a new name: its output goes to a new folder. Events, directions and significances are not written.

Each run's pointing, shared by every telescope, is the mount positions at its start, read from each telescope's mount log (`<night>/<telescope>/mount/`, else hk.pff), since every run starts right after Ekos aligns the telescopes; later positions drift as they follow the guide corrections (`heap.events.run_pointings()`). Give runs that didn't start with an alignment (e.g. a DAQ restart), or have no mount data, a pointing in the config's `source.pointing_overrides`; runs with neither are left out (with a printout). When the mount has no usable target name, a run's source comes from `configs/analysis/source_run_map.yaml` if listed there (add runs by folder name), else by matching the pointing against the catalog in `heap/sources.py`, so add new sources there. The timing reference (`events.reference`) is only used for timing: in each 120 s of each run, the first telescope in that list with data.

To analyze a new source, copy `configs/analysis/crab.yaml`, set `source.name`, and point `configs/local.yaml`'s `analysis` at it.

## simulation_tools
CORSIKA -> sim_telarray simulations of the array with the [PANOSETI sim_telarray config](https://github.com/orelgueta/PANOSETI_sim_telarray_config), run as SLURM job arrays. Configured by `configs/simulation.example.yaml` (tool and data paths, runs, CORSIKA overrides, sim_telarray settings, SLURM resources):
```
python simulation_tools/run_simulation.py submit --config my_simulation.yaml
```
See [simulation_tools/README.md](simulation_tools/README.md) for requirements, the config, stages and outputs.
