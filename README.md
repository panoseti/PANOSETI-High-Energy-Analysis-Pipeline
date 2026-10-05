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

### Layout
* `heap/` - the analysis package. All pipeline logic lives here; notebooks and scripts only call into it.
* `configs/` - YAML configs. Paths in them are resolved relative to the project root unless absolute.
* `NextDay/` - the next-day quick-look pipeline.
* `analysis_tools/` - pff processing and analysis notebooks (see pff analysis) and `pcap_to_root_analysis.ipynb`.
* `example_notebooks/` - step-by-step examples of individual `heap` modules (pre-cleaning, pedestals, gain, image cleaning, coincidences, meridian flip correction).
* `test_data/` - the data the example notebooks and `pcap_to_root_analysis.ipynb` run on.
* `simulation_tools/` - CORSIKA simulation scripts and `panodisplay.C` (see simulation_tools below).
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
Configured by `configs/nextday.yaml` (data paths, module number -> telescope name, reference telescope and image cleaning thresholds).

### pff analysis
Four notebooks in `analysis_tools/` process and analyze one source over one or more nights of pff data. All four read the same config, a YAML like `configs/crab.yaml` (set by `CONFIG` at the top of each notebook), and share its `output_dir`:

1. `process_pff.ipynb` - image processing per telescope, per night: pedestals, pedvars and gain maps, cleaning and Hillas parameters (`heap.events.process_night()`, `heap.process_dataset`, `heap.image_cleaning`, `heap.parameterize`). Shows diagnostics of the inputs and outputs: each run's source/flip side/pointing, spike cuts, gain maps and which fallback made them, pedestals and pedvars over time, rates, and timing corrections and coincidences.
2. `inspect_data_products.ipynb` - what processing writes to `output_dir` and how to load it for further analysis (`heap.diagnostics.load_products()`, `heap.events.load_camera_frame()`, `heap.events.build_array_events()`).
3. `explore_events.ipynb` - Hillas parameter distributions, image centroid histograms and reconstructed events (direction histogram and event displays) for a selection of telescopes, nights and cuts, overridable in the notebook without editing the config.
4. `pff_analysis.ipynb` - the source analysis. Processes any nights not already in `output_dir`, then:
    * applies each camera's pointing correction, then corrects timestamps against the timing reference and builds array events from coincident images (`heap.coincidences`, `heap.events`)
    * applies image cuts, then reconstructs each event's arrival direction from the intersection of image axes (`heap.reconstruction`)
    * counts on/off regions (reflected or ring), computes Li & Ma significance and makes sky maps (`heap.significance`)

Processing writes, per night, telescope and source:
```
<output_dir>/<YYYYMMDD>/<telescope>/processing.json             settings used
<output_dir>/<YYYYMMDD>/<telescope>/<source>/<source>.npz        cleaned images + Hillas parameters, one row per frame
<output_dir>/<YYYYMMDD>/<telescope>/<source>/calibrations.npz    pedestals and pedvars per frame, gain map per night
```
Nights already in `output_dir` are not reprocessed, and processing them with different settings raises an error, so give configs with different telescopes or cleaning settings their own `output_dir`. Events, directions and significances are not written.

Each run's pointing, shared by every telescope, is its commanded wobble position: the source ±`source.wobble_offset` in Dec, or the source itself, whichever is closest to the hk mount positions (`heap.events.wobble_pointings()`). The hk mount positions aren't used as the pointing, since they drift away from it while the telescopes guide. Give runs commanded somewhere else, or with no hk mount data, a pointing in the config's `source.pointing_overrides`; runs with neither are left out (with a printout). When hk `target_name` is blank, the source is identified by matching the pointing against the catalog in `heap/sources.py`, so add new sources there. The timing reference (`events.reference`) is only used for timing: in each 120 s of each run, the first telescope in that list with data.

To analyze a new source, copy `configs/crab.yaml`, set `source.name` and `paths`, and point the notebooks' `CONFIG` at it.

## simulation_tools
### Dependencies
* [ROOT](http://root.cern.ch/). Verified for version 6.28/04
* [CORSIKA 7](https://www.iap.kit.edu/corsika/index.php). Verified for version 77410. Compiled with the following options enabled:
    * IACT
    * CHERENKOV
    * VOLUMEDET
    * SLANT
    * ATMEXT
    * QGSJET-II-04
    * UrQMD 1.3.1
* [This version](https://github.com/nkorzoun/corsikaIOreader) of corsikaIOreader

### Installation
* Install dependencies and be sure to compile corsikaIOreader with `make corsikaIOreader`

### Examples
Examples can be found on the [repo wiki](https://github.com/nkorzoun/PANOSETI-High-Energy-Analysis-Pipeline/wiki)
