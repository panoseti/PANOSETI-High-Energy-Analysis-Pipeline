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
* `analysis_tools/` - notebooks, including `pff_analysis.ipynb` for source analyses.

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

### Source analysis
`analysis_tools/pff_analysis.ipynb` analyzes one source over one or more nights, configured by a YAML like `configs/crab.yaml`:
1. Clean images and compute Hillas parameters per telescope, per night (`heap.process_dataset`, `heap.image_cleaning`, `heap.parameterize`). Nights already in `output_dir` are not reprocessed.
2. Correct timestamps against the reference telescope and build array events from coincident images (`heap.coincidences`, `heap.events`).
3. Apply image cuts, then reconstruct each event's arrival direction from the intersection of image axes (`heap.reconstruction`).
4. Count on/off regions (reflected or ring), compute Li & Ma significance and make sky maps (`heap.significance`).

Each run's pointing comes from the reference telescope's hk mount table. When hk `target_name` is blank, the source is identified by matching the pointing against the catalog in `heap/sources.py`, so add new sources there. Runs with missing or wrong mount data can be given a pointing in the config's `source.pointing_overrides`.

To analyze a new source, copy `configs/crab.yaml`, set `source.name` and `paths`, and point the notebook at it. Give configs with different telescopes or cleaning settings their own `output_dir`, since processed nights are reused.

## simulation-tools
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
