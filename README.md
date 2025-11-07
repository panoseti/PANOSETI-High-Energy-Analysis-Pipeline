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
