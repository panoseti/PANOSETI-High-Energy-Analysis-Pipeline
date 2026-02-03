# PANOSETI Docker Container

This directory contains a complete Docker environment for the PANOSETI High Energy Analysis Pipeline. The container includes all necessary tools to run CORSIKA air shower simulations and analyze them using ROOT.

## Overview

The Docker container provides:
- **ROOT v6.28.10**: Data analysis framework
- **CORSIKA 7.7550**: Air shower simulation with QGSJET-II-04
- **corsikaIOreader Tool**: Convert CORSIKA output to ROOT format
- **Example simulation**: Pre-configured CORSIKA example ready to visualize

## Quick Start

### Building the Container

**Build with CORSIKA credentials (required):**

CORSIKA requires registration at https://www.iap.kit.edu/corsika/. Once you have credentials:

```bash
cd /path/to/PANOSETI-High-Energy-Analysis-Pipeline/container
docker build \
  --build-arg CORSIKA_USER=your_username \
  --build-arg CORSIKA_PASSWORD=your_password \
  -t panoseti:latest \
  .
```

**Using the Makefile:**
```bash
# Set credentials as environment variables
export CORSIKA_USER=your_username
export CORSIKA_PASSWORD=your_password

# Build
make build

# Or inline:
make CORSIKA_USER=user CORSIKA_PASSWORD=pass build
```

**Security Note:** CORSIKA credentials are NOT stored in this repository. They are only passed as build-time arguments and not embedded in any files or the final image layers.

**Alternative: Pre-download CORSIKA manually:**

If you want to avoid passing credentials to Docker:

1. Download CORSIKA manually:
   ```bash
   cd container
   wget --user=YOUR_USER --password=YOUR_PASS \
     https://web.iap.kit.edu/corsika/download/corsika-v7750/corsika-77550.tar.gz
   ```

2. Build the container (will use the existing tarball):
   ```bash
   docker build -t panoseti:latest .
   ```

### Testing the Container

After building, run the test suite to verify all components work:

```bash
make test
```

This checks:
- ROOT installation
- CORSIKA executable
- corsikaIOreader
- Example files
- Runs a quick CORSIKA simulation

### Quick Start: Run Example Simulation

The Makefile provides convenient one-command targets for running simulations and visualizations.

**Run example simulation and view visualization:** (requires X11 forwarding)
```bash
xhost +local:docker
make example
```
*Note: The default example uses ATMOD 6 (US Standard Atmosphere). For ATMOSPHERE 61 T (external profile), additional setup is needed.*

**Run example simulation only (batch mode):**
```bash
make example-sim
```

**Visualize pre-built example.root (requires X11):**
```bash
xhost +local:docker
make example-vis
```
*This visualizes the pre-computed 10-event gamma shower simulation.*

**Interactive shell:**

```bash
docker run -it --rm panoseti:latest
```

**With volume mounts for data persistence:**

```bash
docker run -it --rm \
  -v $(pwd)/runs:/opt/panoseti/data \
  -v $(pwd)/results:/opt/panoseti/results \
  panoseti:latest
```

## Using the Example Simulation

The container includes a pre-built example simulation demonstrating CORSIKA air showers.

### Method 1: View with ROOT (Interactive)

```bash
docker run -it --rm \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  panoseti:latest

# Inside container:
cd /opt/panoseti/simulation-tools
root -l run_example.C
```

This displays event 0 from the 10-event simulation.

### Method 2: View Event Information

```bash
docker run --rm panoseti:latest cat /opt/panoseti/example/example.inp
```

### Method 3: Run ROOT commands interactively

```bash
docker run -it --rm panoseti:latest

# Inside container:
cd /opt/panoseti/simulation-tools
root -l
> gSystem->Load("libGenVector")
> .L panodisplay.C+
> readFile("example/example.root")
> panodisplay(0)  # Display events 0-9
> .q
```

## Running New CORSIKA Simulations

### CORSIKA Atmosphere Models

The container CORSIKA is built with ATMEXT (external atmosphere) enabled. Use **internal** atmosphere models (ATMOD 1-6) for simple simulations:

| ATMOD | Description | Example |
|-------|-------------|---------|
| 1 | CUS (Summer) | `ATMOD 1` |
| 2 | CUS (Winter) | `ATMOD 2` |
| 3 | MSIS (Summer) | `ATMOD 3` |
| 4 | MSIS (Winter) | `ATMOD 4` |
| 5 | US Standard (Summer) | `ATMOD 5` |
| 6 | **US Standard (Annual) - Recommended** | `ATMOD 6` |

**Note:** The `ATMOSPHERE 61 T` line in the original example references an external atmosphere profile. This requires additional data files beyond the standard CORSIKA installation. For most use cases, use `ATMOD 6` instead of `ATMOSPHERE 61 T`.

### Step 1: Create a CORSIKA Input File

```bash
# Start container with volume mount
docker run -it --rm \
  -v $(pwd)/my_simulation:/opt/panoseti/data \
  panoseti:latest

# Inside container, create input file
cat > /opt/panoseti/data/my_simulation.inp << 'EOF'
RUNNR 1
EVTNR 1
NSHOW 10
PRMPAR 1
ERANGE 3E4 1E5
ESLOPE 0
THETAP 0. 30.
PHIP 0. 360.
SEED 805 0 0
SEED 3234 0 0
SEED 305 0 0
SEED 2220 0 0
ATMOD 1
MAGNET 25.2 40.88
ARRANG 12.77
ELMFLG F T
RADNKG 200.E2
FIXCHI 0.
HADFLG 0 0 0 0 0 2
QGSJET T 0
QGSSIG T
HILOW 100.
ECUTS 0.30 0.05 0.02 0.02
MUADDI F
MUMULT T
LONGI T 20. F F
MAXPRT 50
PAROUT F F
ECTMAP 1.E6
DEBUG F 6 F 1000000
DIRECT ./
USER user
HOST host
ATMOSPHERE 61 T
TELFIL ./example.telescope
OBSLEV 1239.E2
CSCAT 1 100.E2 100.E2
CERFIL 0
CERSIZ 5.
CWAVLG 200. 700.
TELESCOPE 53.59E2 73.52E2 1E2 0.25E2
TELESCOPE 53.59E2 -80.48E2 1E2 0.25E2
TELESCOPE -107.18E2 6.95E2 1E2 0.25E2
EXIT
EOF
```

### Step 2: Run CORSIKA

```bash
cd /opt/panoseti/corsika-77550/run
./corsika77550Linux_QGSII_urqmd < /opt/panoseti/data/my_simulation.inp > /opt/panoseti/data/output.log
```

Output files:
- `DAT*`: Particle detector data
- `CER*`: Cherenkov photon data
- `TEC*`: Telescope configuration
- `TLU*`: Longitudinal profile

### Step 3: Convert CORSIKA Output to ROOT

```bash
corsikaIOreader -cors DAT000001 -histo /opt/panoseti/data/output.root -abs CORSIKA
```

### Step 4: Visualize Simulation

```bash
cd /opt/panoseti/simulation-tools
root -l
> gSystem->Load("libGenVector")
> .L panodisplay.C+
> readFile("/opt/panoseti/data/output.root")
> panodisplay(0)
```

## CORSIKA Input File Parameters

Important parameters in `.inp` files:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `RUNNR` | Run number | `1` |
| `NSHOW` | Number of showers to simulate | `10` |
| `PRMPAR` | Primary particle type (1=gamma, 14=proton) | `1` |
| `ERANGE` | Energy range in GeV | `3E4 1E5` (30-100 TeV) |
| `THETAP` | Zenith angle range in degrees | `0. 30.` |
| `PHIP` | Azimuth angle range in degrees | `0. 360.` |
| `TELESCOPE` | Telescope position (x, y, z, radius) | `53.59E2 73.52E2 1E2 0.25E2` |
| `OBSLEV` | Observation altitude in cm | `1239.E2` (1239 m) |
| `CWAVLG` | Cherenkov wavelength range in nm | `200. 700.` |

## Container Structure

```
/opt/panoseti/
├── corsika-77550/             # CORSIKA source and executable
│   └── run/
│       └── corsika77550Linux_QGSII_urqmd
├── corsikaIOreader/           # Source code for reader tool
├── simulation-tools/          # Analysis and visualization scripts
│   ├── panodisplay.C         # Visualization macro
│   ├── run_example.C         # Example runner
│   └── example/              # Pre-built example simulation
│       ├── example.root      # Simulation data (10 events)
│       ├── example.inp       # CORSIKA input configuration
│       └── example.telescope # Telescope configuration
├── data/                      # Volume mount for CORSIKA inputs/outputs
├── results/                   # Volume mount for analysis results
└── activate_panoseti.sh       # Environment activation script
```

## Makefile Targets

| Target | Description |
|--------|-------------|
| `make build` | Build the Docker image |
| `make test` | Run all container tests |
| `make run` / `make shell` | Run interactive container |
| `make clean` | Remove Docker image |
| `make example` | Run simulation + visualization (X11 required) |
| `make example-sim` | Run simulation only (batch mode) |
| `make example-vis` | Visualize existing example.root (X11 required) |
| `make help` | Show help message |

## Environment Variables

The container can be configured with the following environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| `CORSIKA_USER` | CORSIKA download username | (none) |
| `CORSIKA_PASSWORD` | CORSIKA download password | (none) |
| `IMAGE_NAME` | Docker image name | `panoseti` |
| `IMAGE_TAG` | Docker image tag | `latest` |

## CORSIKA Credentials

**Security Notice:** CORSIKA credentials are NOT stored in this repository.

- Credentials are only passed during `docker build` as build arguments
- They are not stored in any files committed to Git
- The final Docker image does not contain the credentials (only the downloaded software)
- Recommended: Build with credentials using `--build-arg`, then delete the shell history

To get CORSIKA credentials:
1. Register at https://www.iap.kit.edu/corsika/
2. Receive username and password via email
3. Use during Docker build with `--build-arg CORSIKA_USER=... --build-arg CORSIKA_PASSWORD=...`

**Alternative without credentials:**
Pre-download `corsika-77550.tar.gz` to the container directory before building:
```bash
cd container
wget --user=YOUR_USER --password=YOUR_PASS \
  https://web.iap.kit.edu/corsika/download/corsika-v7750/corsika-77550.tar.gz
docker build -t panoseti:latest .
```

## Troubleshooting

### CORSIKA Download Fails

**Solution 1: Use credentials**
```bash
docker build \
  --build-arg CORSIKA_USER=user \
  --build-arg CORSIKA_PASSWORD=pass \
  -t panoseti:latest .
```

**Solution 2: Pre-download**
```bash
cd container
wget --user=YOUR_USER --password=YOUR_PASS \
  https://web.iap.kit.edu/corsika/download/corsika-v7750/corsika-77550.tar.gz
docker build -t panoseti:latest .
```

**Solution 3: Get credentials**
Register at https://www.iap.kit.edu/corsika/ and rebuild

### ROOT JIT Warnings

You may see `libssl.so.1.1` warnings in interactive ROOT. These are non-critical. All PANOSETI operations work correctly:
- Compiling C++ with ROOT
- Running ROOT macros in batch mode
- Using ROOT libraries in your code

### Graphics Display Issues

If ROOT doesn't open windows:

```bash
# Allow X11 connections
xhost +local:docker

# Run with X11 forwarding
docker run -it --rm \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  panoseti:latest
```

### Permission Issues with Volume Mounts

```bash
# Run as current user
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v $(pwd)/data:/opt/panoseti/data \
  panoseti:latest
```

## Performance Notes

- **CORSIKA compilation**: Takes ~30-60 minutes
- **CORSIKA simulation**: ~5-10 seconds per event on typical hardware
- **Image size**: ~4GB after all installations
- **Memory requirement**: 4GB minimum, 8GB recommended

## Testing

Run the full test suite:
```bash
make test
```

This performs:
1. Check ROOT installation
2. Verify CORSIKA executable
3. Test corsikaIOreader
4. Verify example files
5. Run a quick CORSIKA simulation

## Advanced Usage

### Custom CORSIKA Configuration

To use custom CORSIKA build options, modify `container/scripts/install_corsika.sh`:
- Edit the coconut configuration inputs
- Rebuild the container

### Adding Additional ROOT Packages

The container includes a standard ROOT installation. To add packages:

1. Install system dependencies in Dockerfile
2. Recompile or install additional ROOT plugins
3. Rebuild the container

### Batch Processing

Run multiple simulations in parallel inside the container:

```bash
docker run --rm -d \
  -v $(pwd)/batch:/opt/panoseti/data \
  panoseti:latest \
  bash -c "cd /opt/panoseti/corsika-77550/run && for file in /opt/panoseti/data/*.inp; do ./corsika77550Linux_QGSII_urqmd < \$file > \${file%.inp}.log; done"
```

## Container Build Stages

The Dockerfile uses sequential build stages:

1. **Base system**: Install Debian dependencies
2. **CORSIKA**: Download, configure, and compile simulation engine
3. **corsikaIOreader**: Download and compile IO tool
4. **Examples**: Copy example simulation data
5. **Configuration**: Set up environment and paths

## References

- **CORSIKA**: https://www.iap.kit.edu/corsika/
- **ROOT**: https://root.cern/
- **PANOSETI**: https://panoseti.ucsd.edu/
- **Repository**: https://github.com/nkorzoun/PANOSETI-High-Energy-Analysis-Pipeline

## Support

For issues or questions:
- Check troubleshooting sections above
- Review main repository documentation at `AGENTS.md` and `README.md`
- Run `make test` to verify container setup
- Report bugs on GitHub issues
- Contact PANOSETI team directly for project-specific questions