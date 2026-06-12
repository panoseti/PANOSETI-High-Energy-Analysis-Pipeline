# pypanodecoder

A collection of Python tools for decoding and analyzing PANOSETI data stored in PCAP/PCAPNG and PFF formats.

Author: Stephen Fegan <sfegan@llr.in2p3.fr>
Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

## Overview

This directory contains the pure Python tools for loading and analyzing PANOSETI data. At the low level a raw packet decoding and event merging tool is provided. At a higher level there are utilities for pedestal calculation and data quality monitoring. The tools support both raw PCAP/PCAPNG network captures and the PANOSETI File Format (PFF). The tools aspire to be simple and extensible, so that they can be used to build more complex analysis tasks or as part of a larger pipelines. They are largely free from non-standard dependencies (except NumPy, SciPy, and Matplotlib).

## Modules

### `pcapdecoder.py`
The core decoder for PANOSETI Science packets from PCAP/PCAPNG files.
- **`PcapDecoder`**: A pure Python class that reads `.pcap` and `.pcapng` files without requiring external libraries like `libpcap` or `scapy`.
- **`GTIFilter`**: Manages Good Time Intervals (GTI) for filtering packets into runs during decoding.
- **`get_panoseti_packets`**: A generator function to stream packets from multiple files or glob patterns, optionally filtering based on GTIs.

### `eventbuilder.py`
Tools for reconstructing full camera events from individual Quabo packets.
- **`CameraEvent`**: Groups packets from the four Quabos (quadrants) of a telescope that share the same hardware timestamp.
- **Image Reconstruction**: Provides a `get_image()` method that rotates and mosaics Quabo data into a 32x32 camera image.
- **`CameraImages`**: A container for holding stacks of camera images and associated metadata in memory, supporting GTI filtering and pedestal correction. Now supports tracking the data source (PCAP or PFF) for each GTI individually.
- **`get_camera_events`**: A generator function to stream merged camera events from multiple files, with GTI filtering.
- **`load_camera_images`**: Unified high-level utility that handles both PCAP and PFF files, supporting priority-based merging of overlapping GTIs.
- **`load_pcap_camera_images`**: Specific utility to load data directly from PCAP/PCAPNG files.
- **`load_pff_camera_images`**: Specific utility to load data directly from PFF files, including automatic image alignment.

### `pedestals.py`
Utilities for charge spectra analysis and pedestal management.
- **`ChargeHistogram`**: Calculates robust statistics (mean, median, quantiles, Winsorized mean/variance) from charge distributions.
- **`ChargeSpectra`**: Generates histograms at the pixel, SiPM (8x8), Quabo (16x16), and full camera (32x32) levels.
- **Pedestal Correction**: Implements time-dependent polynomial pedestal fitting and subtraction to handle drift in detector baselines.

### `plotting.py`
Standard plotting functions for PANOSETI data.
- **`plot_image`**: Core function for plotting 32x32 (pixel), 4x4 (SiPM), or 2x2 (Quabo) images with automatic component delineation.
- **`plot_star_field`**: Plots a star field around a target or specific coordinates.
- **`overlay_stars`**: Overlays star markers and labels on an image axes using a pointing solution.

### `pointing.py`
Calculate pointing solutions from reference stars and manage target coordinates.
- **`PointingSolution`**: Translates between sky (RA/Dec) and image (x, y) coordinates. Can be solved from two reference stars.
- **`get_target_coordinates`**: Utility to retrieve coordinates for common astronomical targets (e.g., Crab, Mrk 421).

### `stars.py`
Utilities to query and filter the Yale Bright Star Catalog (YBSC5).
- **`get_bright_stars`**: Fetches stars near a given RA/Dec within a specified radius and magnitude limit.
- **Automatic Catalog Fetching**: Automatically downloads and caches the YBSC5 catalog if not found locally.

### `dqm/`
Data Quality Monitoring (DQM) subpackage and visualization tools.
- **`plot_event_rate`**: Plots the event rate (Hz) over time, with support for multiple GTIs, subplots, and absolute UT time.
- **`plot_delta_t`**: Analyzes the distribution of time intervals between consecutive events. Includes exponential fitting to measure random trigger rates.

## Example: decode Quabo packets from one file

```python
from pypanodecoder.pcapdecoder import get_panoseti_packets
# Stream packets from a single PCAP file
for packet in get_panoseti_packets("data/20260110/onsky_20260110_050000.pcapng"):
    print(f"Quabo: {packet.quabo_id}, timestamp: {packet.event_time}")
```

This example demonstrates how to use the `get_panoseti_packets` generator to read packets from a single PCAP file. Each packet is decoded into a structured format `SciencePacket` (described in `pcapdecoder.py`), with attributes like `quabo_id`, `event_time` and `pix_data`. Note that multiple files can be read using glob patterns, and GTI filtering can be applied to select only packets within specified time intervals (see below).

## Example: merge Quabo packets from multiple files into camera events

```python
from pypanodecoder.eventbuilder import get_camera_events
GTIs = [ 
    { "start": "2026-01-10 05:00:00", "end": "2026-01-10 08:00:00" },
    { "start": "2026-01-11 05:30:00", "end": "2026-01-11 07:45:00" }
]
# Stream camera events from multiple PCAP files
for camera_event in get_camera_events("data/202601*/onsky_*.pcapng", gtis=GTIs):
    print(f"Camera event at {camera_event.event_time} with {len(camera_event.packets)} packets")
```
Here we demonstrate how to use the `get_camera_events` generator to read and merge Quabo packets into full camera events. The function takes a glob pattern to read multiple files  and an optional list of GTIs to filter the events. Here we illustrate how to specify the whole January 2026 dataset (assuming a heirarchical organization by date) and select only events that occur during two hypothetical  runs, each specified by a start and end time.

Each call to `get_camera_events` yields a `CameraEvent` containing the merged data from the four Quabos that correspond to the same hardware timestamp, allowing for easy manipulation of the full camera image.

## Example: load camera images and plot DQM metrics

```python
from pypanodecoder.eventbuilder import load_camera_images
from pypanodecoder.dqm import plot_event_rate, plot_delta_t
import matplotlib.pyplot as plt

# Define Good Time Intervals (GTIs)
GTIs = [ 
    { "start": "2026-01-10 05:00:00", "end": "2026-01-10 08:00:00" },
    { "start": "2026-01-11 05:30:00", "end": "2026-01-11 07:45:00" }
]

# Load images from multiple files (PCAP or PFF) with GTI filtering
data = load_camera_images("data/202601*/*", gtis=GTIs)

# Check the data source
print(f"Overall source: {data.data_source}") # PCAP, PFF, or MIXED
print(f"Source per GTI: {data.source}")      # {0: 'PCAP', 1: 'PFF', ...}

# Plot the event rate
fig_rate, _ = plot_event_rate(data, subplots=True, uttime=True)
plt.show()

# Plot the inter-event time distribution
fig_dt, _ = plot_delta_t(data, fit=True)
plt.show()
```

Here we show how to load the camera images from multiple runs into memory and then process them. This is the preferred way to work with the data. The `load_camera_images` function reads the specified files (automatically detecting PCAP or PFF format), applies GTI filtering to extract the on-source data, and returns a `CameraImages` object containing the image data and metadata. Once the data is loaded, it can easily be accessed to retrieve the event data or process it. The `data.source` attribute tracks the origin of each GTI.

In this example we use the DQM plotting functions to plot the event rate histogram and inter-event (Delta-T) time distribution. The `plot_event_rate` is configured to create subplots for each GTI and display time in absolute UT time. The `plot_delta_t` function fits an exponential model to the inter-event time distribution to estimate the random trigger rate.

## Dependencies
- **`pcapdecoder.py`**: No external dependencies (Standard Library only).
- **Others**: `numpy`, `matplotlib`, `scipy`.
