# pypanodecoder

A collection of Python tools for decoding and analyzing PANOSETI data stored in PCAP/PCAPNG formats.

Author: Stephen Fegan <sfegan@llr.in2p3.fr>
Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

## Overview

This directory contains the pure Python tools for loading and analyzing PANOSETI data. At the low level a raw packet decoding and event merging tool is provided. At a higher level there are utilities for pedestal calculation and data quality monitoring. The tools aspire to be simple and extensible, so that they can be used to build more complex analysis tasks or as part of a larger pipelines. They are largely free from non-standard dependencies (except NumPy, SciPy, and Matplotlib).

## Modules

### `panodecoder.py`
The core decoder for PANOSETI Science packets from PCAP/PCAPNG files.
- **`PanosetiPcapDecoder`**: A pure Python class that reads `.pcap` and `.pcapng` files without requiring external libraries like `libpcap` or `scapy`.
- **`GTIFilter`**: Manages Good Time Intervals (GTI) for filtering packets into runs during decoding.
- **`get_panoseti_packets`**: A generator function to stream packets from multiple files or glob patterns, optionally filtering based on GTIs.

### `eventbuilder.py`
Tools for reconstructing full camera events from individual Quabo packets.
- **`PanosetiCameraEvent`**: Groups packets from the four Quabos (quadrants) of a telescope that share the same hardware timestamp.
- **Image Reconstruction**: Provides a `get_image()` method that rotates and mosaics Quabo data into a 32x32 camera image.
- **`PanosetiCameraImages`**: A container for holding stacks of camera images and associated metadata in memory, supporting GTI filtering and pedestal correction.
- **`load_camera_images`**: High-level utility to load data directly into a 3D NumPy array (32x32xN).

### `pedestals.py`
Utilities for charge spectra analysis and pedestal management.
- **`PanosetiChargeHistogram`**: Calculates robust statistics (mean, median, quantiles, Winsorized mean/variance) from charge distributions.
- **`PanosetiChargeSpectra`**: Generates histograms at the pixel, SiPM (8x8), Quabo (16x16), and full camera (32x32) levels.
- **Pedestal Correction**: Implements time-dependent polynomial pedestal fitting and subtraction to handle drift in detector baselines.

### `dqm.py`
Data Quality Monitoring (DQM) and visualization tools.
- **`plot_event_rate`**: Plots the event rate (Hz) over time, with support for multiple GTIs, subplots, and absolute UT time.
- **`plot_delta_t`**: Analyzes the distribution of time intervals between consecutive events. Includes exponential fitting to measure random trigger rates.

## Usage Example: decode Quabo packets from one file

```python
from pypanodecoder.panodecoder import get_panoseti_packets
# Stream packets from a single PCAP file
for packet in get_panoseti_packets("data/run_001.pcap"):
    print(f"Quabo: {packet.quabo_id}, timestamp: {packet.event_time}")
```

## Usage Example: merge Quabo packets from multiple files into camera events

```python
from pypanodecoder.eventbuilder import get_camera_events
GTIs = [ 
    { "start": "2026-01-10 05:00:00", "end": "2026-01-10 08:00:00" },
    { "start": "2026-01-11 05:30:00", "end": "2026-01-11 07:45:00" }
]
# Stream camera events from multiple PCAP files
for camera_event in get_camera_events("data/run_*.pcap", gtis=GTIs):
    print(f"Camera event at {camera_event.event_time} with {len(camera_event.packets)} packets")
```

## Usage Example: load camera images and plot DQM metrics

```python
from pypanodecoder.eventbuilder import load_camera_images
from pypanodecoder.dqm import plot_event_rate, plot_delta_t
import matplotlib.pyplot as plt

# Define Good Time Intervals (GTIs)
GTIs = [ 
    { "start": "2026-01-10 05:00:00", "end": "2026-01-10 08:00:00" },
    { "start": "2026-01-11 05:30:00", "end": "2026-01-11 07:45:00" }
]

# Load images from multiple PCAP files with GTI filtering
data = load_camera_images("data/run_*.pcap", gtis=GTIs)

# Plot the event rate
fig_rate, _ = plot_event_rate(data, subplots=True, uttime=True)
plt.show()

# Plot the inter-event time distribution
fig_dt, _ = plot_delta_t(data, fit=True)
plt.show()
```

## Dependencies
- **`panodecoder.py`**: No external dependencies (Standard Library only).
- **Others**: `numpy`, `matplotlib`, `scipy`.
