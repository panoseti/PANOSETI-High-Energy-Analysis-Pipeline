# pypanodecoder

A collection of Python tools for decoding and analyzing PANOSETI data stored in PCAP/PCAPNG formats.

Author: Stephen Fegan <sfegan@llr.in2p3.fr>
Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

## Overview

This directory contains a pure Python implementation of the PANOSETI data pipeline, from raw packet decoding to higher-level calibration and data quality monitoring. These tools aspire to be easy to use and free of non-standard dependencies (except NumPy, SciPy, and Matplotlib).

## Modules

### `panodecoder.py`
The core decoder for PANOSETI Science packets.
- **`PanosetiPcapDecoder`**: A pure Python class that reads `.pcap` and `.pcapng` files without requiring external libraries like `libpcap` or `scapy`.
- **`GTIFilter`**: Efficiently manages Good Time Intervals (GTI) for filtering events during decoding.
- **Time Synchronization**: Implements `wr_to_unix` to convert White Rabbit (WR) hardware timestamps (TAI/nanosec) to Unix epoch time.
- **`get_panoseti_events`**: A generator function to stream events from multiple files or glob patterns.

### `eventbuilder.py`
Tools for reconstructing full camera events from individual Quabo packets.
- **`PanosetiCameraEvent`**: Groups packets from the four Quabos (quadrants) of a telescope that share the same hardware timestamp.
- **Image Reconstruction**: Provides a `get_image()` method that correctly rotates and mosaics Quabo data into a 32x32 camera image.
- **`PanosetiCameraImages`**: A container for stacks of camera images and associated metadata, supporting GTI filtering and pedestal correction.
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

## Usage Example

```python
from pypanodecoder.eventbuilder import load_camera_images
from pypanodecoder.dqm import plot_event_rate, plot_delta_t
import matplotlib.pyplot as plt

# Load images from a PCAP file
data = load_camera_images("data/run_001.pcap")

# Plot the event rate
fig_rate, _ = plot_event_rate(data)
plt.show()

# Plot the inter-event time distribution
fig_dt, _ = plot_delta_t(data, fit=True)
plt.show()
```

## Dependencies
- **`panodecoder.py`**: No external dependencies (Standard Library only).
- **Others**: `numpy`, `matplotlib`, `scipy`.
