#!/usr/bin/env python3

# __init__.py: pypanodecoder - PANOSETI Science Packet Decoder and Analysis

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-18)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

__version__ = '0.1.0'

from .pcapdecoder import (
    PcapDecoder,
    GTIFilter,
    SciencePacket,
    get_panoseti_packets
)

from .eventbuilder import (
    CameraEvent,
    CameraImages,
    CameraEventBuilder,
    get_camera_events,
    load_camera_images
)

from .pedestals import (
    ChargeHistogram,
    ChargeSpectra,
    calculate_charge_spectra,
    apply_polynomial_pedestal_correction
)

from .plotting import (
    plot_image
)

from . import dqm

__all__ = [
    'PcapDecoder',
    'GTIFilter',
    'SciencePacket',
    'get_panoseti_packets',
    'CameraEvent',
    'CameraImages',
    'CameraEventBuilder',
    'get_camera_events',
    'load_camera_images',
    'ChargeHistogram',
    'ChargeSpectra',
    'calculate_charge_spectra',
    'apply_polynomial_pedestal_correction',
    'plot_image',
    'dqm'
]
