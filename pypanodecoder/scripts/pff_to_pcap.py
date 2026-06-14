#!/usr/bin/env python3

# pff_to_pcap.py: Convert PANOSETI PFF files to PCAPNG.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-06-14)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import sys
import os
import argparse
import glob
import re

# Add directory containing 'pypanodecoder' to sys.path
# __file__ is .../pypanodecoder/scripts/pff_to_pcap.py
# We need .../ in sys.path
package_parent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if package_parent not in sys.path:
    sys.path.insert(0, package_parent)

from pypanodecoder.eventbuilder import load_pff_camera_images
from pypanodecoder.data_writer import PcapWriter

def extract_module_id(filename):
    """Attempts to extract module ID from PFF filename."""
    match = re.search(r'module_(\d+)', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return None

def main():
    parser = argparse.ArgumentParser(description="Convert PANOSETI PFF files to PCAPNG.")
    parser.add_argument('filenames', nargs='+', help="Input PFF files or glob patterns.")
    parser.add_argument('--module-id', type=int, help="Force module ID. If not provided, extracted from filename.")
    parser.add_argument('--output', '-o', help="Output filename. If not provided, generated from template.")
    parser.add_argument('--template', help="Filename template for output.")
    parser.add_argument('--scope', help="Override scope name for template (e.g. Fern).")
    parser.add_argument('--verbose', '-v', action='store_true', help="Verbose output.")

    args = parser.parse_args()

    # Expand glob patterns
    all_files = []
    for pattern in args.filenames:
        expanded = glob.glob(pattern)
        if not expanded and not any(char in pattern for char in '*?['):
            expanded = [pattern]
        all_files.extend(expanded)
    
    all_files = sorted(all_files)
    if not all_files:
        print("Error: No input files found.", file=sys.stderr)
        sys.exit(1)

    # Use first file to guess module_id if not provided
    initial_module_id = args.module_id
    if initial_module_id is None:
        initial_module_id = extract_module_id(all_files[0])
        if initial_module_id is not None and args.verbose:
            print(f"Defaulting to module ID {initial_module_id} (extracted from {os.path.basename(all_files[0])})")

    writer_kwargs = {
        'filename': args.output,
        'module_id': initial_module_id,
        'scope': args.scope,
        'app_name': 'pypanodecoder/pff_to_pcap.py',
        'comment': f"Converted from: {', '.join([os.path.basename(f) for f in all_files])}"
    }
    if args.template:
        writer_kwargs['filename_template'] = args.template

    writer = PcapWriter(**writer_kwargs)

    try:
        total_events = 0
        for filename in all_files:
            if args.verbose:
                print(f"Reading {filename}...")
            
            # Extract module_id for this specific file if not forced
            file_module_id = args.module_id
            if file_module_id is None:
                file_module_id = extract_module_id(filename)
            
            # Load images with events stored
            images = load_pff_camera_images(filename, store_camera_events=True, module_id=file_module_id)
            
            if images.events:
                if args.verbose:
                    print(f"  Writing {len(images.events)} events...")
                for event in images.events:
                    writer.write_camera_event(event)
                total_events += len(images.events)
            elif args.verbose:
                print(f"  No events found in {filename}")

    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        writer.close()
        if writer.filename and writer.packet_count > 0:
            print(f"Successfully converted {total_events} events ({writer.packet_count} packets) to {writer.filename}")
        elif total_events == 0:
            print("No events were converted.")

if __name__ == "__main__":
    main()
