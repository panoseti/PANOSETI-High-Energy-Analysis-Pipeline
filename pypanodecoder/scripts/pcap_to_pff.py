#!/usr/bin/env python3

# pcap_to_pff.py: Convert PANOSETI PCAP/PCAPNG files to PFF.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-06-14)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import sys
import os
import argparse
import glob

# Add directory containing 'pypanodecoder' to sys.path
# __file__ is .../pypanodecoder/scripts/pcap_to_pff.py
# We need .../ in sys.path
package_parent = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if package_parent not in sys.path:
    sys.path.insert(0, package_parent)

from pypanodecoder.eventbuilder import load_pcap_camera_images
from pypanodecoder.data_writer import PffWriter

def main():
    parser = argparse.ArgumentParser(description="Convert PANOSETI PCAP/PCAPNG files to PFF.")
    parser.add_argument('filenames', nargs='+', help="Input PCAP/PCAPNG files or glob patterns.")
    parser.add_argument('--module-id', type=int, help="Module ID for metadata and filename.")
    parser.add_argument('--output', '-o', help="Output filename. If not provided, generated from template.")
    parser.add_argument('--template', help="Filename template for output.")
    parser.add_argument('--scope', help="Override scope name for template (e.g. Fern).")
    parser.add_argument('--min-quabos', type=int, default=1, help="Minimum quabos required to write an event (default 1).")
    parser.add_argument('--use-packet-num', action='store_true', help="Merge packets into events based on packet_num.")
    parser.add_argument('--no-sort', action='store_true', help="Do not sort events chronologically; preserve file order.")
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

    writer_kwargs = {
        'filename': args.output,
        'module_id': args.module_id,
        'scope': args.scope,
    }
    if args.template:
        writer_kwargs['filename_template'] = args.template

    writer = PffWriter(**writer_kwargs)

    try:
        if args.verbose:
            print(f"Loading events from {len(all_files)} files...")
        
        # load_pcap_camera_images handles merging packets into CameraEvent instances
        images = load_pcap_camera_images(
            all_files, 
            store_camera_events=True, 
            min_quabos=args.min_quabos,
            module_id=args.module_id,
            use_packet_num=args.use_packet_num,
            no_sort=args.no_sort,
            verbose=args.verbose
        )

        if images.events:
            if args.verbose:
                print(f"Writing {len(images.events)} events to PFF...")
            
            for event in images.events:
                writer.write_camera_event(event)
            
            if writer.filename:
                print(f"Successfully converted {len(images.events)} events to {writer.filename}")
        else:
            print("No events found to convert.")

    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        writer.close()

if __name__ == "__main__":
    main()
