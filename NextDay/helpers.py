from pathlib import Path
import os
import re
import numpy as np
import heap


def _start_key(p: Path):
    # Sort files by timestamp encoded in filename: start_YYYY-...Z
    m = re.search(r"start_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", p.name)
    return m.group(1) if m else p.name


def _timestamps_from_metadata(metadata):
    # Fallback for when pkt-loss cut is skipped
    if isinstance(metadata, dict):
        for k in ("timestamps", "timestamp", "ts"):
            if k in metadata:
                return np.asarray(metadata[k])

    # Structured array / table-like metadata
    if hasattr(metadata, "dtype") and metadata.dtype is not None and metadata.dtype.names:
        for k in ("timestamps", "timestamp", "ts"):
            if k in metadata.dtype.names:
                return np.asarray(metadata[k])

    raise KeyError("Could not find timestamps in metadata when apply_pkt_loss=False.")


def read_module_pffs(
    data_path,
    module,
    apply_pkt_loss=True,
    apply_spike_cut=True,
    require_token="ph1024",
):
    """
    Read all module files in data_path matching:
      start_*{require_token}*.module_{module}.seqno_*.pff

    Then optionally apply:
      - packet-loss cut
      - spike-rate cut

    Returns:
      data_cat: concatenated data array
      ts_cat: concatenated timestamps array
      files_used: list of file paths used
    """
    data_dir = Path(data_path).expanduser()
    pattern = f"start_*{require_token}*.module_{int(module)}.seqno_*.pff"
    files = sorted(data_dir.glob(pattern), key=_start_key)

    if not files:
        raise FileNotFoundError(f"No files matched pattern {pattern} in {data_dir}")

    all_data = []
    all_ts = []

    # Workaround for pypff path parsing issues with dotted directory names:
    # call read_pff on basename after chdir.
    old_cwd = os.getcwd()
    try:
        os.chdir(data_dir)

        for f in files:
            data, metadata = heap.pre_cleaning.read_pff(f.name)

            if apply_pkt_loss:
                data, ts = heap.pre_cleaning.cut_pkt_loss(data, metadata)
            else:
                ts = _timestamps_from_metadata(metadata)

            if apply_spike_cut:
                data, ts = heap.pre_cleaning.spike_cut(data, ts)

            all_data.append(np.asarray(data))
            all_ts.append(np.asarray(ts))

    finally:
        os.chdir(old_cwd)

    data_cat = np.concatenate(all_data, axis=0)
    ts_cat = np.concatenate(all_ts, axis=0)

    return data_cat, ts_cat, [str(p) for p in files]