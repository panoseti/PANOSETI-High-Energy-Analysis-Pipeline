#!/usr/bin/env python3

# stars.py: Functions to query the Yale Bright Star Catalog for bright stars near a given position.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-28)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import math
import urllib.request
import gzip
import io

YBSC5_URLS = [
    "http://tdc-www.harvard.edu/catalogs/bsc5.dat.gz",
    "https://github.com/sfegan/PANOSETI-High-Energy-Analysis-Pipeline/raw/refs/heads/main/pypanodecoder/resources/ybsc5.gz"
]

def _rad(x):
    return math.radians(x)

def _angular_sep(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(_rad, [ra1, dec1, ra2, dec2])
    cosv = (
        math.sin(dec1) * math.sin(dec2)
        + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    )
    cosv = max(-1.0, min(1.0, cosv))
    return math.degrees(math.acos(cosv))

_ybsc_cache = None

def _load_ybsc():
    global _ybsc_cache
    if _ybsc_cache is not None:
        return _ybsc_cache

    raw = None
    for url in YBSC5_URLS:
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Mozilla/5.0"}
            )
            with urllib.request.urlopen(req, timeout=10) as r:
                raw = r.read()
            break
        except Exception:
            continue

    if raw is None:
        raise RuntimeError(f"Failed to download Yale Bright Star Catalog from all sources: {YBSC5_URLS}")

    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
        text = gz.read().decode("latin-1", errors="ignore")

    ybsc = []
    for line in text.splitlines():
        # Lines shorter than 114 chars can't have all requested fields
        if len(line) < 114:
            continue
        try:
            hr   = line[0:4].strip()
            name = line[4:14].strip()

            # --- J2000 RA: bytes 75-76 (H), 77-78 (M), 79-82 (S, F4.1) ---
            rah_s  = line[75:77].strip()
            ram_s  = line[77:79].strip()
            ras_s  = line[79:83].strip()

            # --- J2000 Dec: byte 83 (sign), 84-85 (D), 86-87 (M), 88-89 (S) ---
            de_sign = line[83]
            ded_s   = line[84:86].strip()
            dem_s   = line[86:88].strip()
            des_s   = line[88:90].strip()

            # --- Vmag: bytes 102-106 (F5.2) ---
            vmag_s = line[102:107].strip()

            # --- B-V: bytes 110-114 (F5.2) ---
            bv_s = line[109:114].strip()

            if not rah_s or not ded_s or not vmag_s:
                continue

            rah  = float(rah_s)
            ram  = float(ram_s) if ram_s else 0.0
            ras  = float(ras_s) if ras_s else 0.0
            ra   = 15.0 * (rah + ram / 60.0 + ras / 3600.0)

            ded  = float(ded_s)
            dem  = float(dem_s) if dem_s else 0.0
            des  = float(des_s) if des_s else 0.0
            dec  = ded + dem / 60.0 + des / 3600.0
            if de_sign == '-':
                dec = -dec

            vmag = float(vmag_s)
            bv   = float(bv_s) if bv_s else None

            ybsc.append({
                "hr":     hr,
                "name":   name,
                "ra_deg": ra,
                "dec_deg": dec,
                "vmag":   vmag,
                "bv":     bv,
            })

        except (ValueError, IndexError):
            continue

    _ybsc_cache = ybsc
    return ybsc

def get_bright_stars(ra_deg, dec_deg, radius_deg=8.0, max_magnitude=8.0):
    ybsc = _load_ybsc()
    results = []
    for star in ybsc:
        sep = _angular_sep(ra_deg, dec_deg, star["ra_deg"], star["dec_deg"])
        if sep <= radius_deg and star["vmag"] <= max_magnitude:
            results.append({
                "hr":             star["hr"],
                "name":           star["name"],
                "ra_deg":         star["ra_deg"],
                "dec_deg":        star["dec_deg"],
                "vmag":           star["vmag"],
                "bv":             star["bv"],
                "separation_deg": sep,
            })
    results.sort(key=lambda x: x["vmag"])
    return results

def print_stars(stars):
    """
    Prints a formatted table of stars.

    Args:
        stars (list): List of star dictionaries from get_bright_stars().
    """
    if not stars:
        print("No stars found.")
        return

    # Table header
    header = f"{'Idx':>3} {'HR':>4} {'Name':<12} {'RA [deg]':>10} {'Dec [deg]':>10} {'Vmag':>6} {'B-V':>6} {'Sep [deg]':>10}"
    print(header)
    print("-" * len(header))

    for i, star in enumerate(stars):
        hr = star.get('hr', '')
        # Merge multiple whitespaces in name for cleaner printing
        name = " ".join(star.get('name', '').split())
        ra = star.get('ra_deg', 0.0)
        dec = star.get('dec_deg', 0.0)
        vmag = star.get('vmag', 0.0)
        bv = star.get('bv')
        bv_str = f"{bv:6.2f}" if bv is not None else f"{'':>6}"
        sep = star.get('separation_deg', 0.0)

        print(f"{i:3d} {hr:>4} {name:<12} {ra:10.4f} {dec:10.4f} {vmag:6.2f} {bv_str} {sep:10.4f}")