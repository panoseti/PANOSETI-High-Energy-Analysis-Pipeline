"""
Traces CORSIKA's Cherenkov photon bunches through the PANOSETI optical model
(optical_model.py and resources/optical_model_data_pack.json) on their way to sim_telarray, which
runs with bypass_optics = 1 and so takes each bunch's position as where it lands on the camera.

A CORSIKA IACT stream (eventio) goes in and the same stream comes out, every object the same
length; only the photon bunches (1205) change, for each bunch:
    wavelength  drawn here per bunch from the Cherenkov spectrum (1/lambda^2) over CORSIKA's
                wavelength range (event header words 95, 96), as sim_telarray would otherwise
                draw it, so the tracing, the lens transmission and sim_telarray's QE and
                atmosphere all use it
    x, y        where the bunch lands on the camera, through the lens at that wavelength, written
                so that sim_telarray's bypass ends up there (see camera_to_bunch())
    photons     scaled by the lens's (PMMA) transmission at that wavelength, 0 if the lens
                vignettes the bunch
    cx, cy, time and emission height are kept: sim_telarray's airmass uses the real direction.

Camera coordinates are sim_telarray's: a photon from (x, y) [rad] off the pointing, +x towards
higher altitude and +y towards lower azimuth (heap.process_simulation.camera_position()), lands
near (F x, F y) [cm].

Usage:
    python simulation_tools/trace_photons.py --theta 20 --phi 180 --seed 1 IN.corsika.zst OUT
with IN or OUT '-' for standard input/output (IN is zstd-decompressed if it ends in .zst). The
pointing (--theta, --phi) is sim_telarray's: zenith angle, azimuth from north, clockwise.
"""
import argparse
import struct
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import optical_model # noqa: E402

SYNC = b"\x37\x8a\x1f\xd4" # eventio sync marker, little-endian
TELESCOPE_DEFINITION, EVENT_HEADER, TELESCOPE_DATA, PHOTONS = 1201, 1202, 1204, 1205
BUNCH_COLUMNS = 8 # x, y, cx, cy, time, zem, photons, wavelength (float32)
HC = 1239.841984 # eV nm
# Rays start this far [cm] back along their direction from CORSIKA's bunch position, so they all
# start in front of the lens (optical_model only propagates forward).
RAY_START = 200.
DONE = "trace_photons done"


class PhotonTracer:
    """
    CORSIKA photon bunches -> camera positions through optical_model, for sim_telarray's
    bypass_optics.

    Parameters:
        theta, phi: telescope pointing [deg], sim_telarray's convention (zenith angle; azimuth
            from north, clockwise)
        thick_lens: trace the lens's grooves (optical_model.trace_telescope_thick) instead of
            the thin lens with the ideal aspheric profile (trace_telescope_thin)
        focal_offset: [cm] camera plane at the lens's focal length + focal_offset
        seed: for the wavelengths
        fixed_n: one refractive index for every wavelength instead of PMMA's n(lambda), for
            comparing with sim_telarray's own (achromatic) lens
    """
    def __init__(self, theta, phi, thick_lens=False, focal_offset=0.0, seed=None, fixed_n=None):
        datapack = optical_model.load_datapack()
        self.optics = datapack["optical_model"]
        self.focal_length = self.optics["F"]
        self.trace = optical_model.trace_telescope_thick if thick_lens else optical_model.trace_telescope_thin
        self.focal_offset = focal_offset
        if fixed_n is None:
            self.n = optical_model.create_n_interpolator(datapack)
        else:
            self.n = lambda energy: np.full(np.shape(energy), fixed_n, dtype=float)
        pmma = datapack["pmma_transmission"]
        order = np.argsort(pmma["ev"])
        self.transmission_ev = np.asarray(pmma["ev"], dtype=float)[order]
        self.transmission = np.asarray(pmma["eff"], dtype=float)[order]
        self.rng = np.random.default_rng(seed)
        np.random.seed(None if seed is None else seed % 2**32) # optical_model's own draws (scattering, reflections) if enabled

        # Ground frame: CORSIKA's, x north, y west, z up. Lens frame (optical_model): optical axis
        # along y with light travelling towards -y, camera plane (x, z). Rows: lens axes in ground
        # coordinates, chosen so that a photon from towards higher altitude (+e_alt) lands at +x
        # and one from towards lower azimuth (+e_az) at +z. (A reflection, which the lens's
        # rotational symmetry doesn't see.)
        alt, az = np.radians(90. - theta), np.radians(phi)
        axis = np.array([np.cos(alt) * np.cos(az), -np.cos(alt) * np.sin(az), np.sin(alt)]) # towards the sky
        e_alt = np.array([-np.sin(alt) * np.cos(az), np.sin(alt) * np.sin(az), np.cos(alt)])
        e_az = np.array([np.sin(az), np.cos(az), 0.]) # towards lower azimuth
        self.to_lens = np.array([-e_alt, axis, -e_az])

        self.stats = {"bunches": 0, "vignetted": 0, "photons_in": 0., "photons_out": 0.}

    def wavelengths(self, n, wavelength_range):
        """n wavelengths [nm] from Cherenkov's 1/lambda^2 over wavelength_range, as sim_telarray draws them."""
        low, high = wavelength_range
        return 1. / (1. / low - self.rng.random(n) * (1. / low - 1. / high))

    def camera_positions(self, bunches, wavelength):
        """
        Where bunches (N, 8) land on the camera [cm], through the lens at wavelength [nm].

        Returns:
            (x, y, valid): sim_telarray camera coordinates, and False for bunches the lens vignettes
        """
        cx, cy = bunches[:, 2].astype(float), bunches[:, 3].astype(float)
        direction = np.stack([cx, cy, -np.sqrt(np.clip(1. - cx**2 - cy**2, 0., None))], axis=1)
        position = np.stack([bunches[:, 0], bunches[:, 1], np.zeros(len(bunches))], axis=1) # cm, telescope center at 0
        rays = optical_model.RayBundle((position - RAY_START * direction) @ self.to_lens.T, direction @ self.to_lens.T,
                                       energy_eV=HC / wavelength)
        rays = self.trace(rays, self.n, self.optics, focal_offset=self.focal_offset)
        return rays.pos[:, 0], rays.pos[:, 2], rays.valid & np.isfinite(rays.pos).all(axis=1)

    def __call__(self, bunches, z_telescope, wavelength_range):
        """
        Traced copy of one telescope's bunches (N, 8) float32, z_telescope [cm] its height above
        the observation level (CORSIKA's telescope definition).
        """
        out = bunches.copy()
        if not len(bunches):
            return out
        wavelength = self.wavelengths(len(bunches), wavelength_range)
        x, y, valid = self.camera_positions(bunches, wavelength)
        bunch_x, bunch_y = camera_to_bunch(x, y, bunches[:, 2], bunches[:, 3], z_telescope)
        out[valid, 0] = bunch_x[valid]
        out[valid, 1] = bunch_y[valid]
        transmission = np.interp(HC / wavelength, self.transmission_ev, self.transmission, left=0., right=0.)
        out[:, 6] = bunches[:, 6] * np.where(valid, transmission, 0.)
        out[:, 7] = wavelength

        self.stats["bunches"] += len(bunches)
        self.stats["vignetted"] += int((~valid).sum())
        self.stats["photons_in"] += float(np.abs(bunches[:, 6]).sum())
        self.stats["photons_out"] += float(np.abs(out[:, 6]).sum())
        return out


def camera_to_bunch(x, y, cx, cy, z_telescope):
    """
    Bunch (x, y) [cm] that sim_telarray's bypass_optics puts at camera (x, y) [cm]: it moves a
    bunch from the telescope's level along its direction to the observation level, z_telescope
    below, and takes (x, -y) there as the camera position (camera rotation 0).
    """
    cx, cy = np.asarray(cx, dtype=float), np.asarray(cy, dtype=float)
    distance = z_telescope / np.sqrt(np.clip(1. - cx**2 - cy**2, 0., None)) # z_telescope * airmass
    return x - distance * cx, -y - distance * cy


def _read(stream, n):
    data = stream.read(n)
    while len(data) < n:
        more = stream.read(n - len(data))
        if not more:
            raise EOFError(f"stream ended inside an eventio object ({len(data)} of {n} bytes)")
        data += more
    return data


def _header(data, offset=0):
    """eventio object header at data[offset:] -> (type, version, length, header size)."""
    type_word, _, length_word = struct.unpack_from("<IiI", data, offset)
    length = length_word & 0x3FFFFFFF
    size = 12
    if type_word & (1 << 17): # extended: 12 more length bits follow
        length |= (struct.unpack_from("<I", data, offset + 12)[0] & 0xFFF) << 30
        size = 16
    return type_word & 0xFFFF, (type_word >> 20) & 0xFFF, length, size


def _patch_photons(content, version, tracer, z_telescopes, wavelength_range):
    array, telescope, _, n_bunches = struct.unpack_from("<hhfi", content)
    if array == 999 and telescope == 999: # IACTEXT particles at the observation level
        return content
    if version // 1000 == 1:
        raise ValueError("compact photon bunches (CORSIKA IACT 'compact' output) are not supported")
    if len(content) != 12 + 4 * BUNCH_COLUMNS * n_bunches:
        raise ValueError(f"photon bunch object of {len(content)} bytes for {n_bunches} bunches")
    if z_telescopes is None or wavelength_range is None:
        raise ValueError("photon bunches before the telescope definition or event header")
    bunches = np.frombuffer(content, dtype="<f4", offset=12).reshape(n_bunches, BUNCH_COLUMNS)
    traced = tracer(bunches, z_telescopes[telescope], wavelength_range).astype("<f4")
    # The header's photon sum, recomputed only if the weights changed (summing in another order
    # than CORSIKA's can change its last bit).
    header = content[:12]
    if not np.array_equal(traced[:, 6], bunches[:, 6]):
        header = struct.pack("<hhfi", array, telescope, float(np.abs(traced[:, 6]).sum()), n_bunches)
    return header + traced.tobytes()


def patch_stream(fin, fout, tracer):
    """
    Copy a CORSIKA IACT eventio stream from fin to fout, tracing its photon bunches with tracer.
    """
    z_telescopes = wavelength_range = None
    while True:
        sync = fin.read(4)
        if not sync:
            break
        if len(sync) < 4 or sync != SYNC:
            raise ValueError(f"expected an eventio sync marker, got {sync!r}")
        head = _read(fin, 12)
        if struct.unpack_from("<I", head)[0] & (1 << 17):
            head += _read(fin, 4)
        object_type, version, length, _ = _header(head)
        content = _read(fin, length)

        if object_type == TELESCOPE_DEFINITION:
            n = struct.unpack_from("<i", content)[0]
            z_telescopes = np.frombuffer(content, dtype="<f4", count=4 * n, offset=4).reshape(4, n)[2].astype(float)
        elif object_type == EVENT_HEADER:
            n = struct.unpack_from("<i", content)[0]
            event_header = np.frombuffer(content, dtype="<f4", count=n, offset=4)
            wavelength_range = (float(event_header[95]), float(event_header[96]))
        elif object_type == TELESCOPE_DATA:
            patched = bytearray()
            offset = 0
            while offset < len(content):
                sub_type, sub_version, sub_length, sub_size = _header(content, offset)
                sub_content = content[offset + sub_size:offset + sub_size + sub_length]
                if sub_type == PHOTONS:
                    sub_content = _patch_photons(sub_content, sub_version, tracer, z_telescopes, wavelength_range)
                patched += content[offset:offset + sub_size] + sub_content
                offset += sub_size + sub_length
            content = bytes(patched)
        fout.write(SYNC + head + content)
    fout.flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="CORSIKA IACT file (.zst is decompressed), or - for standard input")
    parser.add_argument("output", help="output file, or - for standard output")
    parser.add_argument("--theta", type=float, required=True, help="pointing zenith angle [deg]")
    parser.add_argument("--phi", type=float, required=True, help="pointing azimuth [deg], from north, clockwise")
    parser.add_argument("--thick-lens", action="store_true", help="trace the lens's grooves (thick lens)")
    parser.add_argument("--focal-offset", type=float, default=0.0, help="[cm] camera plane at F + focal_offset")
    parser.add_argument("--seed", type=int, default=None, help="for the wavelengths")
    parser.add_argument("--fixed-n", type=float, default=None,
                        help="one refractive index for all wavelengths (testing: compare with sim_telarray's own lens)")
    args = parser.parse_args()

    tracer = PhotonTracer(args.theta, args.phi, args.thick_lens, args.focal_offset, args.seed, args.fixed_n)
    fin = sys.stdin.buffer if args.input == "-" else open(args.input, "rb")
    if args.input.endswith(".zst"):
        import zstandard
        fin = zstandard.ZstdDecompressor().stream_reader(fin)
    fout = sys.stdout.buffer if args.output == "-" else open(args.output, "wb")
    print(f"trace_photons: pointing theta={args.theta:g}, phi={args.phi:g} deg, "
          f"{'thick' if args.thick_lens else 'thin'} lens, focal offset {args.focal_offset:g} cm"
          + (f", fixed n = {args.fixed_n:g}" if args.fixed_n is not None else ""), file=sys.stderr, flush=True)
    patch_stream(fin, fout, tracer)
    stats = tracer.stats
    print(f"{DONE}: {stats['bunches']} bunches, {stats['vignetted']} vignetted, photons {stats['photons_in']:.6g} -> "
          f"{stats['photons_out']:.6g} after the lens's transmission", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
