#!/usr/bin/env python3
"""
Quick visualizer for simplified binary files produced by ParticleAndHitHostBinaryDumper.

It can show:
  - hits positions (x,y) or (r,z)
  - tracks (from TrackingParticles) overlaid

Usage:
  python plot-binary --hits hits.bin [--particles particles.bin] [--map map.bin] [--view xy|rz]
"""

import struct
import argparse
import numpy as np
import matplotlib.pyplot as plt

def read_header(f, expected_magic: bytes):
    """Check header and return number of events."""
    magic = f.read(4)
    if magic != expected_magic:
        raise ValueError(f"Invalid magic header: got {magic}, expected {expected_magic}")
    version, endian, nEvents = struct.unpack("III", f.read(12))
    if endian != 0x01020304:
        raise ValueError("Endianness mismatch")
    return nEvents


def read_particles(fname):
    """Read ParticleSimpleHost from file."""
    with open(fname, "rb") as f:
        nEvents = read_header(f, b"PAR1")
        particles = []
        for _ in range(nEvents):
            (nPart,) = struct.unpack("I", f.read(4))
            if nPart == 0:
                particles.append({})
                continue
            # each block contains 14 arrays
            cols = {}
            float_arrays = ["vx","vy","vz","px","py","pz","energy","pt","eta","phi","mass"]
            int16_arrays = ["charge"]
            int32_arrays = ["pdgID"]
            uint32_arrays = ["partInd"]

            def read_array(dtype, n):
                return np.frombuffer(f.read(n * np.dtype(dtype).itemsize), dtype=dtype)

            for k in float_arrays:
                cols[k] = read_array(np.float32, nPart)
            for k in int16_arrays:
                cols[k] = read_array(np.int16, nPart)
            for k in int32_arrays:
                cols[k] = read_array(np.int32, nPart)
            for k in uint32_arrays:
                cols[k] = read_array(np.uint32, nPart)
            particles.append(cols)
    return particles


def read_hits(fname):
    """Read TrackingRecHitHost from file."""
    with open(fname, "rb") as f:
        nEvents = read_header(f, b"TRH1")
        hits = []
        for _ in range(nEvents):
            nHits, nModules = struct.unpack("II", f.read(8))
            _ = f.read((nModules + 1) * 4)  # skip moduleStart
            cols = {}
            float_arrays = [
                "xLocal","yLocal","xerrLocal","yerrLocal",
                "xGlobal","yGlobal","zGlobal","rGlobal"
            ]
            int16_arrays = ["iphi","chargeAndStatus","clusterSizeX","clusterSizeY"]
            uint32_arrays = ["detectorIndex"]
            def read_array(dtype, n):
                return np.frombuffer(f.read(n * np.dtype(dtype).itemsize), dtype=dtype)
            for k in float_arrays:
                cols[k] = read_array(np.float32, nHits)
            for k in int16_arrays:
                cols[k] = read_array(np.int16, nHits)
            for k in uint32_arrays:
                cols[k] = read_array(np.uint32, nHits)
            hits.append(cols)
    return hits


def read_map(fname):
    """Read SimpleMapHost from file."""
    with open(fname, "rb") as f:
        nEvents = read_header(f, b"MAP1")
        maps = []
        for _ in range(nEvents):
            (nEntries,) = struct.unpack("I", f.read(4))
            keys = np.frombuffer(f.read(4*nEntries), dtype=np.uint32)
            vals = np.frombuffer(f.read(4*nEntries), dtype=np.uint32)
            maps.append({"key": keys, "value": vals})
    return maps

def plot_event(hits, particles = None, maps=None, view="xy"):
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.set_aspect("equal")
    if view == "xy":
        x = hits["xGlobal"]
        y = hits["yGlobal"]
        ax.set_xlabel("x [cm]")
        ax.set_ylabel("y [cm]")
        ax.set_title("Hits in XY view")
        ax.scatter(x, y, s=4, color="gray", alpha=0.6, label="Hits")
    else:  # rz view
        r = hits["rGlobal"]
        z = hits["zGlobal"]
        ax.set_xlabel("z [cm]")
        ax.set_ylabel("r [cm]")
        ax.set_title("Hits in RZ view")
        ax.scatter(z, r, s=4, color="gray", alpha=0.6, label="Hits")

    # overlay particle tracks (approx. straight lines)
    if particles:
        p = particles
        colors = plt.cm.tab10(np.linspace(0, 1, min(10, len(p["pt"]))))
        for i, (eta, phi, pt) in enumerate(zip(p["eta"], p["phi"], p["pt"])):
            color = colors[i % len(colors)]
            if view == "xy":
                # direction from phi
                vx, vy = p["vx"][i], p["vy"][i]
                px, py = p["px"][i], p["py"][i]
                ax.arrow(vx, vy, px*0.05, py*0.05, color=color, head_width=0.2)
            else:
                z0, r0 = p["vz"][i], np.hypot(p["vx"][i], p["vy"][i])
                pz = p["pz"][i]
                pr = np.hypot(p["px"][i], p["py"][i])
                ax.arrow(z0, r0, pz*0.05, pr*0.05, color=color, head_width=0.2)
    if maps:
        ax.set_title(ax.get_title() + " (truth-colored)")
    ax.legend()
    plt.show()
    plt.savefig("bin_plots.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot hits and tracks from binary SoA files")
    parser.add_argument("--hits", required=True)
    parser.add_argument("--particles", default=None)
    parser.add_argument("--map", default=None)
    parser.add_argument("--event", type=int, default=0)
    parser.add_argument("--view", choices=["xy", "rz"], default="xy")
    args = parser.parse_args()

    print(f"Reading hits from {args.hits}")
    hits_all = read_hits(args.hits)
    
    if args.particles:
        print(f"Reading particles from {args.particles}")
        parts_all = read_particles(args.particles)
    maps_all = None
    if args.map:
        print(f"Reading map from {args.map}")
        maps_all = read_map(args.map)

    ev = args.event
    if ev >= len(hits_all):
        raise ValueError(f"Event {ev} out of range.")

    hits = hits_all[ev]

    particles = None
    if args.particles:
        particles = parts_all[ev]

    maps = maps_all[ev] if maps_all else None

    plot_event(hits, particles, maps, view=args.view)
