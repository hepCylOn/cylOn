#!/usr/bin/env python3
import numpy as np
import json
from pathlib import Path
import math
import argparse

mappingjson = "inputhits.json"

def phi2short(phi):
    p2i = (np.iinfo(np.int16).max + 1) / math.pi
    return np.round(phi * p2i).astype(np.int16)

def generate_hits(n_hits, n_modules):
    # module IDs spread between 0 and n_modules-1
    detector_index = np.random.randint(0, n_modules, size=n_hits, dtype=np.int32)
    # Local positions (cm)
    x_local = np.random.normal(0, 0.05, size=n_hits).astype(np.float32)
    y_local = np.random.normal(0, 0.05, size=n_hits).astype(np.float32)
    xerr_local = np.abs(np.random.normal(0.01, 0.002, size=n_hits)).astype(np.float32)
    yerr_local = np.abs(np.random.normal(0.01, 0.002, size=n_hits)).astype(np.float32)

    # Global positions (cm)
    phi = np.random.uniform(-math.pi, math.pi, size=n_hits)
    r = np.random.uniform(20, 100, size=n_hits)
    z = np.random.uniform(-300, 300, size=n_hits)

    x_global = (r * np.cos(phi)).astype(np.float32)
    y_global = (r * np.sin(phi)).astype(np.float32)
    z_global = z.astype(np.float32)
    r_global = r.astype(np.float32)
    iphi = phi2short(phi)

    # Cluster & charge info
    charge = np.random.randint(1000, 20000, size=n_hits, dtype=np.int32)
    cluster_size_x = np.random.randint(1, 4, size=n_hits, dtype=np.int16)
    cluster_size_y = np.random.randint(1, 4, size=n_hits, dtype=np.int16)

    # Compose line strings (comma separated)
    cols = np.stack([
        x_local, y_local, xerr_local, yerr_local,
        x_global, y_global, z_global, r_global,
        iphi, charge,
        cluster_size_x, cluster_size_y, detector_index
    ], axis=1)

    lines = [",".join(str(v) for v in row) for row in cols]
    return lines

def main():

    parser = argparse.ArgumentParser(description="Generate dummy hit text files for SoA loader testing")
    parser.add_argument("--output", "-o", type=str, default="dummy_hits",
                        help="Output directory for generated files (default: dummy_hits)")
    parser.add_argument("--events", "-n", type=int, default=10,
                        help="Number of events to generate (default: 10)")
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(exist_ok=True)

    mapping = {
        "xLocal": 0,
        "yLocal": 1,
        "xerrLocal": 2,
        "yerrLocal": 3,
        "xGlobal": 4,
        "yGlobal": 5,
        "zGlobal": 6,
        "rGlobal": 7,
        "iphi": 8,
        "charge": 9,
        "clusterSizeX": 10,
        "clusterSizeY": 11,
        "detectorIndex": 12
    }

    with open(outdir / mappingjson, "w") as f:
        json.dump(mapping, f, indent=2)

    for ev in range(args.events):
        n_hits = np.random.randint(100, 500)
        n_modules = np.random.randint(5, 15)
        lines = generate_hits(n_hits, n_modules)
        # Add blank line as event delimiter
        text = "\n".join(lines) + "\n\n"
        with open(outdir / f"event_{ev+1}.txt", "w") as f:
            f.write(text)

    print(f"Generated {args.events} dummy events in {outdir}/")

if __name__ == "__main__":
    main()
