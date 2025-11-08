#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import struct

FORMAT_VERSION = 1
ENDIANNESS = 0x01020304
MAGIC = {
    "particles": b"PAR1",
    "hits": b"TRH1",
    "map": b"MAP1",
}

def write_header(f, tag: str):
    """Write 4B magic + uint32 version + uint32 endian marker"""
    f.write(MAGIC[tag])
    f.write(struct.pack("<I", FORMAT_VERSION))
    f.write(struct.pack("<I", ENDIANNESS))


# ============================================================
# Utility functions
# ============================================================
def compute_kinematics(df):
    """Compute pt, eta, phi, mass, energy."""
    px, py, pz = df["px"], df["py"], df["pz"]
    pt = np.sqrt(px**2 + py**2)
    p = np.sqrt(px**2 + py**2 + pz**2)
    phi = np.arctan2(py, px)
    eta = np.zeros_like(pt)
    mask = (p - pz) > 0
    eta[mask] = 0.5 * np.log((p[mask] + pz[mask]) / (p[mask] - pz[mask]))

    df["pt"], df["eta"], df["phi"] = pt, eta, phi
    df["energy"] = p
    df["mass"] = 0.0


def build_module_index(geom):
    geom = geom.copy()
    geom["module_index"] = np.arange(len(geom), dtype=np.uint32)
    module_map = {
        (int(r.volume_id), int(r.layer_id), int(r.module_id)): int(r.module_index)
        for _, r in geom.iterrows()
    }
    return geom, module_map


def assign_modules(hits, module_map, n_modules):
    hits = hits.copy()
    keys = list(zip(hits.volume_id, hits.layer_id, hits.module_id))
    hits["module_index"] = [module_map[k] for k in keys]
    if hits["module_index"].max() >= n_modules:
        raise RuntimeError("module_index overflow")
    return hits


def build_module_start(hits, n_modules):
    counts = np.bincount(hits["module_index"], minlength=n_modules)
    module_start = np.zeros(n_modules + 1, dtype=np.uint32)
    module_start[1:] = np.cumsum(counts)
    return module_start

def dump_particles(df, out):
    df = df.sort_values("particle_id").reset_index(drop=True)
    n = len(df)
    compute_kinematics(df)

    with open(out, "wb") as f:
        write_header(f, "particles")
        f.write(struct.pack("<I", n))

        # float32 columns
        cols_f32 = [
            "vx","vy","vz",
            "px","py","pz","energy",
            "pt","eta","phi","mass"
        ]
        for c in cols_f32:
            arr = df[c].to_numpy(dtype=np.float32)
            arr.tofile(f)

        # charge
        df["q"].to_numpy(dtype=np.int16).tofile(f)
        # pdgID = 0
        np.zeros(n, np.int32).tofile(f)
        # partInd = 0..n-1
        np.arange(n, dtype=np.uint32).tofile(f)

    print(f"[particles] {out} ({n})")
    return np.arange(n, dtype=np.uint32), df["particle_id"].to_numpy()


def dump_hits(hits, n_modules, out):
    hits = hits.sort_values(["module_index", "hit_id"]).reset_index(drop=True)

    x = hits["x"].to_numpy(dtype=np.float32)
    y = hits["y"].to_numpy(dtype=np.float32)
    z = hits["z"].to_numpy(dtype=np.float32)
    r = np.sqrt(x**2 + y**2).astype(np.float32)
    det = hits["module_index"].to_numpy(dtype=np.uint16)
    n_hits = len(hits)
    module_start = build_module_start(hits, n_modules)

    zeros_f = np.zeros(n_hits, np.float32)
    zeros_i16 = np.zeros(n_hits, np.int16)

    with open(out, "wb") as f:
        write_header(f, "hits")
        f.write(struct.pack("<II", n_hits, n_modules))
        module_start.tofile(f)

        # Write columns in the order expected by the reader
        for _ in range(4): zeros_f.tofile(f)         # xLocal, yLocal, xerrLocal, yerrLocal
        for arr in (x, y, z, r): arr.tofile(f)       # xGlobal, yGlobal, zGlobal, rGlobal
        for _ in range(4): zeros_i16.tofile(f)       # iphi, chargeAndStatus, clusterSizeX, clusterSizeY
        det.tofile(f)                                # detectorIndex

    print(f"[hits] {out} ({n_hits} hits, {n_modules} modules)")
    return hits

def dump_map(truth, hits_sorted, part_ids, part_ind, out):
    pid_to_idx = {int(pid): int(i) for i, pid in zip(part_ind, part_ids)}
    merged = hits_sorted[["hit_id"]].merge(
        truth[["hit_id", "particle_id"]], on="hit_id", how="left"
    )

    ids = np.full(len(merged), 0xFFFFFFFF, np.uint32)
    for i, pid in enumerate(merged["particle_id"].fillna(0).astype(np.int64)):
        if pid in pid_to_idx:
            ids[i] = np.uint32(pid_to_idx[pid])

    with open(out, "wb") as f:
        write_header(f, "map")
        f.write(struct.pack("<I", len(ids)))
        ids.tofile(f)

    print(f"[map] {out} ({len(ids)})")


def dump_map(truth, hits_sorted, part_ids, part_ind, out):
    pid_to_idx = {int(pid): int(i) for i, pid in zip(part_ind, part_ids)}
    merged = hits_sorted[["hit_id"]].merge(truth[["hit_id", "particle_id"]], on="hit_id", how="left")
    ids = np.full(len(merged), 0xFFFFFFFF, np.uint32)
    for i, pid in enumerate(merged["particle_id"].fillna(0).astype(np.int64)):
        if pid in pid_to_idx:
            ids[i] = np.uint32(pid_to_idx[pid])

    with open(out, "wb") as f:
        write_header(f, "map")
        f.write(struct.pack("<I", len(ids)))
        ids.tofile(f)

    print(f"[map] {out} ({len(ids)})")


# ============================================================
# Main driver
# ============================================================
def main():
    p = argparse.ArgumentParser(description="Create SoA-compatible binaries with headers.")
    p.add_argument("--indir", default = "data/trackml/")
    p.add_argument("--outdir", default = "data/trackml/")
    p.add_argument("--detector", default= "data/trackml/detectors.csv")
    p.add_argument("--pattern", default= "event*")
    a = p.parse_args()

    indir, outdir = Path(a.indir), Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    geom = pd.read_csv(Path(a.detector))
    geom, mmap = build_module_index(geom)

    for truth_path in sorted(indir.glob(f"{a.pattern}-truth.csv")):
        prefix = truth_path.stem.replace("-truth", "")
        part_path = indir / f"{prefix}-particles.csv"
        hit_path = indir / f"{prefix}-hits.csv"
        if not part_path.exists() or not hit_path.exists():
            print(f"Skipping {prefix}: missing inputs")
            continue

        truth = pd.read_csv(truth_path)
        parts = pd.read_csv(part_path)
        hits = pd.read_csv(hit_path)
        n_mod = len(geom)
        hits = assign_modules(hits, mmap, n_mod)

        part_ind, part_ids = dump_particles(parts, outdir / f"{prefix}-particles.bin")
        hits_sorted = dump_hits(hits, n_mod, outdir / f"{prefix}-hits.bin")
        dump_map(truth, hits_sorted, part_ids, part_ind, outdir / f"{prefix}-map.bin")


if __name__ == "__main__":
    main()
