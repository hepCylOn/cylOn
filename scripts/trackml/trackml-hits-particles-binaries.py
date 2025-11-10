#!/usr/bin/env python3
"""
Build SoA-compatible merged binaries (particles / hits / map) from TrackML CSVs.

Unit policy:
- CSV inputs (TrackML): spatial coordinates (vx, vy, vz, x, y, z) are in **mm**.
- This script converts those spatial coordinates to **cm** on write.
- Momenta (px, py, pz) are left unchanged (GeV/c).
- Geometry pitches (pitch_u, pitch_v) are read in **mm** and converted to **cm** for per-hit errors.
- Existing binary format is assumed cm-consistent; we just ensure conversions happen once here.

Outputs:
- Particles -> PAR1 header + uint32 nEvents + [ per-event: uint32 n + columns (vx,vy,vz,... in cm) ]
- Hits      -> TRH1 header + uint32 nEvents + [ per-event: (nHits,u32) (nModules,u32) moduleStart[...] + 13 columns ]
- Map       -> MAP1 header + uint32 nEvents + [ per-event: uint32 n + ids ]

Usage example:
  python3 trackml_hits_particles_binaries_cm.py \
    --indir data/trackml \
    --outdir data/trackml \
    --modules data/trackml/detectors.csv \
    --pattern "event000001*" \
    --volumes pixel strip_short
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import struct
import sys
import math
from typing import Dict, Tuple

from tqdm import tqdm
FORMAT_VERSION = 1
ENDIANNESS = 0x01020304
MAGIC = {
    "particles": b"PAR1",
    "hits": b"TRH1",
    "map": b"MAP1",
}

MM_TO_CM = np.float32(0.1)

def write_header(f, tag: str):
    f.write(MAGIC[tag])
    f.write(struct.pack("<I", FORMAT_VERSION))
    f.write(struct.pack("<I", ENDIANNESS))
    f.write(struct.pack("<I", 0))  # nEvents placeholder

def patch_nevents(f, n_events: int):
    f.seek(12)
    f.write(struct.pack("<I", np.uint32(n_events)))
    f.seek(0, 2)

# ============================================================
# Volume selection
# ============================================================
VOLUME_SETS = {
    "pixel":       {7, 8, 9},
    "strip_short": {12, 13, 14},
    "strip_long":  {16, 17, 18},
}

def resolve_volumes(volume_names):
    if not volume_names:
        return None
    vol_ids = set()
    for name in volume_names:
        if name not in VOLUME_SETS:
            raise ValueError(
                f"Unknown volume name '{name}'. Allowed: {list(VOLUME_SETS.keys())}"
            )
        vol_ids |= VOLUME_SETS[name]
    return vol_ids

# ============================================================
# Geometry helpers
# ============================================================
def build_module_index(geom_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[Tuple[int,int,int], int]]:
    """
    Build mapping: (volume_id, layer_id, module_id) -> module_index (0..nModules-1)
    and return a geometry table carrying per-module pitches.
    Required columns: volume_id, layer_id, module_id, pitch_u, pitch_v (pitches in mm).
    """
    needed = {"volume_id", "layer_id", "module_id", "pitch_u", "pitch_v"}
    if not needed.issubset(geom_df.columns):
        missing = needed - set(geom_df.columns)
        raise RuntimeError(f"detectors.csv missing columns: {sorted(missing)}")

    g = geom_df.sort_values(["volume_id", "layer_id", "module_id"]).reset_index(drop=True)
    g["module_index"] = np.arange(len(g), dtype=np.uint32)

    module_map = {(int(r.volume_id), int(r.layer_id), int(r.module_id)): int(r.module_index)
                  for _, r in g.iterrows()}
    return g, module_map

def per_module_local_errors_cm(geom_df_sorted: pd.DataFrame):
    """
    Compute per-module (sigmaU, sigmaV) in cm from pitch_u/pitch_v (in mm):
      sigma = (pitch / 2) / sqrt(12)
    Then convert mm -> cm (×0.1).
    Returns two float32 arrays of length nModules: sigmaU_cm, sigmaV_cm.
    """
    pitch_u_mm = geom_df_sorted["pitch_u"].to_numpy(dtype=np.float32)
    pitch_v_mm = geom_df_sorted["pitch_v"].to_numpy(dtype=np.float32)

    # sigma = half pitch / sqrt(12)
    denom = np.float32(math.sqrt(12.0))
    sigmaU_mm = (pitch_u_mm * 0.5) / denom
    sigmaV_mm = (pitch_v_mm * 0.5) / denom

    sigmaU_cm = (sigmaU_mm * MM_TO_CM).astype(np.float32)
    sigmaV_cm = (sigmaV_mm * MM_TO_CM).astype(np.float32)

    # Guard: any NaN/inf -> 0
    sigmaU_cm = np.nan_to_num(sigmaU_cm, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    sigmaV_cm = np.nan_to_num(sigmaV_cm, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    return sigmaU_cm, sigmaV_cm

def assign_modules_to_hits(hits_df: pd.DataFrame, module_map: dict, n_modules: int):
    needed = {"hit_id", "x", "y", "z", "volume_id", "layer_id", "module_id"}
    if not needed.issubset(hits_df.columns):
        missing = needed - set(hits_df.columns)
        raise RuntimeError(f"hits.csv missing columns: {sorted(missing)}")

    h = hits_df.copy()
    keys = list(
        zip(
            h["volume_id"].astype(int).to_numpy(),
            h["layer_id"].astype(int).to_numpy(),
            h["module_id"].astype(int).to_numpy(),
        )
    )
    idx = np.empty(len(h), dtype=np.int32)
    for i, k in enumerate(keys):
        try:
            idx[i] = module_map[k]
        except KeyError:
            raise RuntimeError(
                f"Hit references unknown module {k}. "
                "Ensure detectors.csv matches the selected volumes and the event."
            )
    if (idx < 0).any() or (idx >= n_modules).any():
        raise RuntimeError("Computed detectorIndex out of range [0, nModules).")
    h["module_index"] = idx.astype(np.uint32)
    return h

def build_module_start(hits_sorted_by_module: pd.DataFrame, n_modules: int):
    counts = np.bincount(
        hits_sorted_by_module["module_index"].to_numpy(), minlength=n_modules
    ).astype(np.uint32)
    module_start = np.zeros(n_modules + 1, dtype=np.uint32)
    module_start[1:] = np.cumsum(counts)
    return module_start

# ============================================================
# Physics helpers
# ============================================================
def compute_kinematics(df: pd.DataFrame):
    px = df["px"].to_numpy()
    py = df["py"].to_numpy()
    pz = df["pz"].to_numpy()

    pt = np.sqrt(px**2 + py**2)
    p = np.sqrt(px**2 + py**2 + pz**2)
    phi = np.arctan2(py, px)

    eta = np.zeros_like(pt, dtype=np.float32)
    mask = (p - pz) > 0
    safe = mask & np.isfinite(p) & (p > 0)
    eta[safe] = 0.5 * np.log((p[safe] + pz[safe]) / (p[safe] - pz[safe]))

    df["pt"] = pt.astype(np.float32)
    df["eta"] = eta.astype(np.float32)
    df["phi"] = phi.astype(np.float32)
    df["energy"] = p.astype(np.float32)
    df["mass"] = np.float32(0.0)

# ============================================================
# Streaming writers (per-event chunks into open files)
# ============================================================
def write_particles_chunk(f, parts_df: pd.DataFrame):
    df = parts_df.sort_values("particle_id").reset_index(drop=True)
    n = len(df)
    need = {"particle_id", "vx", "vy", "vz", "px", "py", "pz", "q"}
    if not need.issubset(df.columns):
        missing = need - set(df.columns)
        raise RuntimeError(f"particles.csv missing columns: {sorted(missing)}")

    # Convert spatial coordinates from mm -> cm in-memory
    for c in ("vx", "vy", "vz"):
        df[c] = df[c].astype(np.float32) * MM_TO_CM

    compute_kinematics(df)
    f.write(struct.pack("<I", np.uint32(n)))

    cols_f32 = [
        "vx", "vy", "vz",               # now in cm
        "px", "py", "pz", "energy",     # momenta/energy unchanged
        "pt", "eta", "phi", "mass",
    ]
    for c in cols_f32:
        df[c].to_numpy(dtype=np.float32).tofile(f)

    df["q"].to_numpy(dtype=np.int16).tofile(f)
    np.zeros(n, dtype=np.int32).tofile(f)           # pdgId (unknown -> 0)
    np.arange(n, dtype=np.uint32).tofile(f)         # partInd (0..n-1)

    return np.arange(n, dtype=np.uint32), df["particle_id"].to_numpy(dtype=np.int64)

def write_hits_chunk(
    f,
    hits_df: pd.DataFrame,
    n_modules: int,
    sigmaU_cm_per_module: np.ndarray,
    sigmaV_cm_per_module: np.ndarray,
):
    """
    Append one event's hits block:
      uint32 nHits, uint32 nModules, moduleStart[nModules+1], then 13 column arrays.
    xerrLocal / yerrLocal are filled from per-module pitch-derived sigmas (in cm).
    """
    h = hits_df.sort_values(["module_index", "hit_id"]).reset_index(drop=True)
    n_hits = len(h)

    # Convert global positions from mm -> cm
    x = (h["x"].to_numpy(dtype=np.float32) * MM_TO_CM).astype(np.float32)
    y = (h["y"].to_numpy(dtype=np.float32) * MM_TO_CM).astype(np.float32)
    z = (h["z"].to_numpy(dtype=np.float32) * MM_TO_CM).astype(np.float32)
    r = np.sqrt(x**2 + y**2).astype(np.float32)

    det = h["module_index"].to_numpy(dtype=np.uint16)
    module_start = build_module_start(h, n_modules)

    zeros_f = np.zeros(n_hits, dtype=np.float32)
    zeros_i16 = np.zeros(n_hits, dtype=np.int16)
    zeros_u32 = np.zeros(n_hits, dtype=np.uint32)

    # Map per-hit via detectorIndex
    xerrLocal = sigmaU_cm_per_module[det.astype(np.int32)].astype(np.float32)
    yerrLocal = sigmaV_cm_per_module[det.astype(np.int32)].astype(np.float32)

    # per-event header
    f.write(struct.pack("<II", np.uint32(n_hits), np.uint32(n_modules)))
    module_start.tofile(f)

    # xLocal, yLocal (unknown): keep zeros
    zeros_f.tofile(f)  # xLocal
    zeros_f.tofile(f)  # yLocal

    # xerrLocal, yerrLocal (now filled, in cm)
    xerrLocal.tofile(f)
    yerrLocal.tofile(f)

    # xGlobal, yGlobal, zGlobal, rGlobal (float32, in cm)
    for arr in (x, y, z, r):
        arr.tofile(f)

    # iphi (int16)
    zeros_i16.tofile(f)

    # chargeAndStatus (uint32)
    zeros_u32.tofile(f)

    # clusterSizeX, clusterSizeY (int16)
    zeros_i16.tofile(f)
    zeros_i16.tofile(f)

    # detectorIndex (uint16)
    det.tofile(f)

    return h  # sorted (for map merge order)

def write_map_chunk(f, truth_df: pd.DataFrame,
                    hits_sorted: pd.DataFrame,
                    particle_ids_sorted: np.ndarray,
                    partInd_sorted: np.ndarray):
    need = {"hit_id", "particle_id"}
    if not need.issubset(truth_df.columns):
        missing = need - set(truth_df.columns)
        raise RuntimeError(f"truth.csv missing columns: {sorted(missing)}")

    pid_to_index = {int(pid): int(idx) for idx, pid in zip(partInd_sorted, particle_ids_sorted)}

    merged = hits_sorted[["hit_id"]].merge(
        truth_df[["hit_id", "particle_id"]],
        on="hit_id", how="left"
    )

    ids = np.full(len(merged), 0xFFFFFFFF, dtype=np.uint32)
    col = merged["particle_id"].astype("Int64")
    for i in range(len(ids)):
        val = col.iat[i]
        if pd.notna(val):
            pid = int(val)
            if pid in pid_to_index:
                ids[i] = np.uint32(pid_to_index[pid])

    f.write(struct.pack("<I", np.uint32(len(ids))))
    ids.tofile(f)

# ============================================================
# Main
# ============================================================
def main():
    p = argparse.ArgumentParser(
        description="Create merged SoA-compatible binaries (particles, hits, map) from TrackML CSVs; converts CSV spatial coords mm->cm on write."
    )
    p.add_argument("--indir", default="data/trackml/",
                   help="Input directory with *-truth.csv, *-particles.csv, *-hits.csv")
    p.add_argument("--outdir", default="data/trackml/",
                   help="Output directory for merged *.bin files")
    p.add_argument("--modules", default="data/trackml/detectors.csv",
                   help="Detector geometry CSV (detectors.csv)")
    p.add_argument("--pattern", default="event*",
                   help="Event filename prefix pattern (default: event*)")
    p.add_argument("--volumes", nargs="*",
                   help="Subset of detector volumes to use: pixel, strip_short, strip_long. "
                        "Multiple allowed. If omitted, all volumes are used.")
    p.add_argument("--max-events", type=int, default=100)
    args = p.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    selected_vols = resolve_volumes(args.volumes)
    if selected_vols is None:
        print("[INFO] Using ALL volumes (no volume filter).")
    else:
        print(f"[INFO] Using volumes: {sorted(selected_vols)}")

    # --- Load geometry and apply volume filter (if any) ---
    geom_df = pd.read_csv(Path(args.modules))
    if selected_vols is not None:
        geom_df = geom_df[geom_df["volume_id"].isin(selected_vols)].copy()
    if geom_df.empty:
        raise RuntimeError("No modules left after applying volume filter!")

    geom_df_sorted, module_map = build_module_index(geom_df)
    n_modules = int(len(geom_df_sorted))
    print(f"[INFO] Geometry: {n_modules} modules after volume filter. (pitches in mm, errors converted to cm)")

    # Per-module local errors (cm) from pitch_u/pitch_v [mm]
    sigmaU_cm_per_module, sigmaV_cm_per_module = per_module_local_errors_cm(geom_df_sorted)

    # --- Prepare merged output files and write headers with placeholder nEvents ---
    parts_path = outdir / "particles.bin"
    hits_path  = outdir / "hits.bin"
    map_path   = outdir / "map.bin"

    parts_f = open(parts_path, "wb")
    hits_f  = open(hits_path, "wb")
    map_f   = open(map_path, "wb")
    try:
        write_header(parts_f, "particles")
        write_header(hits_f,  "hits")
        write_header(map_f,   "map")

        n_evt = 0
        sorted_files = sorted(indir.glob(f"{args.pattern}-truth.csv"))
        maxEv = min(len(sorted_files),args.max_events)
        for ip in tqdm(range(maxEv),desc="Events conversion:"):
            truth_path = sorted_files[ip]
            if n_evt > args.max_events:
                continue
            prefix = truth_path.stem.replace("-truth", "")
            part_path = indir / f"{prefix}-particles.csv"
            hit_path  = indir / f"{prefix}-hits.csv"

            if not part_path.exists() or not hit_path.exists():
                print(f"[WARN] Skipping {prefix}: missing particles or hits CSV")
                continue

            truth_df = pd.read_csv(truth_path)
            parts_df = pd.read_csv(part_path)
            hits_df  = pd.read_csv(hit_path)

            if selected_vols is not None:
                hits_df = hits_df[hits_df["volume_id"].isin(selected_vols)].copy()

            if hits_df.empty:
                print(f"[WARN] Event {prefix}: no hits left after volume filter, skipping.")
                continue

            # Particles (convert spatial coords mm->cm)
            partInd, part_ids = write_particles_chunk(parts_f, parts_df)

            # Assign module indices to hits via geometry mapping
            hits_with_idx = assign_modules_to_hits(hits_df, module_map, n_modules)

            # Hits + moduleStart (+ local errors from pitches, in cm). Converts x,y,z mm->cm.
            hits_sorted = write_hits_chunk(
                hits_f,
                hits_with_idx,
                n_modules,
                sigmaU_cm_per_module,
                sigmaV_cm_per_module,
            )

            # Map (truth association -> particle index)
            write_map_chunk(map_f, truth_df, hits_sorted, part_ids, partInd)

            n_evt += 1

        if n_evt == 0:
            print("[WARN] No events processed — check your --pattern / directory.")
        else:
            patch_nevents(parts_f, n_evt)
            patch_nevents(hits_f,  n_evt)
            patch_nevents(map_f,   n_evt)
            print(f"\n[OK] Finished. Events merged: {n_evt}")
            print(f"[OUT] {parts_path}")
            print(f"[OUT] {hits_path}")
            print(f"[OUT] {map_path}")
    finally:
        parts_f.close()
        hits_f.close()
        map_f.close()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        sys.exit(1)
