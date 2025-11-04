#!/usr/bin/env python3
import json
import struct
import numpy as np
from pathlib import Path
import math

# ====================== Constants ==========================
MAGIC = b"TRH1"
VERSION = 1
ENDIANNESS_MARKER = 0x01020304

# ====================== Helpers ============================

def phi2short(phi: np.ndarray) -> np.ndarray:
    """Convert float phi [radians] to short (same as constexpr short phi2short)."""
    p2i = (np.iinfo(np.int16).max + 1) / math.pi
    return np.round(phi * p2i).astype(np.int16)


def pack_charge_and_status(event):
    """Pack status + charge fields into uint32 bitfield (SiPixelHitStatusAndCharge)."""
    n = len(event["xGlobal"])
    zeros_u32 = lambda: np.zeros(n, dtype=np.uint32)

    def get_uint32(name):
        arr = event.get(name, zeros_u32())
        return arr.astype(np.uint32, copy=False)

    isBigX = get_uint32("isBigX") & 1
    isOneX = get_uint32("isOneX") & 1
    isBigY = get_uint32("isBigY") & 1
    isOneY = get_uint32("isOneY") & 1
    qBin   = get_uint32("qBin")   & 0x7
    charge = get_uint32("charge") & 0xFFFFFF

    event["chargeAndStatus"] = (
        (isBigX << 0)
        | (isOneX << 1)
        | (isBigY << 2)
        | (isOneY << 3)
        | (qBin  << 4)
        | (charge << 7)
    ).astype(np.uint32)


def ensure_and_fill_missing(event):
    """Ensure required fields exist; fill optional ones with computed/zero values."""
    required = [
        "xLocal", "yLocal", "xerrLocal", "yerrLocal",
        "xGlobal", "yGlobal", "zGlobal", "detectorIndex"
    ]
    for key in required:
        if key not in event:
            raise ValueError(f"Missing required field: {key}")

    n = len(event["xGlobal"])
    zeros_f = lambda: np.zeros(n, dtype=np.float32)
    zeros_i16 = lambda: np.zeros(n, dtype=np.int16)
    zeros_u16 = lambda: np.zeros(n, dtype=np.uint16)
    zeros_i32 = lambda: np.zeros(n, dtype=np.int32)

    # Derived quantities
    if "rGlobal" not in event:
        event["rGlobal"] = np.sqrt(event["xGlobal"] ** 2 + event["yGlobal"] ** 2).astype(np.float32)

    if "iphi" not in event:
        phi = np.arctan2(event["yGlobal"], event["xGlobal"]).astype(np.float32)
        event["iphi"] = phi2short(phi)

    # Optional numeric defaults
    for key, func in [
        ("clusterSizeX", zeros_i16),
        ("clusterSizeY", zeros_i16),
        ("chargeAndStatus", zeros_i32),
    ]:
        if key not in event:
            event[key] = func()


def remap_and_sort_hits(event):
    """Remap detector indices, sort hits, compute moduleStart cumulative vector."""
    dets = event["detectorIndex"].astype(np.int32)
    unique_dets, inverse = np.unique(dets, return_inverse=True)
    event["detectorIndex"] = inverse.astype(np.uint16)
    nModules = len(unique_dets)

    order = np.argsort(event["detectorIndex"])
    for key in event.keys():
        event[key] = event[key][order]

    counts = np.bincount(event["detectorIndex"], minlength=nModules)
    moduleStart = np.zeros(nModules + 1, dtype=np.uint32)
    moduleStart[1:] = np.cumsum(counts)
    return nModules, moduleStart


def load_events_from_txt(file_path, mapping, delimiter=""):
    """Split file into events by delimiter and map columns."""
    with open(file_path) as f:
        lines = f.read().strip().splitlines()

    events, cur = [], []
    for line in lines:
        if line.strip() == delimiter:
            if cur:
                events.append(cur)
                cur = []
        else:
            cur.append(line)
    if cur:
        events.append(cur)

    parsed = []
    for lines in events:
        data = np.array([list(map(float, l.split(','))) for l in lines], dtype=np.float32)
        event = {}
        for key, idx in mapping.items():
            if idx < data.shape[1]:
                event[key] = data[:, idx]

        ensure_and_fill_missing(event)
        pack_charge_and_status(event)
        parsed.append(event)

    return parsed


def save_events_to_bin(events, output_path):
    """
    Binary layout:
    [magic TRH1][uint32 version][uint32 endianness][uint32 nEvents]
    For each event:
      [uint32 nHits][uint32 nModules]
      [uint32 (nModules+1) moduleStart]
      [columns: xLocal..detectorIndex]
    """
    with open(output_path, "wb") as f:
        # Header
        f.write(MAGIC)
        f.write(struct.pack("III", VERSION, ENDIANNESS_MARKER, len(events)))

        for ev in events:
            nModules, moduleStart = remap_and_sort_hits(ev)
            nHits = len(next(iter(ev.values())))
            f.write(struct.pack("II", nHits, nModules))
            f.write(moduleStart.tobytes())

            columns = [
                ("xLocal", np.float32), ("yLocal", np.float32),
                ("xerrLocal", np.float32), ("yerrLocal", np.float32),
                ("xGlobal", np.float32), ("yGlobal", np.float32),
                ("zGlobal", np.float32), ("rGlobal", np.float32),
                ("iphi", np.int16), ("chargeAndStatus", np.uint32),
                ("clusterSizeX", np.int16), ("clusterSizeY", np.int16),
                ("detectorIndex", np.uint16)
            ]
            for name, dtype in columns:
                f.write(ev[name].astype(dtype).tobytes())

    print(f"Wrote {len(events)} events to {output_path}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert hit txt → binary SoA for TrackingRecHitHost")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--delimiter", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.mapping) as f:
        mapping = json.load(f)

    all_events = []
    for fn in args.inputs:
        evs = load_events_from_txt(fn, mapping, args.delimiter)
        hits = np.sum([len(ev["xGlobal"]) for ev in evs])
        print(f"Loaded {len(evs)} events from {fn} with {hits} hits in total.")
        all_events.extend(evs)

    save_events_to_bin(all_events, args.output)
