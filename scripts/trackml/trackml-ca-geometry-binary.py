#!/usr/bin/env python3
import argparse, json, struct
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------
# Volume selection
# ------------------------------
VOLUME_SETS = {
    "pixel":       {7, 8, 9},
    "strip_short": {12, 13, 14},
    "strip_long":  {16, 17, 18},
}

def resolve_volumes(names):
    if not names:
        return None
    out = set()
    for n in names:
        if n not in VOLUME_SETS:
            raise ValueError(f"Unknown volume name '{n}'. Allowed: {list(VOLUME_SETS.keys())}")
        out |= VOLUME_SETS[n]
    return out

# ------------------------------
# Helpers
# ------------------------------
def pick(df, names):
    for n in names:
        if n in df.columns:
            return n
    raise KeyError(f"None of {names} found in columns: {list(df.columns)}")

def eta_to_theta_max_rad(eta_max):
    # |eta| = -ln tan(theta/2)  => theta = 2*atan(exp(-|eta|))
    return abs(2.0 * np.arctan(np.exp(-abs(float(eta_max)))))

def ray_intersects_box(angle_rad, box_row):
    """
    Intersect a ray from origin in (Z,R) with an axis-aligned box [minZ,maxZ]x[minR,maxR].
    angle_rad is w.r.t. +Z axis: dz=cos(theta), dr=sin(theta)
    Returns (hit?, tmin>=0)
    """
    dz, dr = np.cos(angle_rad), np.sin(angle_rad)
    tmin, tmax = -np.inf, np.inf

    # Z slab
    if abs(dz) < 1e-12:
        if not (box_row["minZ"] <= 0 <= box_row["maxZ"]):
            return False, np.inf
    else:
        tz1, tz2 = box_row["minZ"] / dz, box_row["maxZ"] / dz
        tmin, tmax = max(tmin, min(tz1, tz2)), min(tmax, max(tz1, tz2))

    # R slab
    if abs(dr) < 1e-12:
        if not (box_row["minR"] <= 0 <= box_row["maxR"]):
            return False, np.inf
    else:
        tr1, tr2 = box_row["minR"] / dr, box_row["maxR"] / dr
        tmin, tmax = max(tmin, min(tr1, tr2)), min(tmax, max(tr1, tr2))

    if tmax < max(tmin, 0.0):
        return False, np.inf
    return True, max(tmin, 0.0)

def line_sweep_layer_pairs(
    boxes_df,
    theta_max_rad,
    angle_step_deg=1.0,
    allow_skip=False,
    min_cells=2
):
    """
    Perform ray sweep for theta in [-theta_max, +theta_max] around +Z axis.
    boxes_df must have: box_id, minZ,maxZ,minR,maxR, layer_key, layer_idx

    Returns:
      edges_layers: set of (layer_idx_from, layer_idx_to)
      sequences: list of (theta_deg, [box_id sequence])  # only if len(seq) >= min_cells
    """
    theta_m = np.rad2deg(theta_max_rad)
    theta_M = 180 - np.rad2deg(theta_max_rad)
    if theta_m < theta_M:
        angles = np.deg2rad(np.arange(theta_m, theta_M, angle_step_deg))
    else:
        angles = np.deg2rad(np.arange(theta_M, theta_m, angle_step_deg))

    sequences = []
    edges = set()

    rows = list(boxes_df.itertuples(index=False))

    for a in angles:
        hits = []
        for row in rows:
            ok, t = ray_intersects_box(a, {
                "minZ": row.minZ, "maxZ": row.maxZ,
                "minR": row.minR, "maxR": row.maxR
            })
            if ok:
                hits.append((t, row.box_id, row.layer_idx))
        hits.sort(key=lambda x: x[0])
        seq_boxes = [bid for _, bid, _ in hits]
        seq_layers = [lid for _, _, lid in hits]

        # Only build edges/sequences if we have at least min_cells intersections
        if len(seq_layers) >= min_cells:
            # edges between consecutive layers
            for u, v in zip(seq_layers[:-1], seq_layers[1:]):
                if u != v:
                    edges.add((u, v))

            if allow_skip and len(seq_layers) >= (min_cells + 1):
                # skip exactly one layer: i -> k using (i,*,k)
                for i in range(len(seq_layers) - 2):
                    u, w = seq_layers[i], seq_layers[i+2]
                    if u != w:
                        edges.add((u, w))

            sequences.append((np.rad2deg(a), seq_boxes))

    return edges, sequences

# ------------------------------
# Main
# ------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Build CAModules binary, layerStarts, layer ZR boxes and layer-pair graph from detectors.csv and hit files."
    )
    ap.add_argument("--detector", default="data/trackml/detectors.csv")
    ap.add_argument("--hits", nargs="+", default=["data/trackml/event000001000-hits.csv"],
                    help="One or more *-hits.csv files")
    ap.add_argument("--outbin", default="data/trackml/CAModulesTrackML.bin",
                    help="Output CAModules binary")
    ap.add_argument("--outjson", default="data/trackml/CAPairsTrackML.json",
                    help="Output combined JSON (layerStarts + pair graph + meta)")
    ap.add_argument("--graph-png", default="graph.png", help="Output graph plot (Z–R with boxes)")
    ap.add_argument("--zr-png", default="zr_layers.png", help="Output Z–R plot colored by layer")
    ap.add_argument("--volumes", nargs="*", help="pixel, strip_short, strip_long (multiple allowed)")
    ap.add_argument("--eta", type=float, default=4.0, help="Absolute eta for sweep (sweep from -eta..+eta)")
    ap.add_argument("--angle-step", type=float, default=1.0, help="Angle step in degrees for sweep")
    ap.add_argument("--skip-connections", action="store_true", help="Allow edges that skip one layer in a sweep")
    ap.add_argument("--min-cells", type=int, default=2,
                    help="Minimum number of cells to form a tuple (default: 2, so a triplet)")
    # ESProducer config output
    ap.add_argument("--es-config-out", default=None,
                    help="Write a config JSON containing CAGeometryHostESProducerGenericUpgrade block")
    ap.add_argument("--base-config", default=None,
                    help="Optional existing config to merge the ES block into (other keys preserved)")
    args = ap.parse_args()

    outbin = args.outbin
    outjson = args.outjson
    selected_vols = resolve_volumes(args.volumes)
    if selected_vols is None:
        print("[INFO] No volume filter. Using all volumes.")
        chosen_volume_names = []
    else:
        print(f"[INFO] Volume filter: {sorted(selected_vols)}")
        chosen_volume_names = list(args.volumes)
        # suffix filenames with volume selection
        vols_suffix = '_'.join([v.capitalize() for v in args.volumes])
        outbin = outbin.replace(".bin", f"_{vols_suffix}.bin")
        outjson = outjson.replace(".json", f"_{vols_suffix}.json")

    # ------------------------------
    # Load detector modules
    # ------------------------------
    det = pd.read_csv(args.detector)
    det = det.sort_values(["volume_id", "layer_id", "module_id"]).reset_index(drop=True)

    # Optional volume filter on geometry (to be consistent with hits)
    if selected_vols is not None:
        det = det[det["volume_id"].isin(selected_vols)].reset_index(drop=True)

    # Column detection: centers and rotation (cx,cy,cz common in TrackML)
    cx = pick(det, ["center_x", "x", "pos_x", "cx"])
    cy = pick(det, ["center_y", "y", "pos_y", "cy"])
    cz = pick(det, ["center_z", "z", "pos_z", "cz"])
    rot_cols = [
        pick(det, ["rot_xu","rot11","rxx"]),
        pick(det, ["rot_xv","rot12","rxy"]),
        pick(det, ["rot_xw","rot13","rxz"]),
        pick(det, ["rot_yu","rot21","ryx"]),
        pick(det, ["rot_yv","rot22","ryy"]),
        pick(det, ["rot_yw","rot23","ryz"]),
        pick(det, ["rot_zu","rot31","rzx"]),
        pick(det, ["rot_zv","rot32","rzy"]),
        pick(det, ["rot_zw","rot33","rzz"]),
    ]

    # ------------------------------
    # Write CAModules binary
    # ------------------------------
    n_modules = len(det)
    with open(outbin, "wb") as f:
        f.write(struct.pack("<i", int(n_modules)))
        for _, r in det.iterrows():
            vals = [r[cx], r[cy], r[cz]] + [r[c] for c in rot_cols]
            np.asarray(vals, dtype=np.float32).tofile(f)
    print(f"[OK] CAModules: {n_modules} modules → {outbin}")

    outdir = Path(outbin).parent
    bs_path = outdir / "BSZero.bin"

    # the BS POD has 11 floats
    with open(bs_path, "wb") as bf:
        np.zeros(11, dtype=np.float32).tofile(bf)
    print(f"[OK] BeamSpotPOD (all zeros) → {bs_path}")

    # ------------------------------
    # Build layerStarts (layers = unique (volume_id, layer_id))
    # ------------------------------
    layer_counts = det.groupby(["volume_id", "layer_id"], sort=False).size().to_numpy()
    layerStarts = np.zeros(len(layer_counts) + 1, dtype=int)
    layerStarts[1:] = np.cumsum(layer_counts)
    nLayers = int(len(layer_counts))
    print(f"[OK] layerStarts built (nLayers={nLayers})")

    # ------------------------------
    # Load hits (one or more files), apply same volume filter
    # ------------------------------
    hits_list = []
    hits_files = args.hits
    print(f"[INFO] Loading hits from {len(hits_files)} file(s)...")
    for hf in hits_files:
        h = pd.read_csv(hf)
        if selected_vols is not None:
            h = h[h["volume_id"].isin(selected_vols)].copy()
        hits_list.append(h)
    hits = pd.concat(hits_list, ignore_index=True)
    # coordinates are in mm
    hits["R"] = np.sqrt(hits["x"]**2 + hits["y"]**2)

    # ------------------------------
    # Build layer boxes from hits per (volume_id, layer_id) with padding ±5 mm
    # ------------------------------
    pad = 5.0  # mm
    layer_keys = sorted(hits.groupby(["volume_id", "layer_id"]).groups.keys())
    layer_to_idx = {key: i for i, key in enumerate(layer_keys)}

    box_rows = []
    for key in layer_keys:
        v, L = key
        sel = hits[(hits["volume_id"] == v) & (hits["layer_id"] == L)]
        if sel.empty:
            continue
        minZ, maxZ = sel["z"].min() - pad, sel["z"].max() + pad
        minR, maxR = sel["R"].min() - pad, sel["R"].max() + pad
        box_rows.append({
            "box_id": len(box_rows),
            "volume_id": v,
            "layer_id": L,
            "layer_key": key,
            "layer_idx": layer_to_idx[key],
            "minZ": float(minZ), "maxZ": float(maxZ),
            "minR": float(minR), "maxR": float(maxR),
            "centerZ": float(0.5*(minZ+maxZ)),
            "centerR": float(0.5*(minR+maxR)),
        })
    boxes_df = pd.DataFrame(box_rows)

    # ------------------------------
    # Sweep and build pairs
    # ------------------------------
    theta_max = eta_to_theta_max_rad(args.eta)  # radians
    edges_layers, sequences = line_sweep_layer_pairs(
        boxes_df,
        theta_max,
        angle_step_deg=args.angle_step,
        allow_skip=args.skip_connections,
        min_cells=args.min_cells
    )

    # Unique pairs and startingPair
    pairs = sorted(edges_layers)
    nPairs = int(len(pairs))
    incoming = {u: set() for u in range(len(layer_keys))}
    for (u, v) in pairs:
        incoming[v].add(u)
    startingPair = [i for i, (u, v) in enumerate(pairs) if len(incoming[u]) == 0]

    # ------------------------------
    # Write combined JSON: layers + pairs + meta + (optional) sweep
    # ------------------------------
    combined = {
        "nModules": int(n_modules),
        "nLayers": int(nLayers),
        "layerStarts": layerStarts.tolist(),
        "layers": [{"index": int(layer_to_idx[k]), "volume_id": int(k[0]), "layer_id": int(k[1])}
                   for k in layer_keys],
        "pairs": [{"from": int(u), "to": int(v)} for (u, v) in pairs],
        "startingPair": [int(i) for i in startingPair],
        "meta": {
            "eta": float(args.eta),
            "angleStepDeg": float(args.angle_step),
            "volumes": chosen_volume_names,
            "skipConnections": bool(args.skip_connections),
            "minSeqBoxes": int(args.min_cells),
        },
        "sweep": [{"theta_deg": float(a), "seq_box_ids": [int(b) for b in seq]}
                  for (a, seq) in sequences]
    }
    with open(outjson, "w") as jf:
        json.dump(combined, jf, indent=2)
    print(f"[OK] Combined layer+pairs JSON → {outjson} "
          f"(layers={combined['nLayers']}, pairs={len(pairs)}, starting={combined['startingPair']})")

    # ------------------------------
    # ESProducer config block (optional)
    # ------------------------------
    if args.es_config_out:
        # Flattened pairGraph [u0,v0,u1,v1,...]
        pairGraph = []
        for (u, v) in pairs:
            pairGraph.extend([int(u), int(v)])

        # Per-layer bounds (ordered by layer_idx)
        boxes_sorted = boxes_df.sort_values("layer_idx")
        layer_minZ = boxes_sorted["minZ"].to_numpy()
        layer_maxZ = boxes_sorted["maxZ"].to_numpy()
        layer_maxR = boxes_sorted["maxR"].to_numpy()

        # Derive per-pair windows (envelopes of connected boxes)
        minZ_arr = np.empty(nPairs, dtype=np.float32)
        maxZ_arr = np.empty(nPairs, dtype=np.float32)
        maxR_arr = np.empty(nPairs, dtype=np.float32)
        for i, (u, v) in enumerate(pairs):
            minZ_arr[i] = float(min(layer_minZ[u], layer_minZ[v]))
            maxZ_arr[i] = float(max(layer_maxZ[u], layer_maxZ[v]))
            maxR_arr[i] = float(max(layer_maxR[u], layer_maxR[v]))

        # Loose phi, theta, dca cuts
        phiCuts = np.full(nPairs, 1000, dtype=np.int16)
        thetaCuts = np.full(nLayers, 1.0, dtype=np.float)
        dcaCuts = np.full(nLayers, 1.0, dtype=np.float)
        assert(len(pairGraph) == 2 * nPairs)
        assert(len(startingPair) <= nPairs)
        assert(max(startingPair) < nPairs)
        assert(len(minZ_arr) == nPairs  and len(maxZ_arr) == nPairs  and len(maxR_arr) == nPairs)
        assert(len(layerStarts) == nLayers + 1)
        
        es_block = {}
        es_block["BeamSpotESProducer"] = {
            "data": str(bs_path)
        }

        es_block["CAHitNtupletUpgrade"] = {
            "maxNumberOfDoublets": 10000000,
            "maxNumberOfTuples": 500000,

            "avgHitsPerTrack": 7.0,
            "avgCellsPerHit": 6.0,
            "avgCellsPerCell": 0.151,
            "avgTracksPerCell": 0.04,

            "minHitsPerNtuplet": 3,
            "minHitsForSharingCut": 10,
            "ptmin": 0.9,
            "hardCurvCut": 0.0328407225,
            "cellZ0Cut": 7.5,
            "cellPtCut": 0.85,

            "dzdrFact": 15.2,
            "minYsizeB1": -1,
            "minYsizeB2": -1,
            "maxDYsize12": -1,
            "maxDYsize": -1,
            "maxDYPred": -1,

            "useRiemannFit": False,
            "fitNas4": False,
            "earlyFishbone": True,
            "lateFishbone": False,
            "doStats": False,
            "doSharedHitCut": True,
            "dupPassThrough": False,
            "useSimpleTripletCleaner": True,

            "maxChi2": 5.0,
            "minPtCut": 0.9,
            "maxZip": 0.4,
            "maxTip": 12.0
        }

        es_block["CAGeometryHostESProducerGenericUpgrade"] = {
                "data": str(outbin),
                "nLayers": int(nLayers),
                "nPairs": int(nPairs),
                "nModules": int(n_modules),

                "layerStarts": [int(x) for x in layerStarts.tolist()],
                "pairGraph": [int(x) for x in pairGraph],
                "startingPairs": [int(x) for x in startingPair],
                "phiCuts": [int(x) for x in phiCuts.tolist()],
                "minZ": [float(x) for x in minZ_arr.tolist()],
                "maxZ": [float(x) for x in maxZ_arr.tolist()],
                "maxR": [float(x) for x in maxR_arr.tolist()],
                "thetaCuts": [float(x) for x in thetaCuts.tolist()],
                "dcaCuts": [float(x) for x in dcaCuts.tolist()]
            }

        if args.base_config:
            try:
                with open(args.base_config) as f:
                    base = json.load(f)
            except Exception:
                base = {}
            base.update(es_block)
            with open(args.es_config_out, "w") as f:
                json.dump(base, f, indent=2)
        else:
            with open(args.es_config_out, "w") as f:
                json.dump(es_block, f, indent=2)
        print(f"[OK] ESProducer config → {args.es_config_out}")

    # ------------------------------
    # Plots
    # ------------------------------
    # 1) Z–R with hits color-coded by layer
    plt.figure(figsize=(8, 8))
    for key, grp in hits.groupby(["volume_id", "layer_id"]):
        plt.scatter(grp["z"], grp["R"], s=1, alpha=0.35, label=f"V{key[0]}-L{key[1]}")
    for _, b in boxes_df.iterrows():
        plt.plot([b.minZ, b.maxZ, b.maxZ, b.minZ, b.minZ],
                 [b.minR, b.minR, b.maxR, b.maxR, b.minR],
                 "r-", linewidth=1.0, alpha=0.7)
    plt.xlabel("Z [mm]")
    plt.ylabel("R = sqrt(X^2 + Y^2) [mm]")
    plt.title("Hits in Z–R colored by layer")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(args.zr_png, dpi=160)
    plt.close()
    print(f"[OK] Saved ZR plot → {args.zr_png}")

    # 2) Graph plot in Z–R: draw layer boxes + edges between box centers
    plt.figure(figsize=(10, 8))
    # boxes
    for _, b in boxes_df.iterrows():
        plt.plot([b.minZ, b.maxZ, b.maxZ, b.minZ, b.minZ],
                 [b.minR, b.minR, b.maxR, b.maxR, b.minR],
                 "k-", linewidth=0.8, alpha=0.4)
    # centers
    centersZ = boxes_df.sort_values("layer_idx")["centerZ"].to_numpy()
    centersR = boxes_df.sort_values("layer_idx")["centerR"].to_numpy()
    # edges
    for (u, v) in pairs:
        zu, ru = centersZ[u], centersR[u]
        zv, rv = centersZ[v], centersR[v]
        plt.plot([zu, zv], [ru, rv], "C1-", alpha=0.7, linewidth=1.5)
    # highlight starting pairs
    for idx in startingPair:
        u, v = pairs[idx]
        zu, ru = centersZ[u], centersR[u]
        zv, rv = centersZ[v], centersR[v]
        plt.plot([zu, zv], [ru, rv], "C3-", alpha=0.9, linewidth=2.5)
    # nodes
    for i, (z0, r0) in enumerate(zip(centersZ, centersR)):
        k = layer_keys[i]
        plt.scatter([z0], [r0], s=40, c="C0", zorder=3)
        plt.text(z0, r0, f"V{k[0]}-L{k[1]}", fontsize=7, ha="left", va="bottom")

    plt.xlabel("Z [mm]")
    plt.ylabel("R [mm]")
    plt.title(f"Layer pair graph in Z–R (η sweep ±{args.eta}, step={args.angle_step}°, "
              f"min-seq={args.min_cells}" +
              (", skip-1" if args.skip_connections else "") + ")")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(args.graph_png, dpi=160)
    plt.close()
    print(f"[OK] Saved graph plot with boxes → {args.graph_png}")

if __name__ == "__main__":
    main()
