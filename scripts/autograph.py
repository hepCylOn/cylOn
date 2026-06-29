#!/usr/bin/env python3
"""
Adaptive histogram-based boxing and connectivity finder in Z–R space.

Steps:
1. Load input file (X, Y, Z columns)
2. Compute R = sqrt(X^2 + Y^2)
3. Detect peaks in Z and R histograms
4. Derive tolerances = (min peak spacing) / 10
5. Build boxes between consecutive peak intervals (expanded by ± tol/5)
6. Assign points to boxes
7. Compute box connections via ray sweep (-180° → 180°)
8. Output CSV/JSON results and plots
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import json
from collections import Counter

def ray_intersects_box(angle_deg, box_row):
    a = np.radians(angle_deg)
    dz, dr = np.cos(a), np.sin(a)
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

def line_sweep(boxes_df,
               angle_step=1.0,
               min_boxes=3,
               max_skip=1):

    angles = np.arange(-180, 180 + 1e-9, angle_step)
    sequences = []

    for a in angles:
        hits = []

        for _, b in boxes_df.iterrows():
            ok, t = ray_intersects_box(a, b)
            if ok:
                hits.append((t, int(b["box_id"])))

        hits.sort(key=lambda x: x[0])

        seq = [bid for _, bid in hits]

        if len(seq) >= min_boxes:
            sequences.append((a, seq))

    edges = set()

    for _, seq in sequences:

        n = len(seq)

        for i in range(n):

            for jump in range(1, max_skip + 2):

                j = i + jump

                if j >= n:
                    break

                u = seq[i]
                v = seq[j]

                if u != v:
                    edges.add((u, v))

    edges_df = pd.DataFrame(
        sorted(edges),
        columns=["box_from", "box_to"]
    )

    return sequences, edges_df

def hemisphere(row, eps):
    if row["maxZ"] < -eps:
        return 1

    return 0

def estimate_tolerances(z_region, r_region,
                        tolerance_factor=0.1,
                        bins_z=200,
                        bins_r=200):

    counts_z, edges_z = np.histogram(
        z_region,
        bins=bins_z,
        density=True
    )

    smooth_z = gaussian_filter1d(counts_z, sigma=2)

    peaks_z, _ = find_peaks(
        smooth_z,
        prominence=0.0005,
        distance=5
    )

    z_centers = 0.5*(edges_z[:-1] + edges_z[1:])

    if len(peaks_z) > 0:
        max_val = np.max(smooth_z[peaks_z])
        peaks_z = np.array([
            p for p in peaks_z
            if smooth_z[p] >= 0.5*max_val
        ])

    if len(peaks_z) < 2:
        toleranceZ = np.std(z_region)/10.
    else:
        toleranceZ = (
            tolerance_factor *
            np.min(np.diff(z_centers[peaks_z]))
        )

    counts_r, edges_r = np.histogram(
        r_region,
        bins=bins_r,
        density=True
    )

    smooth_r = gaussian_filter1d(counts_r, sigma=2)

    peaks_r, _ = find_peaks(
        smooth_r,
        prominence=0.02,
        distance=5
    )

    r_centers = 0.5*(edges_r[:-1] + edges_r[1:])

    if len(peaks_r) > 0:
        max_val = np.max(smooth_r[peaks_r])
        peaks_r = np.array([
            p for p in peaks_r
            if smooth_r[p] >= 0.5*max_val
        ])

    if len(peaks_r) < 2:
        toleranceR = np.std(r_region)/10.
    else:
        toleranceR = (
            tolerance_factor *
            np.min(np.diff(r_centers[peaks_r]))
        )

    return toleranceZ, toleranceR

def get_tolerances(r_value, regional_tolerances):

    for reg in regional_tolerances:

        if reg["rmin"] <= r_value < reg["rmax"]:
            return reg["tolZ"], reg["tolR"]

    return (
        regional_tolerances[-1]["tolZ"],
        regional_tolerances[-1]["tolR"]
    )

def main():
    ap = argparse.ArgumentParser(description="Adaptive boxing and connectivity in Z–R plane.")
    ap.add_argument("--input", required=True, help="Input CSV or TXT file with X,Y,Z columns.")
    ap.add_argument("--x-col", type=int, default=4, help="0-based column index for X.")
    ap.add_argument("--y-col", type=int, default=5, help="0-based column index for Y.")
    ap.add_argument("--z-col", type=int, default=6, help="0-based column index for Z.")
    ap.add_argument("--delimiter", default=",", help="Delimiter (default ','). Use 'auto' for autodetect.")
    ap.add_argument("--out-prefix", default="results/adaptive", help="Output prefix for CSV and plots.")
    ap.add_argument("--plot", action="store_true", help="Show plots.")
    ap.add_argument("--graph", action="store_true", help="Do graph.")
    ap.add_argument("--tolerance-factor", type=float, default=0.1, help="Factor to derive tolerances from peak spacing.")
    ap.add_argument("--box-expansion-factor", type=float, default=0.1, help="Factor to expand box boundaries.")
    args = ap.parse_args()

    # --- Load input ---
    df = pd.read_csv(args.input, header=None, sep=None if args.delimiter == "auto" else args.delimiter, engine="python", on_bad_lines="skip")
    x, y, z = df[args.x_col], df[args.y_col], df[args.z_col]
    r = np.sqrt(x**2 + y**2)
    print(f"[INFO] Loaded {len(z)} points")

    # Limits in R that need to be input by the user. This follows the geometry of the detector, so it should be straightforward to find from data, but it is like this for the moment
    region_edges = [3,20,72,np.inf]

    regions = [(region_edges[i], region_edges[i+1]) for i in range(0,len(region_edges)-1)]

    regional_tolerances = []

    for rmin, rmax in regions:

        mask = (r >= rmin) & (r < rmax)

        if np.sum(mask) < 10:
            continue

        tolZ, tolR = estimate_tolerances(
            z[mask],
            r[mask],
            tolerance_factor=args.tolerance_factor
        )

        regional_tolerances.append({
            "rmin": rmin,
            "rmax": rmax,
            "tolZ": tolZ,
            "tolR": tolR
        })

        print(
            f"[INFO] R=[{rmin},{rmax}] "
            f"-> tolZ={tolZ:.3f}, "
            f"tolR={tolR:.3f}"
        )

    # --- Histogram and peaks ---
    bins_z = 200
    bins_r = 200
    counts_z, edges_z = np.histogram(z, bins=bins_z, density=True)
    counts_r, edges_r = np.histogram(r, bins=bins_r, density=True)

    smooth_z = gaussian_filter1d(counts_z, sigma=2)
    smooth_r = gaussian_filter1d(counts_r, sigma=2)

    peaks_z, _ = find_peaks(smooth_z, prominence=0.0005, distance=5)
    peaks_r, _ = find_peaks(smooth_r, prominence=0.02, distance=5)

    z_centers = 0.5 * (edges_z[:-1] + edges_z[1:])
    r_centers = 0.5 * (edges_r[:-1] + edges_r[1:])

    # Remove spurious peaks (below 50% of max)
    if len(peaks_z):
        max_z_val = np.max(smooth_z[peaks_z])
        peaks_z = np.array([p for p in peaks_z if smooth_z[p] >= 0.5 * max_z_val])
    if len(peaks_r):
        max_r_val = np.max(smooth_r[peaks_r])
        peaks_r = np.array([p for p in peaks_r if smooth_r[p] >= 0.5 * max_r_val])

    if len(peaks_z) < 2:
        peaks_z = np.array([0, len(smooth_z) - 1])
    if len(peaks_r) < 2:
        peaks_r = np.array([0, len(smooth_r) - 1])
    
    if args.plot:
        plt.figure(figsize=(12, 5))

        # Z histogram
        plt.subplot(1, 2, 1)
        plt.plot(z_centers, smooth_z, color="steelblue", label="Smoothed Z histogram")
        plt.scatter(z_centers[peaks_z], smooth_z[peaks_z], color="red", s=40, label="Peaks")
        for p in z_centers[peaks_z]:
            plt.axvline(p, color="k", linestyle="--", alpha=0.3)
        plt.xlabel("Z")
        plt.ylabel("Normalized density")
        plt.title("Detected Peaks in Z Histogram")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)

        # R histogram
        plt.subplot(1, 2, 2)
        plt.plot(r_centers, smooth_r, color="darkorange", label="Smoothed R histogram")
        plt.scatter(r_centers[peaks_r], smooth_r[peaks_r], color="red", s=40, label="Peaks")
        for p in r_centers[peaks_r]:
            plt.axvline(p, color="k", linestyle="--", alpha=0.3)
        plt.xlabel("R = sqrt(X²+Y²)")
        plt.ylabel("Normalized density")
        plt.title("Detected Peaks in R Histogram")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)

        plt.tight_layout()
        plt.show()
        plt.savefig("peaks.png")
        plt.close()

    boxes = []

    for reg in regional_tolerances:

        rmin_reg = reg["rmin"]
        rmax_reg = reg["rmax"]

        toleranceZ = reg["tolZ"]
        toleranceR = reg["tolR"]

        z_step = toleranceZ/10.
        r_step = toleranceR/10.

        mask_region = (
            (r >= rmin_reg) &
            (r < rmax_reg)
        )

        z_region = z[mask_region]
        r_region = r[mask_region]

        if len(z_region) == 0:
            continue

        z_min = z_region.min()
        z_max = z_region.max()

        z_cuts = [z_min - 5*toleranceZ]

        current_z = z_min
        last_point_z = current_z

        while current_z <= z_max + 5*toleranceZ:

            mask_z = (
                (z_region >= current_z)
                &
                (z_region < current_z + z_step)
            )

            if mask_z.any():
                last_point_z = current_z

            elif current_z - last_point_z >= 2*toleranceZ:
                z_cuts.append(current_z)
                last_point_z = current_z

            current_z += z_step

        z_cuts.append(z_max + 5*toleranceZ)
        z_cuts = sorted(set(z_cuts))

        for i in range(len(z_cuts)-1):

            z_low = z_cuts[i]
            z_high = z_cuts[i+1]

            mask_z = (
                (z_region >= z_low)
                &
                (z_region < z_high)
            )

            if not mask_z.any():
                continue

            r_in_area = r_region[mask_z]

            r_min = r_in_area.min()
            r_max = r_in_area.max()

            r_cuts = [r_min - 5*toleranceR]

            current_r = r_min
            last_point_r = current_r

            while current_r <= r_max + 5*toleranceR:

                mask_r = (
                    (r_in_area >= current_r)
                    &
                    (r_in_area < current_r + r_step)
                )

                if mask_r.any():
                    last_point_r = current_r

                elif current_r - last_point_r >= 2*toleranceR:
                    r_cuts.append(current_r)
                    last_point_r = current_r

                current_r += r_step

            r_cuts.append(r_max + 5*toleranceR)
            r_cuts = sorted(set(r_cuts))

            iz = np.digitize(
                z_region[mask_z],
                z_cuts
            ) - 1

            ir = np.digitize(
                r_region[mask_z],
                r_cuts
            ) - 1

            for z_idx, r_idx, zz, rr in zip(
                    iz,
                    ir,
                    z_region[mask_z],
                    r_region[mask_z]):

                boxes.append(
                    (zz, rr, (z_idx, r_idx, rmin_reg))
                )

    min_occupancy = 20

    zz = [b[0] for b in boxes]
    rr = [b[1] for b in boxes]
    labels = [b[2] for b in boxes]

    occupancy = Counter(labels)

    valid_labels = {
        lbl for lbl, count in occupancy.items()
        if count >= min_occupancy
    }

    filtered = [
        (zv, rv, lbl)
        for zv, rv, lbl in zip(zz, rr, labels)
        if lbl in valid_labels
    ]

    zz = [f[0] for f in filtered]
    rr = [f[1] for f in filtered]
    labels = [f[2] for f in filtered]

    unique_labels = list(set(labels))
    color_map = {lbl: i for i, lbl in enumerate(unique_labels)}
    colors = [color_map[lbl] for lbl in labels]
    
    # Map tuple labels to integer cluster IDs for proper masking
    label_to_int = {lbl: i for i, lbl in enumerate(unique_labels)}
    int_labels = np.array([label_to_int[lbl] for lbl in labels])

    zz_arr = np.array(zz)
    rr_arr = np.array(rr)

    for i_lbl in np.unique(int_labels):

        mask = int_labels == i_lbl

    cluster_boxes = []

    # Compute one bounding box per integer label
    for i_lbl in np.unique(int_labels):
        mask = int_labels == i_lbl
        if not np.any(mask):
            continue
        z_min, z_max = np.min(zz_arr[mask]), np.max(zz_arr[mask])
        r_min, r_max = np.min(rr_arr[mask]), np.max(rr_arr[mask])

        # Apply tolerance offsets
        z_min_adj = z_min - toleranceZ * args.box_expansion_factor
        z_max_adj = z_max + toleranceZ * args.box_expansion_factor
        r_min_adj = r_min - toleranceR * args.box_expansion_factor
        r_max_adj = r_max + toleranceR * args.box_expansion_factor

        # Store vertices
        vertices = [
            (z_min_adj, r_min_adj),
            (z_min_adj, r_max_adj),
            (z_max_adj, r_min_adj),
            (z_max_adj, r_max_adj)
        ]
        cluster_boxes.append((i_lbl, vertices))

    if args.graph or args.plot:

        # Convert to boxes_df-compatible format
        boxes_data = []
        for box_id, verts in cluster_boxes:
            zs = [v[0] for v in verts]
            rs = [v[1] for v in verts]
            z_min, z_max = np.min(zs), np.max(zs)
            r_min, r_max = np.min(rs), np.max(rs)
            centerZ = 0.5 * (z_min + z_max)
            centerR = 0.5 * (r_min + r_max)
            distance = np.sqrt(centerZ**2 + centerR**2)
            boxes_data.append({
                "box_id": box_id,
                "minZ": z_min,
                "maxZ": z_max,
                "minR": r_min,
                "maxR": r_max,
                "centerZ": centerZ,
                "centerR": centerR,
                "distance": distance
            })

    if args.graph:

        boxes_df = pd.DataFrame(boxes_data)

        def r_region(row, region_edges):

            r = row["centerR"]

            for i in range(len(region_edges)-1):
            
                if region_edges[i] <= r < region_edges[i+1]:
                    return i

            return len(region_edges)-1

        # region_edges = [3,20,72,np.inf]

        boxes_df["r_region"] = boxes_df.apply(r_region, axis=1, args=(region_edges,))

        boxes_df["distance"] = np.sqrt(
            boxes_df["centerZ"]**2 + boxes_df["centerR"]**2
        )

        z_deadband = toleranceZ

        boxes_df["z_priority"] = boxes_df.apply(
            hemisphere,
            axis=1,
            eps=z_deadband
        )

        boxes_df = boxes_df.sort_values(
            by=[
                "r_region",
                "z_priority",
                "distance"
            ]
        ).reset_index(drop=True)

        boxes_df["box_id"] = np.arange(len(boxes_df))

        # Checks is boxes are horizontal or vertical to know if it is barrel or not
        boxes_df["isBarrel"] = (
            (boxes_df["maxZ"] - boxes_df["minZ"]) >
            (boxes_df["maxR"] - boxes_df["minR"])
        ).astype(int)

        sequences, edges_df = line_sweep(boxes_df, angle_step=1.0, min_boxes=3, max_skip=1)
        sequences_df = pd.DataFrame(sequences,columns=["id","seq"])
        edges_df.to_csv("line_sweep_edges.csv", index=False)

        with open("line_sweep_edges.txt", "w") as f:
            f.write("numberOfLayers: " + str(edges_df["box_to"][len(edges_df) - 1] + 1))
            f.write("\n\n" + "nPairs: " + str(len(edges_df)))
            f.write("\n\n" + "numberOfModules: " + str(edges_df["box_to"][len(edges_df) - 1] + 1))
            f.write("\n\n" + "layerPairs: ")
            for a in range(len(edges_df)):
                f.write(str(edges_df["box_from"][a]) + "," + str(edges_df["box_to"][a]) + ",")
            f.write("\n\n" + "startingPairs: ")
            for a in range(len(edges_df)): 
                if a < 3: f.write("1,")
                else: f.write("0,")
            f.write("\n\n" + "isBarrel: ")
            for a in boxes_df["isBarrel"]: f.write(str(a) + ",")
            f.write("\n\n" + "ptCuts: ")
            for a in range(len(edges_df)): f.write("0.5,")
            f.write("\n\n" + "layerStart: ")
            for a in range(edges_df["box_to"][len(edges_df) - 1] + 2): f.write(str(a) + ",")
            f.write("\n\n" + "phicuts: ")
            for a in range(len(edges_df)): f.write("100,")
            f.write("\n\n" + "minz: ")
            for a in range(len(edges_df)): f.write("0.5,")
            f.write("\n\n" + "maxz: ")
            for a in range(len(edges_df)): f.write("0.5,")
            f.write("\n\n" + "maxr: ")
            for a in range(len(edges_df)): f.write("0.5,")
            f.write("\n\n" + "dcaCuts: ")
            for a in range(edges_df["box_to"][len(edges_df) - 1] + 1): f.write("0.5,")
            f.write("\n\n" + "thetaCuts: ")
            for a in range(edges_df["box_to"][len(edges_df) - 1] + 1): f.write("0.5,")

        with open("line_sweep_sequences.txt", "w") as f:
            for a, seq in sequences:
                f.write(f"Angle {a:7.2f}°: " + " → ".join(f"box {b}" for b in seq) + "\n")

    if args.plot:
        
        for _, row in boxes_df.iterrows():

            rect = patches.Rectangle(
                (row["minZ"], row["minR"]),
                row["maxZ"] - row["minZ"],
                row["maxR"] - row["minR"],
                linewidth=1.5,
                edgecolor="red",
                facecolor="none"
            )

            plt.gca().add_patch(rect)

            plt.text(
                row["centerZ"],
                row["centerR"],
                str(int(row["box_id"])),
                fontsize=5,
                ha="center",
                va="center",
                fontweight="bold",
                bbox=dict(
                    facecolor="white",
                    edgecolor="black",
                    alpha=0.8
                )
            )

        for _, edge in edges_df.iterrows():

            b1 = boxes_df.loc[
                boxes_df["box_id"] == edge["box_from"]
            ].iloc[0]

            b2 = boxes_df.loc[
                boxes_df["box_id"] == edge["box_to"]
            ].iloc[0]

            plt.plot(
                [b1["centerZ"], b2["centerZ"]],
                [b1["centerR"], b2["centerR"]],
                "-",
                linewidth=1.5,
                alpha=0.3
            )

        plt.xlabel("Z")
        plt.ylabel("R")
        plt.title("Boxes conectadas")

        plt.grid(True, alpha=0.3)

        plt.savefig(
            "boxesconnected.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        for _, row in boxes_df.iterrows():

            rect = patches.Rectangle(
                (row["minZ"], row["minR"]),
                row["maxZ"] - row["minZ"],
                row["maxR"] - row["minR"],
                linewidth=1.5,
                edgecolor="red",
                facecolor="none"
            )

            plt.gca().add_patch(rect)

            plt.text(
                row["centerZ"],
                row["centerR"],
                str(int(row["box_id"])),
                fontsize=5,
                ha="center",
                va="center",
                fontweight="bold",
                bbox=dict(
                    facecolor="white",
                    edgecolor="black",
                    alpha=0.8
                )
            )

        plt.xlabel("Z")
        plt.ylabel("R")
        plt.title("Boxes")

        plt.grid(True, alpha=0.3)

        plt.savefig(
            "boxes.png",
            dpi=300,
            bbox_inches="tight"
        )

if __name__ == "__main__":
    main()
