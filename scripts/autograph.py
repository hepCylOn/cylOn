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

def line_sweep(boxes_df, angle_step=1.0, min_boxes=3):
    angles = np.arange(-180, 180 + 1e-9, angle_step)
    sequences = []  # [(angle, [box ids])]

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
        for u, v in zip(seq[:-1], seq[1:]):
            if u != v:
                edges.add((u, v))

    edges_df = pd.DataFrame(sorted(edges), columns=["box_from", "box_to"])
    return sequences, edges_df


def ray_intersects_box(angle_deg, box):
    a = np.radians(angle_deg)
    dz, dr = np.cos(a), np.sin(a)
    tmin, tmax = -np.inf, np.inf
    for (minv, maxv, d) in [(box["minZ"], box["maxZ"], dz), (box["minR"], box["maxR"], dr)]:
        if abs(d) < 1e-9:
            if 0 < minv or 0 > maxv:
                return False, np.inf
        else:
            t1, t2 = (minv / d), (maxv / d)
            tmin, tmax = max(tmin, min(t1, t2)), min(tmax, max(t1, t2))
    if tmax < max(tmin, 0):
        return False, np.inf
    return True, max(tmin, 0)


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

    # --- Tolerances ---
    toleranceZ = args.tolerance_factor * np.min(np.diff(z_centers[peaks_z]))
    toleranceR = args.tolerance_factor * np.min(np.diff(r_centers[peaks_r]))
    print(f"[INFO] Derived toleranceZ = {toleranceZ:.3f}, toleranceR = {toleranceR:.3f}")

    # --- Step 2: Apply adaptive boxing ---
    z_step = toleranceZ
    r_step = toleranceR

    z_min, z_max = z.min(), z.max()
    z_cuts = [z_min - 5*z_step]
    current_z = z_min
    last_point_z = current_z
    while current_z <= z_max + 5*z_step:
        mask = (z >= current_z) & (z < current_z + z_step)
        if mask.any():
            last_point_z = current_z
        elif current_z - last_point_z >= toleranceZ:
            z_cuts.append(current_z)
            last_point_z = current_z
        current_z += z_step
    z_cuts.append(z_max + 5*z_step)
    z_cuts = sorted(set(z_cuts))

    # Build boxes in R inside each Z region
    boxes = []
    for i in range(len(z_cuts)-1):
        z_low, z_high = z_cuts[i], z_cuts[i+1]
        mask_z = (z >= z_low) & (z < z_high)
        if not mask_z.any():
            continue
        r_in_area = r[mask_z]
        r_min, r_max = r_in_area.min(), r_in_area.max()
        r_cuts = [r_min - 5*r_step]
        current_r = r_min
        last_point_r = current_r
        while current_r <= r_max + 5*r_step:
            mask_r = (r_in_area >= current_r) & (r_in_area < current_r + r_step)
            if mask_r.any():
                last_point_r = current_r
            elif current_r - last_point_r >= toleranceR:
                r_cuts.append(current_r)
                last_point_r = current_r
            current_r += r_step
        r_cuts.append(r_max + 5*r_step)
        r_cuts = sorted(set(r_cuts))
        iz = np.digitize(z[mask_z], z_cuts) - 1
        ir = np.digitize(r[mask_z], r_cuts) - 1
        for z_idx, r_idx, zz, rr in zip(iz, ir, z[mask_z], r[mask_z]):
            boxes.append((zz, rr, (z_idx, r_idx)))

    # --- Step 3: Plot segmentation ---
    zz = [b[0] for b in boxes]
    rr = [b[1] for b in boxes]
    labels = [b[2] for b in boxes]

    unique_labels = list(set(labels))
    color_map = {lbl: i for i, lbl in enumerate(unique_labels)}
    colors = [color_map[lbl] for lbl in labels]
    
    # Map tuple labels to integer cluster IDs for proper masking
    label_to_int = {lbl: i for i, lbl in enumerate(unique_labels)}
    int_labels = np.array([label_to_int[lbl] for lbl in labels])

    zz_arr = np.array(zz)
    rr_arr = np.array(rr)

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
        boxes_df["distance"] = np.sqrt(boxes_df["centerZ"]**2 + boxes_df["centerR"]**2)
        boxes_df["box_id"] = boxes_df.index
        boxes_df = boxes_df.sort_values("distance").reset_index(drop=True)
        sequences, edges_df = line_sweep(boxes_df, angle_step=1.0, min_boxes=3)
        sequences_df = pd.DataFrame(sequences,columns=["id","seq"])
        edges_df.to_csv("line_sweep_edges.csv", index=False)
        with open("line_sweep_sequences.txt", "w") as f:
            for a, seq in sequences:
                f.write(f"Angle {a:7.2f}°: " + " → ".join(f"box {b}" for b in seq) + "\n")

    if args.plot:
        
        # --- Plot clusters and bounding boxes ---
        plt.figure(figsize=(8, 8))
        for n,i in enumerate(boxes_df.box_id):
            zzz = zz_arr[int_labels==n]
            rrr = rr_arr[int_labels==n]
            plt.scatter(zzz, rrr, s=1, alpha=0.6, label=i)

        plt.gca().set_prop_cycle(None)

        # Draw red rectangles for each cluster bounding box
        for i_lbl, vertices in cluster_boxes:
            (z_min_adj, r_min_adj) = vertices[0]
            (z_max_adj, r_max_adj) = vertices[3]
            rect = patches.Rectangle(
                (z_min_adj, r_min_adj),
                z_max_adj - z_min_adj,
                r_max_adj - r_min_adj,
                linewidth=1.5,
                edgecolor='red',
                facecolor='none',
                alpha=0.8
            )
            plt.gca().add_patch(rect)

        plt.xlabel("Z")
        plt.ylabel("R = sqrt(X²+Y²)")
        plt.title(f"Bounding boxes per cluster (expanded by tol/10)\n{len(cluster_boxes)} boxes total")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(markerscale=3.0)
        plt.xlim()
        plt.show()
        plt.savefig("boxes.png")

        plt.gca().set_prop_cycle(None)
        if args.graph:
            #multilines = []
            for _, seq in sequences:
                zs = [boxes_df["centerZ"].values[s] for s in seq]
                rs = [boxes_df["centerR"].values[s]  for s in seq]
                #multilines.append((rs,zs))

                plt.plot(zs,rs,"-o",markersize=10,alpha=0.2)
            plt.savefig("boxesconnected.png")
        # # --- Plot: boxes with sample rays ---
        # plt.figure(figsize=(8,8))
        # plt.scatter(points_df["z"], points_df["r"], s=3, alpha=0.35, label="Points")
        # for _, b in boxes_df.iterrows():
        #     plt.plot([b["minZ"], b["maxZ"], b["maxZ"], b["minZ"], b["minZ"]],
        #             [b["minR"], b["minR"], b["maxR"], b["maxR"], b["minR"]],
        #             'r-', alpha=0.85, linewidth=1.0)

        # sample_every = max(1, len(sequences)//20)
        # for a, _ in sequences[::sample_every]:
        #     dz, dr = np.cos(np.radians(a)), np.sin(np.radians(a))
        #     L = max(z.max()-z.min(), r.max()-r.min()) * 2.0
        #     Zs = np.array([0, L*dz])
        #     Rs = np.array([0, L*dr])
        #     plt.plot(Zs, Rs, 'k--', alpha=0.25)

        # plt.xlabel("Z")
        # plt.ylabel("R = sqrt(X^2 + Y^2)")
        # plt.title(f"Line sweep intersections (angles with ≥ {min_boxes} boxes: {len(sequences)})")
        # plt.grid(True)
        # plt.legend()
        # plt.show()

    # # --- Assign points to boxes ---
    # def assign_box(zi, ri, boxes):
    #     for _, b in boxes.iterrows():
    #         if b["minZ"] <= zi <= b["maxZ"] and b["minR"] <= ri <= b["maxR"]:
    #             return int(b["box_id"])
    #     return -1

    # df_points = pd.DataFrame({"x": x, "y": y, "z": z, "r": r})
    # df_points["box"] = [assign_box(zi, ri, boxes_df) for zi, ri in zip(z, r)]
    # df_points = df_points.sort_values("box").reset_index(drop=True)

    # # --- Connections ---
    # connections = []
    # for a in np.arange(-180, 180, 1.0):
    #     hits = []
    #     for _, box in boxes_df.iterrows():
    #         ok, t = ray_intersects_box(a, box)
    #         if ok:
    #             hits.append((t, int(box["box_id"])))
    #     hits.sort(key=lambda x: x[0])
    #     seq = [h[1] for h in hits]
    #     if len(seq) > 1:
    #         for u, v in zip(seq[:-1], seq[1:]):
    #             connections.append((u, v))

    # connections_df = pd.DataFrame(sorted(set(connections)), columns=["box_from", "box_to"])
    # connections_json = connections_df.to_dict(orient="records")

    # # --- Save results ---
    # boxes_df.to_csv(f"{args.out_prefix}_boxes.csv", index=False)
    # df_points.to_csv(f"{args.out_prefix}_points_labeled.csv", index=False)
    # connections_df.to_csv(f"{args.out_prefix}_connections.csv", index=False)
    # with open(f"{args.out_prefix}_connections.json", "w") as f:
    #     json.dump(connections_json, f, indent=2)
    # print(f"[INFO] Outputs written to {args.out_prefix}_*.csv/json")

    # # --- Plots ---
    # if args.plot:
    #     # Boxes
    #     plt.figure(figsize=(8, 8))
    #     plt.scatter(z, r, c=df_points["box"], s=5, cmap="tab10", alpha=0.7)
    #     for _, b in boxes_df.iterrows():
    #         rect = patches.Rectangle(
    #             (b["minZ"], b["minR"]),
    #             b["maxZ"] - b["minZ"],
    #             b["maxR"] - b["minR"],
    #             linewidth=1.5,
    #             edgecolor="red",
    #             facecolor="none",
    #             alpha=0.8,
    #         )
    #         plt.gca().add_patch(rect)
    #     plt.xlabel("Z")
    #     plt.ylabel("R = sqrt(X²+Y²)")
    #     plt.title(f"Adaptive Z–R Boxes (expanded by tol/5) — {len(boxes_df)} boxes")
    #     plt.grid(True)
    #     plt.show()

    #     # Connections
    #     plt.figure(figsize=(8, 8))
    #     plt.scatter(z, r, c=df_points["box"], s=5, cmap="tab10", alpha=0.5)
    #     for _, row in connections_df.iterrows():
    #         b1 = boxes_df.loc[boxes_df["box_id"] == row["box_from"]].iloc[0]
    #         b2 = boxes_df.loc[boxes_df["box_id"] == row["box_to"]].iloc[0]
    #         plt.plot([b1["centerZ"], b2["centerZ"]], [b1["centerR"], b2["centerR"]], "k-", alpha=0.4)
    #     plt.xlabel("Z")
    #     plt.ylabel("R = sqrt(X²+Y²)")
    #     plt.title("Box Connections (-180° → 180° sweep)")
    #     plt.grid(True)
    #     plt.show()


if __name__ == "__main__":
    main()
