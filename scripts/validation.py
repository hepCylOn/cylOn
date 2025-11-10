#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import uproot
import matplotlib.pyplot as plt

KIND_MAP = {0: "Efficiency", 1: "FakeRate", 2: "Duplicates", 3: "Resolution"}
COORD_MAP = {0: "pT", 1: "eta", 2: "phi", 3: "d0", 4: "dz"}

def binomial_error(p, n):
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.sqrt(np.clip(p * (1.0 - p) / np.maximum(n, 1.0), 0.0, None))
    err[(n <= 0)] = 0.0
    return err

def poisson_error(n):
    return np.sqrt(np.maximum(n, 0.0))

def ensure_outdir(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

def load_validation_tree(filename: Path):
    with uproot.open(filename) as f:
        tree = f["Validation"]
        arrs = tree.arrays(["kind", "coord", "binCenter", "value", "num", "den"], library="np")
    return arrs

def load_counts_tree_or_none(filename: Path):
    with uproot.open(filename) as f:
        if "Counts" not in f:
            return None
        t = f["Counts"]
        # vector<float> branches -> jagged with one entry; take first
        def v(name):
            arr = t[name].array(library="np")
            if len(arr) == 0:
                return np.array([], dtype=float)
            return np.asarray(arr[0], dtype=float)
        return {
            "ptCenters":  v("ptCenters"),
            "ptCounts":   v("ptCounts"),
            "etaCenters": v("etaCenters"),
            "etaCounts":  v("etaCounts"),
            "phiCenters": v("phiCenters"),
            "phiCounts":  v("phiCounts"),
            "nHitsValues": v("nHitsValues"),
            "nHitsCounts": v("nHitsCounts"),
        }

def subset_by_kind_coord(arrs, kind_val=None, coord_vals=None):
    sel = np.ones_like(arrs["kind"], dtype=bool)
    if kind_val is not None:
        sel &= (arrs["kind"] == kind_val)
    if coord_vals is not None:
        coord_vals = np.array(list(coord_vals))
        sel &= np.isin(arrs["coord"], coord_vals)
    return {k: v[sel] for k, v in arrs.items()}

def plot_xy(ax, x, y, yerr=None, title="", xlabel="", ylabel=""):
    if yerr is not None:
        ax.errorbar(x, y, yerr=yerr, fmt="o-", lw=1.25, ms=4, capsize=3)
    else:
        ax.plot(x, y, "o-", lw=1.25, ms=4)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)

def make_eff_fake_dup_plots(arrs, outdir: Path, show: bool):
    for kind_val in (0, 1, 2):  # Eff, Fake, Dup
        kind_name = KIND_MAP[kind_val]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
        coords = [(0, "pT [GeV]"), (1, "η"), (2, "φ [rad]")]
        for ax, (coord_val, xlabel) in zip(axes, coords):
            s = subset_by_kind_coord(arrs, kind_val=kind_val, coord_vals=[coord_val])
            x = s["binCenter"]
            y = np.clip(s["value"], 0.0, 1.0)
            n = s["den"]
            yerr = binomial_error(y, n)
            plot_xy(ax, x, y, yerr, f"{kind_name} vs {COORD_MAP[coord_val]}", xlabel, kind_name)
        outpath = outdir / f"{kind_name.lower()}_pt_eta_phi.png"
        fig.savefig(outpath, dpi=150)
        if show: plt.show()
        plt.close(fig)

def make_resolution_plots(arrs, outdir: Path, show: bool):
    kind_val = 3
    coords = [
        (0, r"pT bin center [GeV]", r"σ(ΔpT) [GeV]"),
        (1, r"η bin center",        r"σ(Δη)"),
        (2, r"φ bin center [rad]",  r"σ(Δφ) [rad]"),
        (3, r"η bin center",        r"σ(Δd0) [cm]"),
        (4, r"η bin center",        r"σ(Δdz) [cm]"),
    ]
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), constrained_layout=True)
    for ax, (coord_val, xlabel, ylabel) in zip(axes, coords):
        s = subset_by_kind_coord(arrs, kind_val=kind_val, coord_vals=[coord_val])
        x = s["binCenter"]
        y = s["value"]
        plot_xy(ax, x, y, None, f"Resolution vs {COORD_MAP[coord_val]}", xlabel, ylabel)
    outpath = outdir / "resolution_all.png"
    fig.savefig(outpath, dpi=150)
    if show: plt.show()
    plt.close(fig)

def make_counts_plots(counts, outdir: Path, show: bool):
    if counts is None:
        print("Counts tree not found — skipping counts plots.")
        return

    # pT / eta / phi counts with Poisson errors
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    triplets = [
        ("ptCenters",  "ptCounts",  "Track counts vs pT",  "pT [GeV]"),
        ("etaCenters", "etaCounts", "Track counts vs η",   "η"),
        ("phiCenters", "phiCounts", "Track counts vs φ",   "φ [rad]"),
    ]
    for ax, (cx, cy, title, xlabel) in zip(axes, triplets):
        x = counts[cx]
        y = counts[cy]
        yerr = poisson_error(y)
        plot_xy(ax, x, y, yerr, title, xlabel, "Counts")
    outpath = outdir / "counts_pt_eta_phi.png"
    fig.savefig(outpath, dpi=150)
    if show: plt.show()
    plt.close(fig)

    # nHits (discrete)
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
    x = counts["nHitsValues"]
    y = counts["nHitsCounts"]
    yerr = poisson_error(y)
    # draw as markers + stems to emphasize discreteness
    markerline = ax.errorbar(x, y, yerr=yerr, fmt="o", ms=4, capsize=3)
    ax.vlines(x, 0, y, lw=1.0, alpha=0.5)
    ax.set_title("Track counts vs nHits")
    ax.set_xlabel("nHits")
    ax.set_ylabel("Counts")
    ax.grid(True, alpha=0.3)
    outpath = outdir / "counts_nhits.png"
    fig.savefig(outpath, dpi=150)
    if show: plt.show()
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Plot PixelTrackValidation results (Validation + Counts trees)")
    ap.add_argument("-i", "--input", required=True, help="path to pixelTrackValidation_ntuple.root")
    ap.add_argument("-o", "--outdir", default="plots_validation", help="output directory for plots")
    ap.add_argument("--show", action="store_true", help="show interactive windows")
    args = ap.parse_args()

    inpath = Path(args.input)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    # Load trees
    val = load_validation_tree(inpath)
    cnt = load_counts_tree_or_none(inpath)

    # Quick summary
    n_entries = len(val["kind"])
    print(f"Loaded {n_entries} rows from 'Validation'.")
    if cnt is None:
        print("No 'Counts' tree found.")
    else:
        print(f"Loaded 'Counts' vectors: "
              f"pt={len(cnt['ptCenters'])}, eta={len(cnt['etaCenters'])}, "
              f"phi={len(cnt['phiCenters'])}, nHits={len(cnt['nHitsValues'])}")

    # Plots
    make_eff_fake_dup_plots(val, outdir, args.show)
    make_resolution_plots(val, outdir, args.show)
    make_counts_plots(cnt, outdir, args.show)

    print(f"Saved plots to: {outdir.resolve()}")

if __name__ == "__main__":
    main()
