#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import uproot
import matplotlib.pyplot as plt

KIND_MAP = {0: "Efficiency", 1: "Duplicates", 2: "FakeRate", 3: "Resolution"}
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
        branch_names = set(t.keys())

        def v(name):
            if name not in branch_names:
                return np.array([], dtype=float)
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
            "trueNHitsValues": v("trueNHitsValues"),
            "trueNHitsCounts": v("trueNHitsCounts"),
        }

def truth_counts_for_coord(arrs, coord):
    selected = subset_by_kind_coord(arrs, kind_val=0, coord_vals=[coord])
    return selected["binCenter"], selected["den"]


def scale_counts(y, yerr, mode, reco_total, label):
    if mode == "none":
        return y, yerr

    denominator = reco_total if mode == "reco" else np.sum(y)
    if denominator <= 0:
        raise ValueError(
            f"Cannot apply '{mode}' normalization to '{label}': "
            "the normalization count is zero."
        )

    return y / denominator, yerr / denominator

def subset_by_kind_coord(arrs, kind_val=None, coord_vals=None):
    sel = np.ones_like(arrs["kind"], dtype=bool)
    if kind_val is not None:
        sel &= (arrs["kind"] == kind_val)
    if coord_vals is not None:
        coord_vals = np.array(list(coord_vals))
        sel &= np.isin(arrs["coord"], coord_vals)
    return {k: v[sel] for k, v in arrs.items()}

def plot_xy(ax, x, y, yerr=None, title="", xlabel="", ylabel="", isLog=False):
    if yerr is not None:
        ax.errorbar(x, y, yerr=yerr, fmt="o-", lw=1.25, ms=4, capsize=3)
    else:
        ax.plot(x, y, "o-", lw=1.25, ms=4)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if isLog: ax.set_xscale('log')
    ax.grid(True, alpha=0.3)

def ratio_with_error(num, num_err, den, den_err):
    ratio = np.full_like(num, np.nan, dtype=float)
    ratio_err = np.full_like(num, np.nan, dtype=float)

    valid = den > 0
    ratio[valid] = num[valid] / den[valid]
    ratio_err[valid] = np.sqrt(
        (num_err[valid] / den[valid]) ** 2
        + (num[valid] * den_err[valid] / den[valid] ** 2) ** 2
    )
    return ratio, ratio_err

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
            # plot_xy(ax, x, y, yerr, f"{kind_name} vs {COORD_MAP[coord_val]}", xlabel, kind_name)
            if "pT" in xlabel: plot_xy(ax, x, y, yerr, f"{kind_name} vs {COORD_MAP[coord_val]}", xlabel, kind_name, isLog=True)
            else: plot_xy(ax, x, y, yerr, f"{kind_name} vs {COORD_MAP[coord_val]}", xlabel, kind_name)
        outpath = outdir / f"{kind_name.lower()}_pt_eta_phi.pdf"
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
    outpath = outdir / "resolution_all.pdf"
    fig.savefig(outpath, dpi=150)
    if show: plt.show()
    plt.close(fig)

def make_counts_plots(arrs, counts, outdir: Path, show: bool, normalization: str, show_ratio: bool):
    if counts is None:
        print("Counts tree not found — skipping counts plots.")
        return

    if normalization == "none":
        ylabel = "Counts"
    elif normalization == "reco":
        ylabel = "Entries / reconstructed entries"
    else:
        ylabel = "Fraction of entries"

    if show_ratio:
        fig, axes = plt.subplots(
            2, 3, figsize=(16, 7), sharex="col",
            gridspec_kw={"height_ratios": [3, 1]},
            constrained_layout=True,
        )
        plot_axes = axes[0, :]
        ratio_axes = axes[1, :]
    else:
        fig, plot_axes = plt.subplots(
            1, 3, figsize=(14, 4), constrained_layout=True,
        )
        ratio_axes = [None, None, None]
    triplets = [
        (0, "ptCenters", "ptCounts", "Track counts vs pT", "pT [GeV]"),
        (1, "etaCenters", "etaCounts", "Track counts vs η", "η"),
        (2, "phiCenters", "phiCounts", "Track counts vs φ", "φ [rad]"),
    ]

    for col, (coord, cx, cy, title, xlabel) in enumerate(triplets):
        ax = plot_axes[col]
        ratio_ax = ratio_axes[col]

        reco_x = counts[cx]
        reco_y = counts[cy]
        reco_total = np.sum(reco_y)
        reco_y_plot, reco_err = scale_counts(
            reco_y, poisson_error(reco_y), normalization, reco_total, title
        )

        truth_x, truth_y = truth_counts_for_coord(arrs, coord)
        truth_y_plot, truth_err = scale_counts(
            truth_y, poisson_error(truth_y), normalization, reco_total, title
        )

        if len(reco_x) != len(truth_x) or not np.allclose(reco_x, truth_x):
            raise ValueError(f"Truth and reconstructed bin centers differ for '{title}'")

        ax.errorbar(
            reco_x, reco_y_plot, yerr=reco_err, fmt="o-", lw=1.25, ms=4,
            capsize=3, label="Reconstructed tracks", color="steelblue"
        )
        ax.errorbar(
            truth_x, truth_y_plot, yerr=truth_err, fmt="s--", lw=1.25, ms=4,
            capsize=3, label="Truth particles", color="orange", mfc="none"
        )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend()

        if show_ratio:
            ratio, ratio_err = ratio_with_error(
                reco_y_plot, reco_err, truth_y_plot, truth_err
            )
            ratio_ax.errorbar(
                reco_x, ratio, yerr=ratio_err, fmt="o", ms=4, capsize=3
            )
            ratio_ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
            ratio_ax.set_xlabel(xlabel)
            ratio_ax.set_ylabel("Reco / Truth")
            ratio_ax.grid(True, alpha=0.3)

        if "pT" in xlabel:
            ax.set_xscale("log")
            # ax.set_yscale("log")
            if show_ratio:
                ratio_ax.set_xscale("log")

    fig.savefig(outdir / "counts_pt_eta_phi.pdf", dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    if show_ratio:
        fig, (ax, ratio_ax) = plt.subplots(
            2, 1, figsize=(6, 7), sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
            constrained_layout=True,
        )
    else:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
        ratio_ax = None

    reco_x = counts["nHitsValues"]
    reco_y = counts["nHitsCounts"]
    reco_total = np.sum(reco_y)
    reco_y_plot, reco_err = scale_counts(
        reco_y, poisson_error(reco_y), normalization, reco_total,
        "reconstructed nHits",
    )
    ax.errorbar(
        reco_x, reco_y_plot, yerr=reco_err, fmt="o", ms=4, capsize=3,
        label="Reconstructed tracks",
    )
    ax.vlines(reco_x, 0, reco_y_plot, lw=1.0, alpha=0.5)

    truth_x = counts["trueNHitsValues"]
    truth_y = counts["trueNHitsCounts"]
    if truth_x.size:
        truth_y_plot, truth_err = scale_counts(
            truth_y, poisson_error(truth_y), normalization, reco_total,
            "truth nHits",
        )
        ax.errorbar(
            truth_x, truth_y_plot, yerr=truth_err, fmt="s", ms=4, capsize=3,
            label="Truth particles",
        )
        ax.vlines(truth_x, 0, truth_y_plot, lw=1.0, alpha=0.5)

        if show_ratio:
            # Align the distributions in case they contain different nHits values.
            all_nhits = np.union1d(reco_x, truth_x)
            reco_aligned = np.zeros(len(all_nhits))
            reco_err_aligned = np.zeros(len(all_nhits))
            truth_aligned = np.zeros(len(all_nhits))
            truth_err_aligned = np.zeros(len(all_nhits))
    
            reco_indices = np.searchsorted(all_nhits, reco_x)
            truth_indices = np.searchsorted(all_nhits, truth_x)
            reco_aligned[reco_indices] = reco_y_plot
            reco_err_aligned[reco_indices] = reco_err
            truth_aligned[truth_indices] = truth_y_plot
            truth_err_aligned[truth_indices] = truth_err
    
            ratio, ratio_err = ratio_with_error(
                reco_aligned, reco_err_aligned,
                truth_aligned, truth_err_aligned,
            )
            ratio_ax.errorbar(
                all_nhits, ratio, yerr=ratio_err, fmt="o", ms=4, capsize=3
            )
            ratio_ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    else:
        print("Truth nHits branches not found; rerun validation with the updated producer.")

    ax.set_title("Track and truth-particle counts vs nHits")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend()

    if show_ratio:
        ratio_ax.set_xlabel("nHits")
        ratio_ax.set_ylabel("Reco / truth")
        # ratio_ax.set_yscale("log")
        ratio_ax.grid(True, alpha=0.3)

    fig.savefig(outdir / "counts_nhits.pdf", dpi=150)
    if show:
        plt.show()
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Plot PixelTrackValidation results (Validation + Counts trees)")
    ap.add_argument("-i", "--input", required=True, help="path to pixelTrackValidation_ntuple.root")
    ap.add_argument("-o", "--outdir", default="plots_validation", help="output directory for plots")
    ap.add_argument("--show", action="store_true", help="show interactive windows")
    ap.add_argument("-n", "--normalization", choices=("none", "reco", "unit"), default="none",
        help=(
            "count-plot normalization: none=raw counts, "
            "reco=divide both distributions by the reconstructed total, "
            "unit=normalize each distribution separately to 1"
        ),
    )
    ap.add_argument("-r", "--ratio-plots", action="store_true", help=("include reconstructed/truth ratio panels in the count plots"))

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
    make_counts_plots(val, cnt, outdir, args.show, args.normalization, show_ratio=args.ratio_plots)

    print(f"Saved plots to: {outdir.resolve()}")

if __name__ == "__main__":
    main()
