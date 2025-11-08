#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def plot_eff_fake_dupl(rootfile, outdir="plots"):
    os.makedirs(outdir, exist_ok=True)
    with uproot.open(rootfile) as f:
        keys = f.keys()
        print("[INFO] File keys:", keys)

        # Typical names in PixelTrackValidator histograms
        # (adjust if your file has different ones)
        h_eff = f["efficPt"] if "efficPt" in keys else None
        h_fake = f["fakeratePt"] if "fakeratePt" in keys else None
        h_dupl = f["duplicatesRatePt"] if "duplicatesRatePt" in keys else None

        def plot_hist(hist, title, xlabel, ylabel, filename):
            if hist is None:
                print(f"[WARN] Histogram {title} not found.")
                return
            h = hist.to_hist()
            x = h.axes[0].centers
            y = h.values()
            yerr = np.sqrt(h.variances()) if h.variances() is not None else None

            plt.figure(figsize=(7,5))
            plt.errorbar(x, y, yerr=yerr, fmt='o', color='blue', label=title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.grid(True, ls='--', alpha=0.6)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, filename))
            plt.close()
            print(f"[OK] Saved {filename}")

        # Plot Efficiency
        plot_hist(h_eff, "Tracking Efficiency vs pT", "pT [GeV]", "Efficiency", "efficiency_pt.png")

        # Plot Fake Rate
        plot_hist(h_fake, "Fake Rate vs pT", "pT [GeV]", "Fake Rate", "fakerate_pt.png")

        # Plot Duplicate Rate
        plot_hist(h_dupl, "Duplicate Rate vs pT", "pT [GeV]", "Duplicate Rate", "duplicates_pt.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot efficiency, fake rate and duplicate rate from PixelTrackValidatorFromHits.root")
    parser.add_argument("input", help="Path to PixelTrackValidatorFromHits.root")
    parser.add_argument("--outdir", default="plots", help="Output directory for plots")
    args = parser.parse_args()
    plot_eff_fake_dupl(args.input, args.outdir)
