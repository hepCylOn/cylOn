import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

parser.add_argument('-i', '--inpFileName', default='SimDoublets')
parser.add_argument('-o', '--outDirName', default='plots_simPixelTracks')
parser.add_argument('-p', '--percentageForCuts', type=float, default=99.0)
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-s', '--summary-only', action='store_true')
parser.add_argument('-m', '--multipage-pdf', action='store_true')

args = parser.parse_args()

inputFile = 'output' + args.inpFileName + '.txt'
outDir = args.outDirName

os.makedirs(outDir, exist_ok=True)

pdf = None
pdf_path = os.path.join(outDir, f"testSimPixelTrackHistos_{args.inpFileName}.pdf")

if args.multipage_pdf and not args.summary_only:
    pdf = PdfPages(pdf_path)

# Read non-empty lines
with open(inputFile, "r") as f:
    lines = [l.strip() for l in f if l.strip()]

# objects to save cut values
dcaCuts = []
thetaCuts = []
maxr = []
hardCurvCut = 0.0
phicuts = []
minz = []
maxz = []
cellZ0cut = []

i = 0
while i < len(lines):
    # Requires header starting with #
    if lines[i].startswith("#"):
        label = lines[i][1:].strip()
        i += 1
        if i + 1 >= len(lines):
            raise ValueError(f"Data missing for the plots '{label}'")

        # Reads y and x lines
        y = np.array([float(v) for v in lines[i].split(",") if v.strip()])
        x = np.array([float(v) for v in lines[i+1].split(",") if v.strip()])
        i += 2

        if len(x) != len(y):
            raise ValueError(f"x and y data have different sizes for plot '{label}'")

        # --- Interval calcuation ---
        order = np.argsort(x)
        x = x[order]
        y = y[order]

        cdf = np.cumsum(y)
        cdf = cdf / cdf[-1]

        percentageForCuts = args.percentageForCuts / 100.0
        modForSymmetricCuts = (1.0 - percentageForCuts) / 2.0

        # Interval depends on the name of the variable
        if "INNER Z" in label.upper():
            q_low, q_high = modForSymmetricCuts, percentageForCuts + modForSymmetricCuts
        elif "DPHI" in label.upper() and "IDPHI" not in label.upper():
            q_low, q_high = modForSymmetricCuts, percentageForCuts + modForSymmetricCuts
        elif "DZ" in label.upper():
            q_low, q_high = modForSymmetricCuts, percentageForCuts + modForSymmetricCuts
        else:
            q_low, q_high = 0.0, percentageForCuts

        low_idx = np.searchsorted(cdf, q_low)
        high_idx = np.searchsorted(cdf, q_high)
        x_low, x_high = x[low_idx], x[high_idx]

        if "DCA" in label.upper(): dcaCuts.append(x_high)
        if "THETA" in label.upper(): thetaCuts.append(x_high)
        if "DR" in label.upper(): maxr.append(x_high)
        if "HARDCURVCUT" in label.upper(): hardCurvCut = x_high
        if "IDPHI" in label.upper(): phicuts.append(x_high)
        if "INNER Z" in label.upper():
            minz.append(x_low - 1.0)
            maxz.append(x_high + 1.0)
        if "Z0" in label.upper(): cellZ0cut.append(x_high)
        

        if args.debug: print(f"📊 {label}: intervalo ({q_low:.3f}, {q_high:.3f}) → ({x_low:.4g}, {x_high:.4g})")

        if not args.summary_only:

            # --- Plotting with original binning ---
            plt.figure(figsize=(6, 6))
            bin_edges = np.empty(len(x) + 1)
            bin_edges[1:-1] = (x[:-1] + x[1:]) / 2
            bin_edges[0] = x[0] - (bin_edges[1] - x[0])
            bin_edges[-1] = x[-1] + (x[-1] - bin_edges[-2])

            plt.stairs(
                y,
                bin_edges,
                color='steelblue',
                alpha=0.7,
                fill=True,
                label=f"{label}\n[{q_low:.3f}, {q_high:.3f}] → [{x_low:.4f}, {x_high:.4f}]"
            )

            # Lines from interval
            plt.axvline(x_low, color='red', linestyle='dotted', label='Inferior limit')
            plt.axvline(x_high, color='green', linestyle='dashed', label='Superior limit')

            # Scale adjustement
            if "PT" in label.upper():
                plt.xscale("log")
            if "IPHI" in label.upper():
                plt.yscale("log")

            plt.xlabel(label)
            plt.ylabel("Entries")
            plt.legend()
            output_name = f"{outDir}/testSimPixelTrackHisto_{label.replace(' ', '_')}_{args.inpFileName}.pdf"
            if pdf is not None:
                pdf.savefig()
            else:
                output_name = (
                    f"{outDir}/testSimPixelTrackHisto_"
                    f"{label.replace(' ', '_')}_{args.inpFileName}.pdf"
                )
                plt.savefig(output_name)
            plt.close()
            if args.debug: print(f"✅ Plot created: {output_name}")
    else:
        raise ValueError(f"Unexpected line (expected '# NameOfPlot'): {lines[i]}")

if pdf is not None:
            pdf.close()

summary = "\n".join([
    "=======================================================================",
    f"\"phiCuts\": [{', '.join('{:g}'.format(value) for value in phicuts)}],\n",
    f"\"minZ\": {minz},\n",
    f"\"maxZ\": {maxz},\n",
    f"\"maxR\": {maxr},\n",
    f"\"thetaCuts\": {thetaCuts},\n",
    f"\"dcaCuts\": {dcaCuts}",
    "=======================================================================",
    f"\"hardCurvCut\": {hardCurvCut},",
    f"\"cellZ0cut\": {max(cellZ0cut)}",
    "=======================================================================",
])

print(summary)

summaryFile = os.path.join(outDir, f"simPixelTracks_summary_{args.inpFileName}.txt")
with open(summaryFile, "w", encoding="utf-8") as f:
    f.write(summary + "\n")
