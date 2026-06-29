import matplotlib.pyplot as plt
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

args = parser.parse_args()

inputFile = 'output' + args.inpFileName + '.txt'
outDir = args.outDirName

os.mkdir(outDir)

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
cellZ0cut = 0.0

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

        # # Interval depends on the name of the variable
        # if "INNER Z" in label.upper():
        #     q_low, q_high = 0.005, 0.995
        # elif "DPHI" in label.upper() and "IDPHI" not in label.upper():
        #     q_low, q_high = 0.005, 0.995
        # elif "DZ" in label.upper():
        #     q_low, q_high = 0.005, 0.995
        # else:
        #     q_low, q_high = 0.0, 0.99

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
        if "Z0" in label.upper(): cellZ0cut = x_high
        

        if args.debug: print(f"📊 {label}: intervalo ({q_low:.3f}, {q_high:.3f}) → ({x_low:.4g}, {x_high:.4g})")

        # --- Plotting with original binning ---
        plt.figure(figsize=(6, 6))
        plt.hist(
            x,
            bins=len(x),
            weights=y,
            color='steelblue',
            alpha=0.7,
            label=f"{label}\n[{q_low:.3f}, {q_high:.3f}] → [{x_low:.4f}, {x_high:.4f}]"
        )

        # Lines from interval
        plt.axvline(x_low, color='red', linestyle='--', label='Inferior limit')
        plt.axvline(x_high, color='red', linestyle='--', label='Superior limit')

        # Scale adjustement
        if "PT" in label.upper():
            plt.xscale("log")
        if "IPHI" in label.upper():
            plt.yscale("log")

        plt.xlabel(label)
        plt.ylabel("Entries")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        output_name = f"{outDir}/testSimPixelTrackHisto_{label.replace(' ', '_')}_{args.inpFileName}.pdf"
        plt.savefig(output_name)
        plt.close()
        if args.debug: print(f"✅ Plot created: {output_name}")
    else:
        raise ValueError(f"Unexpected line (expected '# NameOfPlot'): {lines[i]}")

print("dcaCuts: ",dcaCuts)
print("=======================================================================")
print("thetaCuts: ",thetaCuts)
print("=======================================================================")
print("maxr: ",maxr)
print("=======================================================================")
print("hardCurvCut: ",hardCurvCut)
print("=======================================================================")
print("phicuts: ",phicuts)
print("=======================================================================")
print("minz: ",minz)
print("=======================================================================")
print("maxz: ",maxz)
print("=======================================================================")
# print("cellZ0cut: ",cellZ0cut)
