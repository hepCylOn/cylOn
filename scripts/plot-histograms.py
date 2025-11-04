#!/usr/bin/env python3

import re
import sys
import copy
import glob
import collections
import os

import matplotlib.pyplot as plt
import numpy as np

import argparse

class Plot:
    def __init__(self, xlabel, ylabel):
        self._xlabel = xlabel
        self._ylabel = ylabel

    def xlabel(self):
        return self._xlabel

    def ylabel(self):
        return self._ylabel

plots = dict(
    cluster_n = Plot("Number of clusters", "Events"),
    cluster_per_module_n = Plot("Number of clusters / module", "Modules"),
    digi_adc = Plot("Digi ADC", "Digis"),
    digi_n = Plot("Number of Digis", "Events"),
    hit_charge = Plot("Cluster charge", "Hits"),
    hit_gr = Plot("Hit global r (cm)", "Hits"),
    hit_gx = Plot("Hit global x (cm)", "Hits"),
    hit_gy = Plot("Hit global y (cm)", "Hits"),
    hit_gz = Plot("Hit global z (cm)", "Hits"),
    hit_lex = Plot("Hit local x error (cm)", "Hits"),
    hit_ley = Plot("Hit local y error (cm)", "Hits"),
    hit_lx = Plot("Hit local x (cm)", "Hits"),
    hit_ly = Plot("Hit local y (cm)", "Hits"),
    hit_n = Plot("Number of hits", "Events"),
    hit_sizex = Plot("Cluster size in local x (pixels)", "Hits"),
    hit_sizey = Plot("Cluster size in local y (pixels)", "Hits"),
    module_n = Plot("Number of modules with data", "Events"),
    track_chi2 = Plot("Track chi2/ndof", "Tracks"),
    track_eta = Plot("Track eta", "Tracks"),
    track_n = Plot("Number of tracks", "Events"),
    track_nhits = Plot("Number of hits / track", "Tracks"),
    track_phi = Plot("Track phi", "Tracks"),
    track_pt = Plot("Track pT (GeV)", "Tracks"),
    track_quality = Plot("Track quality", "Tracks"),
    track_tip = Plot("Track transverse IP (cm)", "Tracks"),
    track_tip_zoom = Plot("Track transverse IP (cm)", "Tracks"),
    track_zip = Plot("Track longitudinal IP (cm)", "Tracks"),
    track_zip_zoom = Plot("Track longitudinal IP (cm)", "Tracks"),
    vertex_chi2 = Plot("Vertex chi2", "Vertices"),
    vertex_n = Plot("Number of vertices", "Events"),
    vertex_ndof = Plot("Vertex ndof", "Vertices"),
    vertex_pt2 = Plot("Vertex sum(pT^2) (GeV^2)", "Vertices"),
    vertex_z = Plot("Vertex z (cm)", "Vertices"),
)

class Histo:
    def __init__(self, content):
        self._name = content[0]
        self._allBins = int(content[1])
        self._nbins = int(self._allBins-2)
        self._min = float(content[2])
        self._max = float(content[3])
        data = [int(x) for x in content[4:]]
        self._underflow = data[0]
        self._overflow = data[-1]
        self._data = data[1:-1]

        self._binWidth = (self._max-self._min) / self._nbins

    def __str__(self):
        ret = self._name
        if self._underflow != 0:
            ret += " underflow {}".format(self._underflow)
        if self._overflow != 0:
            ret += " overflow {}".format(self._overflow)
        return ret

    def name(self):
        return self._name

    def binWidth(self):
        return self._binWidth

    def bins(self):
        return np.arange(0, self._nbins) * self._binWidth + self._min
        #return np.arange(0, self._nbins) * self._binWidth

    def binX(self, i):
        return self._min + i*self._binWidth

    def values(self):
        return self._data

def only_tracks(name):
    return name.startswith("track_")

def only_vertices(name):
    return name.startswith("vertex_")

def only_trackvtx(name):
    return only_tracks(name) or only_vertices(name)

def ratioHisto(num, den):
    ret = copy.copy(num)
    def div(n, d):
        if d == 0:
            return 0
        return n/d
    ret._underflow = div(ret._underflow, den._underflow)
    ret._overflow = div(ret._overflow, den._overflow)
    for i, n in enumerate(num._data):
        ret._data[i] = div(n, den._data[i])
    return ret

def ratio(num, den):
    x = []
    y = []
    for i, n in enumerate(num._data):
        if den._data[i] != 0:
            x.append(num.binX(i))
            y.append(n/den._data[i])
    return (x, y)

def fillAxBar(ax, histo, label, log):
    #print(len(h.bins()), len(h.values()))
    #print(h.bins())
    #print(h.values())
    ax.bar(histo.bins(), histo.values(), width=histo.binWidth(), label=label, align="edge", log=log, alpha=0.7)

def fillAxDot(ax, data, label, log):
    ax.plot(data[0], data[1], ".", label=label)

def fillAxes(dataLabels, ax, fillAx=fillAxBar,
             title=None, custom_plt=None, legend=True,
             xlabel=None, ylabel=None,
             xlim=None, ylim=None,
             skipColors=0, skipColorIndices=[],
             log=False):
    if xlabel is None:
        xlabel = ""
    if ylabel is None:
        ylabel = "Events"
    ax.set(xlabel=xlabel, ylabel=ylabel)
    if title is not None:
        ax.set(title=title)
    ax.grid()

    for i in range(skipColors):
        ax.plot([], [])

    ind = 0
    skip = skipColorIndices[:]
    for data, label in dataLabels:
        while len(skip) > 0 and ind == skip[0]:
            ind += 1
            ax.plot([], [])
            del skip[0]
        fillAx(ax, data, label, log)
        ind += 1

    if ylim is None:
        ylim = dict(bottom=0)
    if len(ylim) > 0:
        ax.set_ylim(**ylim)
    if xlim is not None:
        ax.set_xlim(**xlim)
    if legend:
        ax.legend()
    if custom_plt is not None:
        custom_plt(ax)

def makePlot(histos, output=None, *args, **kwargs):
    fig, ax = plt.subplots(1, 1)
    fillAxes(histos, ax, *args, **kwargs)
    #fig.show()
    if output is not None:
        fig.savefig(output+".png")
        #fig.savefig(output+".pdf")
    plt.close(fig)

def makePlots(histoData, labels, outdir, histoFilter=None, **kwargs):
    for histoName, histos in histoData.items():
        if histoFilter is not None and not histoFilter(histoName):
            continue
        histoLabels = []
        for label in labels:
            print("{} {}".format(label, str(histos[label])))
            histoLabels.append( (histos[label], label) )
        p = plots[histoName]
        makePlot(histoLabels, output=outdir + "/" + histoName, xlabel=p.xlabel(), ylabel=p.ylabel(), **kwargs)

def makeRatioPlots(histoData, denomLabel, numLabels, outdir, histoFilter=None, **kwargs):
    for histoName, histos in histoData.items():
        if histoFilter is not None and not histoFilter(histoName):
            continue
        print("ratio for {}".format(histoName))
        denom = histos[denomLabel]
        ratioHistoLabels = [ (ratio(histos[num], denom), num+"/"+denomLabel) for num in numLabels]
        p = plots[histoName]
        makePlot(ratioHistoLabels, output=outdir + "/" + histoName + "_ratio", fillAx=fillAxDot, xlabel=p.xlabel(), ylabel=p.ylabel(), skipColors=1, **kwargs)

def makeManyRatioPlots(histoData, denomLabel, numLabels, **kwargs):
    for histoName, histos in histoData.items():
        print("ratio for {}".format(histoName))
        denom = histos[denomLabel]
        ratioHistoLabels = []
        for numLabel, numFiles in numLabels.items():
            x = []
            y = []
            for n in numFiles:
                tx, ty = ratio(histos[n], denom)
                x.extend(tx)
                y.extend(ty)
            #print(x, y)
            ratioHistoLabels.append( ((x, y), numLabel+"/"+denomLabel) )
        makePlot(ratioHistoLabels, output=histoName+"_ratiomany", fillAx=fillAxDot, **kwargs)

def readHistograms(files):
    label_re = re.compile("histograms_(?P<label>.*)\.txt")
    histos = collections.defaultdict(dict)
    for fname in files:
        m = label_re.search(fname)
        label = m.group("label")
        with open(fname) as f:
            #print(fname)
            for line in f:
                histo = Histo(line.split(" "))
                histos[histo.name()][label] = histo
                #break
    return histos


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot histograms from text files. "
                    "First file is the reference for ratio plots."
    )
    parser.add_argument(
        "files",
        metavar="histogram_file",
        nargs="+",
        help="Histogram text files in the format 'histograms_<label>.txt'. "
             "First one is the reference."
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Use log scale for the y-axis in plots."
    )

    parser.add_argument(
        "--no-local-reco",
        action="store_true",
        help="Plot only track and vertices plots."
    )

    parser.add_argument(
        "--ylim",
        type=float,
        nargs=2,
        metavar=("YMIN", "YMAX"),
        help="Set y-axis limits for ratio plots."
    )

    parser.add_argument(
        "--output-dir",
        default="plots",
        help="Directory where plots will be saved (default: plots/)"
    )

    parser.add_argument(
        "--filter",
        help="Only plot histograms whose names match this regex."
    )
    
    args = parser.parse_args()

    if args.filter and args.no_local_reco:
        parser.error("Options --filter and --no-local-reco cannot be used together.")

    os.makedirs(args.output_dir, exist_ok=True)

    histoData = readHistograms(args.files)

    label_re = re.compile(r"histograms_(?P<label>.*)\.txt")
    labels = []
    for f in args.files:
        m = label_re.search(f)
        if not m:
            parser.error(f"{f} does not match expected format 'histograms_<label>.txt'")
        labels.append(m.group("label"))

    ref_label = labels[0]
    other_labels = labels[1:]

    histoFilter = None if not args.no_local_reco else only_trackvtx
    if args.filter:
        pattern = re.compile(args.filter)
        histoFilter = lambda name: bool(pattern.search(name))
        
    # Make plots for all
    makePlots(histoData, labels, log=args.log, outdir = args.output_dir, histoFilter = histoFilter)

    ylim_dict = dict()
    if args.ylim:
        ylim_dict = dict(ymin=args.ylim[0], ymax=args.ylim[1])

    makeRatioPlots(histoData, ref_label, other_labels, ylim=ylim_dict, outdir = args.output_dir, histoFilter = histoFilter)

# histoData = readHistograms(glob.glob("histograms_*.txt"))
# makePlots(histoData, ["baseline", "baseline_2", "0dd4e43","5a514d0"], log=True)
# makeRatioPlots(histoData, "baseline", ["baseline_2", "0dd4e43","5a514d0"], ylim=dict())


#makeManyRatioPlots(histoData, "cuda", dict(
#    cuda = ["cuda_{}".format(i) for i in range(0,100)],
#    kokkos_cuda = ["kokkos_cuda_{}".format(i) for i in range(0,100)],
#), ylim=dict(ymin=0.8, ymax=1.2))
