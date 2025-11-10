#!/usr/bin/env python3
"""
TrackML → SimDoublets-style plots (cut-tuning aide)

This script ingests TrackML event data (CSV from the original TrackML release,
OR the SoA-style binaries produced by your "trackml-hits-particles-binaries.py"
helper) and produces a set of diagnostic plots analogous to those used to tune
seeding/doublet cuts (Δφ, Δz, Δη, ΔR, curvature proxy, d0 proxy, etc.),
with **true-only SimDoublets** support to emulate the CMSSW `SimDoublets` logic.

Outputs are PNGs in the chosen --outdir; a quick HTML index is generated.

Key ideas
---------
- **True-only (default):** For each TrackingParticle, collect its pixel hits, sort them by distance from the TP vertex,
  and form **only** the doublets between the first occurrence of a layer and the next different layer (pair inner hit
  to **all** hits in that next layer). This mirrors `SimDoublets::getSimDoublets()`.
- **Mixed mode:** Optionally, form all adjacent-layer pairs and tag signal/background using truth labels.
- Variables mirror common seeding studies: Δφ, Δz, Δη, ΔR, z0, and pT-from-radius (R) from the two hits + beamspot.
  (TrackML CSVs lack cluster-size info, so Y-size based cuts/plots are omitted.)

Usage examples
--------------
CSV input (TrackML official):
  python3 trackml_simdoublets_plotter.py \
    --indir data/trackml --pattern event000001000 \
    --outdir plots/evt1000 --volumes pixel strip_short \
    --mode true-only --max-pairs 2_000_000

Binary input (from your SoA builder):
  python3 trackml_simdoublets_plotter.py \
    --binaries data/trackml --pattern event000001000 \
    --outdir plots/evt1000-bin --volumes pixel strip_short --mode true-only

Notes on binaries
-----------------
We attempt to auto-detect the TRH1 / PAR1 layout written by your helper.
If the layout differs, use --binary-map to provide a small JSON describing the
column order; otherwise the script will raise a clear error explaining what's
missing.
"""

from __future__ import annotations
import argparse, json, math, os, textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------
# Volume selections (align with your existing geometry script)
# ------------------------------
VOLUME_SETS = {
    "pixel":       {7, 8, 9},
    "strip_short": {12, 13, 14},
    "strip_long":  {15, 16},
    "all":         set(range(1, 20)),
}

# ------------------------------
# Helpers
# ------------------------------

def _angle_wrap(dphi: np.ndarray) -> np.ndarray:
    return (dphi + np.pi) % (2*np.pi) - np.pi


def _safe_arctan2(y, x):
    return np.arctan2(y, x)


def _eta(x, y, z):
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(r, np.abs(z))  # handle z=0 robustly
    # Avoid infs at exactly z=0
    eta = -np.sign(z) * np.log(np.tan(0.5 * theta + 1e-12))
    eta[np.isinf(eta)] = 0.0
    return eta


def _circle_through_two_points_and_origin(x1, y1, x2, y2):
    """Return curvature proxy (1/R) and d0 proxy for the circle through (0,0), (x1,y1), (x2,y2).
    We use algebraic circle fit with one point fixed at origin -> robust & fast.
    If points are collinear with origin, returns zeros.
    """
    # Using 2x2 determinant-based solution for circle through three points.
    # Determinant D = 2 * |x1 y1 1; x2 y2 1; 0 0 1| = 2*(x1*y2 - x2*y1)
    D = 2.0 * (x1*y2 - x2*y1)
    eps = 1e-12
    mask = np.abs(D) > eps

    cx = np.zeros_like(D)
    cy = np.zeros_like(D)

    # S1 = x1^2 + y1^2; S2 = x2^2 + y2^2
    S1 = x1*x1 + y1*y1
    S2 = x2*x2 + y2*y2

    cx[mask] = ( (S1 * (y2 - 0.0)) - (S2 * (y1 - 0.0)) ) / D[mask]
    cy[mask] = ( (S2 * (x1 - 0.0)) - (S1 * (x2 - 0.0)) ) / D[mask]

    R = np.sqrt(cx*cx + cy*cy)
    k = np.zeros_like(R)
    k[mask] = 1.0 / (R[mask] + eps)

    # d0 proxy: signed distance of closest approach to origin = R - |C|
    # For a circle center C=(cx,cy) and origin (0,0), d0 = R - ||C||
    d0 = R - np.sqrt(cx*cx + cy*cy)
    d0[~mask] = 0.0

    return k, d0


@dataclass
class Inputs:
    hits: pd.DataFrame  # columns: hit_id, x, y, z, volume_id, layer_id, module_id, (optional) particle_id
    particles: Optional[pd.DataFrame]  # columns: particle_id, q, px, py, pz, vx, vy, vz, nhits


# ------------------------------
# Loaders
# ------------------------------

def load_csv(indir: Path, pattern: str, volumes: set[int]) -> Inputs:
    hits_path = indir / f"{pattern}-hits.csv"
    parts_path = indir / f"{pattern}-particles.csv"

    hits = pd.read_csv(hits_path)
    # standard TrackML columns: hit_id, x, y, z, volume_id, layer_id, module_id
    # sometimes truth map is in a separate file; if present, merge
    truth_path = indir / f"{pattern}-truth.csv"
    if truth_path.exists():
        truth = pd.read_csv(truth_path)[["hit_id", "particle_id"]]
        hits = hits.merge(truth, on="hit_id", how="left")
    else:
        hits["particle_id"] = -1

    # Filter volumes
    hits = hits[hits["volume_id"].isin(volumes)].copy()

    particles = None
    if parts_path.exists():
        particles = pd.read_csv(parts_path)
        if "q" not in particles.columns:
            # TrackML uses 'q' implicitly via particle_id sign; if missing, set q=1 for all
            particles["q"] = 1
    return Inputs(hits=hits, particles=particles)


def _read_block(f, count, dtype):
    arr = np.fromfile(f, dtype=dtype, count=count)
    if arr.size != count:
        raise RuntimeError("Unexpected EOF while reading binary block")
    return arr


def load_binaries(bindir: Path, pattern: str, binary_map_json: Optional[Path]) -> Inputs:
    """Read PAR1/TRH1/MAP1 files with minimal header checking.
    Expected filenames: {pattern}.particles.bin, {pattern}.hits.bin, {pattern}.map.bin
    """
    par = bindir / f"{pattern}.particles.bin"
    hit = bindir / f"{pattern}.hits.bin"
    mapp = bindir / f"{pattern}.map.bin"

    if not par.exists() or not hit.exists() or not mapp.exists():
        raise FileNotFoundError("Missing required binary files (.particles.bin/.hits.bin/.map.bin)")

    # Optional layout map
    layout = None
    if binary_map_json is not None:
        with open(binary_map_json) as jf:
            layout = json.load(jf)

    def read_header_and_check(f, magic: bytes):
        hdr = f.read(4)
        if hdr != magic:
            raise RuntimeError(f"Magic mismatch: expected {magic} got {hdr}")
        endian = int.from_bytes(f.read(4), "little")
        if endian != 0x01020304:
            raise RuntimeError("Endianness word mismatch (expect 0x01020304)")
        version = int.from_bytes(f.read(4), "little")
        if version != 1:
            raise RuntimeError(f"Unsupported format version {version}")

    # Particles
    with open(par, "rb") as f:
        read_header_and_check(f, b"PAR1")
        n = int.from_bytes(f.read(8), "little")
        # Default particle layout (match common TrackML order):
        # particle_id (int32), vx,vy,vz, px,py,pz (float32), q (int8 or int16), nhits (int16)
        if layout and "particles" in layout:
            dt = np.dtype(layout["particles"])  # e.g. a list of (name, <numpy-dtype>)
        else:
            dt = np.dtype([
                ("particle_id", np.int32),
                ("vx", np.float32), ("vy", np.float32), ("vz", np.float32),
                ("px", np.float32), ("py", np.float32), ("pz", np.float32),
                ("q", np.int16), ("nhits", np.int16),
            ])
        parr = np.fromfile(f, dtype=dt, count=n)
    particles = pd.DataFrame(parr)

    # Hits
    with open(hit, "rb") as f:
        read_header_and_check(f, b"TRH1")
        n = int.from_bytes(f.read(8), "little")
        if layout and "hits" in layout:
            dt = np.dtype(layout["hits"])  # must include hit_id, x,y,z, volume_id, layer_id, module_id
        else:
            dt = np.dtype([
                ("hit_id", np.int32),
                ("x", np.float32), ("y", np.float32), ("z", np.float32),
                ("volume_id", np.int16), ("layer_id", np.int16), ("module_id", np.int32),
            ])
        harr = np.fromfile(f, dtype=dt, count=n)
    hits = pd.DataFrame(harr)

    # Map: particle index per hit (same order as hits)
    with open(mapp, "rb") as f:
        read_header_and_check(f, b"MAP1")
        n = int.from_bytes(f.read(8), "little")
        idx = np.fromfile(f, dtype=np.int32, count=n)
    # Guard length
    if len(idx) != len(hits):
        raise RuntimeError("Hit→particle map length does not match number of hits")

    # particle_id through index
    if "particle_id" in particles.columns:
        pid = particles["particle_id"].to_numpy()
        particle_id = np.where((idx >= 0) & (idx < pid.size), pid[idx], -1)
    else:
        particle_id = idx  # already particle ids

    hits["particle_id"] = particle_id
    return Inputs(hits=hits, particles=particles)


# ------------------------------
# Doublet construction
# ------------------------------

def _quantized_iphi(phi: np.ndarray, nbits: int = 7) -> np.ndarray:
    # mimic approx_atan2s<7>: quantize phi symmetrically into integer bins
    scale = 1 << nbits  # 128 bins per pi
    return np.round(phi * (scale / np.pi)).astype(np.int32)


def make_true_doublets(hits: pd.DataFrame, particles: Optional[pd.DataFrame], max_pairs: int) -> pd.DataFrame:
    """Emulate SimDoublets: for each TP, form doublets only between first-occurrence layer and the next layer.

    Requires hits.particle_id (>=0). Uses particle vertex (vx,vy,vz) if available to sort hits; otherwise sorts by r.
    """
    if "particle_id" not in hits.columns:
        return pd.DataFrame(columns=["volume","li","lj","hit1","hit2","pid","dphi","dz","dr","innerZ","idphi","z0","R","pT","layerPairId","numSkippedLayers"])  

    # Build quick particle vertex lookup
    vmap = {}
    if particles is not None and {"particle_id","vx","vy","vz"}.issubset(particles.columns):
        vmap = particles.set_index("particle_id")[ ["vx","vy","vz"] ].to_dict("index")

    recs = []
    for pid, dfp in hits[hits.particle_id >= 0].groupby("particle_id"):
        if len(dfp) < 2:
            continue
        vx,vy,vz = (0.0,0.0,0.0)
        if pid in vmap:
            vv = vmap[pid]
            vx,vy,vz = vv["vx"], vv["vy"], vv["vz"]
        # positions relative to vertex
        x = dfp["x"].to_numpy(np.float64) - vx
        y = dfp["y"].to_numpy(np.float64) - vy
        z = dfp["z"].to_numpy(np.float64) - vz
        r = np.sqrt(x*x + y*y)
        phi = np.arctan2(y, x)
        iphi = _quantized_iphi(phi)
        layer = dfp["layer_id"].to_numpy()
        volume = dfp["volume_id"].to_numpy()
        hit_id = dfp["hit_id"].to_numpy()

        # sort by distance from vertex (mag^2)
        order = np.argsort(r*r + z*z)
        x,y,z,r,phi,iphi,layer,volume,hit_id = [arr[order] for arr in (x,y,z,r,phi,iphi,layer,volume,hit_id)]

        # scan like SimDoublets: for each i, find outer layer start j where layer[j] != layer[i], pair i with all j.. while same layer
        n = len(layer)
        for i in range(n):
            li = layer[i]
            # find next layer index
            j = i+1
            while j < n and layer[j] == li:
                j += 1
            if j >= n:
                break
            lj = layer[j]
            outer_start = j
            k = outer_start
            while k < n and layer[k] == lj:
                # features
                inner_r, inner_z, inner_phi, inner_iphi = r[i], z[i], phi[i], iphi[i]
                outer_r, outer_z, outer_phi, outer_iphi = r[k], z[k], phi[k], iphi[k]
                dr = outer_r - inner_r
                dz = outer_z - inner_z
                dphi = _angle_wrap(outer_phi - inner_phi)
                idphi = min(abs(int(outer_iphi - inner_iphi)), abs(int(inner_iphi - outer_iphi)))
                # z0 and R (curvature radius) like in analyzer
                eps = 1e-9
                z0 = abs(inner_r * outer_z - inner_z * outer_r) / max(abs(dr), eps)
                R = 0.5 * math.sqrt((dr / max(abs(dphi), eps))**2 + (inner_r * outer_r))
                pT = R / 87.78  # GeV (3.8T)
                layerPairId = int(li)*100 + int(lj)
                numSkipped = int(lj - li - 1)
                recs.append({
                    "volume": int(volume[i]),
                    "li": int(li), "lj": int(lj),
                    "hit1": int(hit_id[i]), "hit2": int(hit_id[k]),
                    "pid": int(pid),
                    "dr": dr, "dz": dz, "dphi": dphi, "idphi": idphi,
                    "innerZ": inner_z,
                    "z0": z0, "R": R, "pT": pT,
                    "layerPairId": layerPairId,
                    "numSkippedLayers": numSkipped,
                })
                k += 1
    if not recs:
        return pd.DataFrame(columns=["volume","li","lj","hit1","hit2","pid","dphi","dz","dr","innerZ","idphi","z0","R","pT","layerPairId","numSkippedLayers"])  

    df = pd.DataFrame.from_records(recs)
    if len(df) > max_pairs:
        df = df.sample(n=max_pairs, random_state=123, replace=False)
    return df


def make_mixed_doublets(hits: pd.DataFrame,
                  max_pairs: int,
                  pre_dz: float,
                  pre_dphi: float,
                  only_adjacent: bool = True) -> pd.DataFrame:
    """Form doublets between hits of consecutive layers inside each volume (signal/background tagging).

    Returns DataFrame with columns for indices and kinematic deltas, plus is_signal.
    """
    # Compute cylindrical/angles for convenience
    x = hits["x"].to_numpy(np.float64)
    y = hits["y"].to_numpy(np.float64)
    z = hits["z"].to_numpy(np.float64)

    r = np.sqrt(x*x + y*y)
    phi = _safe_arctan2(y, x)
    eta = _eta(x, y, z)

    base = hits[["hit_id", "volume_id", "layer_id", "module_id", "particle_id"]].copy()
    base["r"] = r
    base["phi"] = phi
    base["eta"] = eta
    base["z"] = z

    # Group by volume (keeps pixel vs strips separate if desired)
    groups = base.groupby("volume_id")

    recs = []
    for vol, dfv in groups:
        layers = sorted(dfv["layer_id"].unique())
        lay_pairs = zip(layers, layers[1:]) if only_adjacent else [
            (li, lj) for li in layers for lj in layers if lj > li
        ]
        for li, lj in lay_pairs:
            a = dfv[dfv.layer_id == li]
            b = dfv[dfv.layer_id == lj]
            if a.empty or b.empty:
                continue
            a_ = a.assign(phi_bin=(np.floor((a["phi"].to_numpy()+np.pi)/(pre_dphi+1e-9)).astype(int)))
            b_ = b.assign(phi_bin=(np.floor((b["phi"].to_numpy()+np.pi)/(pre_dphi+1e-9)).astype(int)))
            merged = a_.merge(b_, on="phi_bin", suffixes=("1","2"))
            if merged.empty:
                continue
            dphi = _angle_wrap(merged["phi2"].to_numpy() - merged["phi1"].to_numpy())
            dz   = (merged["z2"].to_numpy()   - merged["z1"].to_numpy())
            sel = (np.abs(dphi) <= pre_dphi) & (np.abs(dz) <= pre_dz)
            merged = merged.loc[sel]
            if merged.empty:
                continue

            # Compute features
            deta = merged["eta2"].to_numpy() - merged["eta1"].to_numpy()
            dr = merged["r2"].to_numpy() - merged["r1"].to_numpy()
            dR = np.sqrt(deta*deta + dphi[sel]*dphi[sel])
            k, d0 = _circle_through_two_points_and_origin(
                merged["r1"].to_numpy()*np.cos(merged["phi1"].to_numpy()),
                merged["r1"].to_numpy()*np.sin(merged["phi1"].to_numpy()),
                merged["r2"].to_numpy()*np.cos(merged["phi2"].to_numpy()),
                merged["r2"].to_numpy()*np.sin(merged["phi2"].to_numpy()),
            )
            rec = pd.DataFrame({
                "volume": vol,
                "li": li, "lj": lj,
                "hit1": merged["hit_id1"].to_numpy(),
                "hit2": merged["hit_id2"].to_numpy(),
                "pid1": merged["particle_id1"].to_numpy(),
                "pid2": merged["particle_id2"].to_numpy(),
                "dphi": dphi[sel],
                "dz": dz[sel],
                "deta": deta,
                "dr": dr,
                "dR": dR,
                "k": k,
                "d0": d0,
            })
            recs.append(rec)
    if not recs:
        return pd.DataFrame(columns=["volume","li","lj","hit1","hit2","pid1","pid2","dphi","dz","deta","dr","dR","k","d0"])  

    doublets = pd.concat(recs, ignore_index=True)
    # Random downsample if too many
    if len(doublets) > max_pairs:
        doublets = doublets.sample(n=max_pairs, random_state=123, replace=False)
    # Signal/background flag
    doublets["is_signal"] = (doublets["pid1"] >= 0) & (doublets["pid1"] == doublets["pid2"])
    return doublets


# ------------------------------
# Plotting
# ------------------------------

def _mkdir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def plot_hist(ax, data_s, data_b, title, xlabel, bins=120, logy=True):
    if data_s.size:
        ax.hist(data_s, bins=bins, histtype="step", label="signal", density=True)
    if data_b.size:
        ax.hist(data_b, bins=bins, histtype="step", label="background", density=True)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    if logy:
        ax.set_yscale("log")
    ax.legend()


def make_plots(hits: pd.DataFrame, doublets: pd.DataFrame, outdir: Path, mode: str):
    _mkdir(outdir)
    # 1) Per-layer hit multiplicity
    fig, ax = plt.subplots(figsize=(7,4))
    counts = hits.groupby(["volume_id","layer_id"]).size().rename("nhits").reset_index()
    labels = [f"V{v}-L{l}" for v,l in zip(counts.volume_id, counts.layer_id)]
    ax.bar(np.arange(len(labels)), counts.nhits)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_ylabel("# hits")
    ax.set_title("Per-layer hit multiplicity")
    fig.tight_layout(); fig.savefig(outdir/"hits_per_layer.png"); plt.close(fig)

    if doublets.empty:
        return

    # 2) Doublet kinematics (true-only has no background overlay)
    if mode == "true-only":
        s = doublets  # all entries are true SimDoublets
        plots = [
            (s.dphi.to_numpy(),    None, r"Δφ", "dphi", (-0.2,0.2)),
            (s.dz.to_numpy(),      None, r"Δz [cm]", "dz", None),
            (s.dr.to_numpy(),      None, r"dr [cm]", "dr", None),
            (s.innerZ.to_numpy() if "innerZ" in s.columns else s.dz.to_numpy()*0,  None, r"z_{inner} [cm]", "innerZ", None),
            (s.idphi.to_numpy() if "idphi" in s.columns else np.zeros(len(s)),   None, r"idφ (int)", "idphi", None),
            (s.z0.to_numpy()  if "z0" in s.columns else np.zeros(len(s)),        None, r"z0 [cm]", "z0", None),
            (s.pT.to_numpy()  if "pT" in s.columns else np.zeros(len(s)),        None, r"p_{T} from R [GeV]", "pT", None),
        ]
        for sdata, bdata, title, name, xlim in plots:
            fig, ax = plt.subplots(figsize=(6,4))
            ax.hist(sdata, bins=120, histtype="step", density=True, label="true doublets")
            ax.set_title(f"Doublets: {title}")
            ax.set_xlabel(title); ax.set_ylabel("density"); ax.set_yscale("log")
            if xlim: ax.set_xlim(*xlim)
            ax.legend(); fig.tight_layout(); fig.savefig(outdir[f"doublets_{name}.png"]); plt.close(fig)
    else:
        sig = doublets[doublets.is_signal]
        bkg = doublets[~doublets.is_signal]
        plots = [
            (sig.dphi.to_numpy(), bkg.dphi.to_numpy(), r"Δφ", "dphi", (-0.2,0.2)),
            (sig.dz.to_numpy(),   bkg.dz.to_numpy(),   r"Δz [cm]", "dz", None),
            (sig.deta.to_numpy(), bkg.deta.to_numpy(), r"Δη", "deta", None),
            (sig.dR.to_numpy(),   bkg.dR.to_numpy(),   r"ΔR", "dR", None),
        ]
        for sdata, bdata, title, name, xlim in plots:
            fig, ax = plt.subplots(figsize=(6,4))
            plot_hist(ax, sdata, bdata, f"Doublets: {title}", title)
            if xlim: ax.set_xlim(*xlim)
            fig.tight_layout(); fig.savefig(outdir/f"doublets_{name}.png"); plt.close(fig)

    # 3) Per-layer-pair yields
    fig, ax = plt.subplots(figsize=(7,4))
    if mode == "true-only":
        pair_counts = doublets.groupby(["volume","li","lj"]).size()
        lab = [f"V{v}-L{li}→L{lj}" for (v,li,lj) in pair_counts.index]
        X = np.arange(len(lab))
        ax.bar(X, pair_counts.values)
        ax.set_xticks(X); ax.set_xticklabels(lab, rotation=90)
        ax.set_ylabel("# true doublets")
        ax.set_title("True doublets per layer-pair")
    else:
        pair_counts = doublets.groupby(["volume","li","lj","is_signal"]).size().unstack(fill_value=0)
        lab = [f"V{v}-L{li}→L{lj}" for v,li,lj in pair_counts.index]
        y_sig = pair_counts.get(True, pd.Series(index=pair_counts.index, data=0))
        y_bkg = pair_counts.get(False, pd.Series(index=pair_counts.index, data=0))
        X = np.arange(len(lab))
        ax.bar(X-0.2, y_sig, width=0.4, label="signal")
        ax.bar(X+0.2, y_bkg, width=0.4, label="background")
        ax.set_xticks(X); ax.set_xticklabels(lab, rotation=90)
        ax.set_ylabel("# doublets")
        ax.set_title("Doublets per layer-pair")
        ax.legend()
    fig.tight_layout(); fig.savefig(outdir/"doublets_per_pair.png"); plt.close(fig)

    # 4) NEW: Per-layer-pair distributions (one PNG per pair per variable)
    #   Files named: doublets_pair_V{v}_L{li}-{lj}_{var}.png
    per_pair_vars = []
    if mode == "true-only":
        per_pair_vars = [
            ("dphi", r"Δφ", (-0.2,0.2)),
            ("dz",   r"Δz [cm]", None),
            ("dr",   r"dr [cm]", None),
            ("idphi", r"idφ (int)", None) if "idphi" in doublets.columns else None,
            ("innerZ", r"z_{inner} [cm]", None) if "innerZ" in doublets.columns else None,
            ("z0", r"z0 [cm]", None) if "z0" in doublets.columns else None,
            ("pT", r"p_{T} from R [GeV]", None) if "pT" in doublets.columns else None,
        ]
        per_pair_vars = [v for v in per_pair_vars if v is not None]
        for (v, dfpair) in doublets.groupby(["volume","li","lj"]):
            vol, li, lj = v
            for col, label, xlim in per_pair_vars:
                arr = dfpair[col].to_numpy()
                if arr.size == 0:
                    continue
                fig, ax = plt.subplots(figsize=(6,4))
                ax.hist(arr, bins=120, histtype="step", density=True, label=f"V{vol} L{li}→L{lj}")
                ax.set_title(f"{label} — pair V{vol} L{li}→L{lj}")
                ax.set_xlabel(label); ax.set_ylabel("density"); ax.set_yscale("log")
                if xlim: ax.set_xlim(*xlim)
                ax.legend(); fig.tight_layout()
                fig.savefig(outdir/f"doublets_pair_V{vol}_L{li}-{lj}_{col}.png"); plt.close(fig)
    else:
        per_pair_vars = [
            ("dphi", r"Δφ", (-0.2,0.2)),
            ("dz",   r"Δz [cm]", None),
            ("deta", r"Δη", None),
            ("dR",   r"ΔR", None),
        ]
        for (v, dfpair) in doublets.groupby(["volume","li","lj"]):
            vol, li, lj = v
            sig = dfpair[dfpair.is_signal]
            bkg = dfpair[~dfpair.is_signal]
            for col, label, xlim in per_pair_vars:
                s = sig[col].to_numpy(); b = bkg[col].to_numpy()
                if s.size == 0 and b.size == 0:
                    continue
                fig, ax = plt.subplots(figsize=(6,4))
                if s.size: ax.hist(s, bins=120, histtype="step", density=True, label="signal")
                if b.size: ax.hist(b, bins=120, histtype="step", density=True, label="background")
                ax.set_title(f"{label} — pair V{vol} L{li}→L{lj}")
                ax.set_xlabel(label); ax.set_ylabel("density"); ax.set_yscale("log")
                if xlim: ax.set_xlim(*xlim)
                ax.legend(); fig.tight_layout()
                fig.savefig(outdir/f"doublets_pair_V{vol}_L{li}-{lj}_{col}.png"); plt.close(fig)

    # 5) Quick HTML index (append per-pair section)
    imgs = [p.name for p in outdir.glob("*.png")]
    html = ["<html><body><h1>SimDoublets-style plots</h1>"]
    # Main images first
    for im in sorted(imgs):
        if not im.startswith("doublets_pair_"):
            html.append(f"<div><h3>{im}</h3><img src='{im}' width='900'></div>")
    # Then per-pair
    html.append("<hr><h2>Per-layer-pair distributions</h2>")
    for im in sorted([i for i in imgs if i.startswith("doublets_pair_")]):
        html.append(f"<div><h3>{im}</h3><img src='{im}' width='900'></div>")
    html.append("</body></html>")
    (outdir/"index.html").write_text("
".join(html))


# ------------------------------
# CLI
# ------------------------------

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--indir", type=Path, help="Directory with TrackML CSVs")
    src.add_argument("--binaries", type=Path, help="Directory with PAR1/TRH1/MAP1 binaries")

    ap.add_argument("--pattern", required=True, help="Event pattern, e.g. event000001000")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory for plots")

    ap.add_argument("--volumes", nargs="+", default=["pixel"],
                    choices=list(VOLUME_SETS.keys()),
                    help="Volume subsets to include")

    ap.add_argument("--include-vol-ids", nargs="*", type=int, default=[],
                    help="Additional volume IDs to include (overrides set)")

    ap.add_argument("--pre-dz", type=float, default=10.0, help="Preselection |Δz| in cm (mixed mode)")
    ap.add_argument("--pre-dphi", type=float, default=0.1, help="Preselection |Δφ| in rad (mixed mode)")
    ap.add_argument("--max-pairs", type=int, default=1_000_000, help="Cap number of doublets")
    ap.add_argument("--non-adjacent", action="store_true", help="Also allow non-adjacent layer pairs (mixed mode)")

    ap.add_argument("--binary-map", type=Path, default=None,
                    help="JSON describing binary dtype layouts if non-standard")

    ap.add_argument("--mode", choices=["true-only","mixed"], default="true-only",
                    help="true-only = build SimDoublets per TP; mixed = all pairs with S/B tagging")

    args = ap.parse_args()

    vols = set().union(*[VOLUME_SETS[v] for v in args.volumes])
    vols |= set(args.include_vol_ids)

    if args.indir:
        data = load_csv(args.indir, args.pattern, vols)
    else:
        data = load_binaries(args.binaries, args.pattern, args.binary_map)
        # Filter volumes after loading
        data.hits = data.hits[data.hits["volume_id"].isin(vols)].copy()

    if args.mode == "true-only":
        doublets = make_true_doublets(data.hits, data.particles, max_pairs=args.max_pairs)
    else:
        doublets = make_mixed_doublets(
            data.hits,
            max_pairs=args.max_pairs,
            pre_dz=args.pre_dz,
            pre_dphi=args.pre_dphi,
            only_adjacent=(not args.non_adjacent),
        )

    make_plots(data.hits, doublets, args.outdir, mode=args.mode)
    print(f"Wrote plots to {args.outdir}")


if __name__ == "__main__":
    main()


# ------------------------------
# Plotting
# ------------------------------

def _mkdir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def plot_hist(ax, data_s, data_b, title, xlabel, bins=120, logy=True):
    if data_s.size:
        ax.hist(data_s, bins=bins, histtype="step", label="signal", density=True)
    if data_b.size:
        ax.hist(data_b, bins=bins, histtype="step", label="background", density=True)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    if logy:
        ax.set_yscale("log")
    ax.legend()


def make_plots(hits: pd.DataFrame, doublets: pd.DataFrame, outdir: Path):
    _mkdir(outdir)
    # 1) Per-layer hit multiplicity
    fig, ax = plt.subplots(figsize=(7,4))
    counts = hits.groupby(["volume_id","layer_id"]).size().rename("nhits").reset_index()
    labels = [f"V{v}-L{l}" for v,l in zip(counts.volume_id, counts.layer_id)]
    ax.bar(np.arange(len(labels)), counts.nhits)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_ylabel("# hits")
    ax.set_title("Per-layer hit multiplicity")
    fig.tight_layout(); fig.savefig(outdir/"hits_per_layer.png"); plt.close(fig)

    if doublets.empty:
        return

    # 2) Doublet kinematics
    sig = doublets[doublets.is_signal]
    bkg = doublets[~doublets.is_signal]

    plots = [
        (sig.dphi.to_numpy(), bkg.dphi.to_numpy(), r"Δφ", "dphi", (-0.2,0.2)),
        (sig.dz.to_numpy(),   bkg.dz.to_numpy(),   r"Δz [cm]", "dz", None),
        (sig.deta.to_numpy(), bkg.deta.to_numpy(), r"Δη", "deta", None),
        (sig.dR.to_numpy(),   bkg.dR.to_numpy(),   r"ΔR", "dR", None),
        (sig.k.to_numpy(),    bkg.k.to_numpy(),    r"k ≈ 1/R [1/cm]", "k", None),
        (sig.d0.to_numpy(),   bkg.d0.to_numpy(),   r"d0 proxy [cm]", "d0", None),
    ]
    for s, b, title, name, xlim in plots:
        fig, ax = plt.subplots(figsize=(6,4))
        plot_hist(ax, s, b, f"Doublets: {title}", title)
        if xlim: ax.set_xlim(*xlim)
        fig.tight_layout(); fig.savefig(outdir/f"doublets_{name}.png"); plt.close(fig)

    # 3) Per-layer-pair yields
    fig, ax = plt.subplots(figsize=(7,4))
    pair_counts = doublets.groupby(["volume","li","lj","is_signal"]).size().unstack(fill_value=0)
    lab = [f"V{v}-L{li}→L{lj}" for v,li,lj in pair_counts.index]
    y_sig = pair_counts.get(True, pd.Series(index=pair_counts.index, data=0))
    y_bkg = pair_counts.get(False, pd.Series(index=pair_counts.index, data=0))
    X = np.arange(len(lab))
    ax.bar(X-0.2, y_sig, width=0.4, label="signal")
    ax.bar(X+0.2, y_bkg, width=0.4, label="background")
    ax.set_xticks(X); ax.set_xticklabels(lab, rotation=90)
    ax.set_ylabel("# doublets")
    ax.set_title("Doublets per layer-pair")
    ax.legend()
    fig.tight_layout(); fig.savefig(outdir/"doublets_per_pair.png"); plt.close(fig)

    # 4) Quick HTML index
    imgs = [p.name for p in outdir.glob("*.png")]
    html = ["<html><body><h1>SimDoublets-style plots</h1>"]
    for im in sorted(imgs):
        html.append(f"<div><h3>{im}</h3><img src='{im}' width='900'></div>")
    html.append("</body></html>")
    (outdir/"index.html").write_text("\n".join(html))


# ------------------------------
# CLI
# ------------------------------

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--indir", type=Path, help="Directory with TrackML CSVs")
    src.add_argument("--binaries", type=Path, help="Directory with PAR1/TRH1/MAP1 binaries")

    ap.add_argument("--pattern", required=True, help="Event pattern, e.g. event000001000")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory for plots")

    ap.add_argument("--volumes", nargs="+", default=["pixel"],
                    choices=list(VOLUME_SETS.keys()),
                    help="Volume subsets to include")

    ap.add_argument("--include-vol-ids", nargs="*", type=int, default=[],
                    help="Additional volume IDs to include (overrides set)")

    ap.add_argument("--pre-dz", type=float, default=10.0, help="Preselection |Δz| in cm")
    ap.add_argument("--pre-dphi", type=float, default=0.1, help="Preselection |Δφ| in rad")
    ap.add_argument("--max-pairs", type=int, default=1_000_000, help="Cap number of doublets")
    ap.add_argument("--non-adjacent", action="store_true", help="Also allow non-adjacent layer pairs")

    ap.add_argument("--binary-map", type=Path, default=None,
                    help="JSON describing binary dtype layouts if non-standard")

    args = ap.parse_args()

    vols = set().union(*[VOLUME_SETS[v] for v in args.volumes])
    vols |= set(args.include_vol_ids)

    if args.indir:
        data = load_csv(args.indir, args.pattern, vols)
    else:
        data = load_binaries(args.binaries, args.pattern, args.binary_map)
        # Filter volumes after loading
        data.hits = data.hits[data.hits["volume_id"].isin(vols)].copy()

    doublets = make_doublets(
        data.hits,
        max_pairs=args.max_pairs,
        pre_dz=args.pre_dz,
        pre_dphi=args.pre_dphi,
        only_adjacent=(not args.non_adjacent),
    )

    make_plots(data.hits, doublets, args.outdir)
    print(f"Wrote plots to {args.outdir}")


if __name__ == "__main__":
    main()

