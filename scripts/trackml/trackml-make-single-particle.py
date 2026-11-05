#!/usr/bin/env python3
"""
Skim TrackML events to a single particle (and its hits), per event.
Converts spatial coordinates read from CSVs (vx, vy, vz, x, y, z) from mm -> cm.
Binary files are assumed already consistent; no conversion is applied to bins.

Default cuts:
  - pt in (10, 100)
  - sqrt(vx^2 + vy^2) <= 0.05   # in cm after conversion
  - q != 0
  - >= 3 hits in volumes {7, 8, 9}

Species selection:
  - If particles.csv has a PDG-like column (pid/pdgId/...), use that to enforce species.
  - Otherwise, only the charge sign is enforced for anti- vs particle.

Binaries:
  - Optionally call an external builder (e.g. trackml-hits-particles-binaries.py) once after
    skimming is done. If --pattern is provided, it is forwarded; otherwise the builder is called
    with no --pattern and it will discover files under --outdir.
"""

import argparse
from pathlib import Path
import sys
import subprocess
import pandas as pd
import numpy as np
import re
from typing import Optional, Tuple, Dict, List, Set

SPECIES_TO_PDG: Dict[str, int] = {
    "electron": 11, "e-": 11, "anti-electron": -11, "positron": -11, "e+": -11,
    "muon": 13, "mu-": 13, "anti-muon": -13, "mu+": -13,
    "pion": 211, "pi+": 211, "anti-pion": -211, "pi-": -211,
}

MM_TO_CM = 0.1  # apply only when reading CSVs

def parse_species(spec: Optional[str]) -> Optional[int]:
    if not spec:
        return None
    key = spec.strip().lower()
    if key not in SPECIES_TO_PDG:
        raise ValueError(f"Unknown species '{spec}'. Valid: {', '.join(sorted(set(SPECIES_TO_PDG.keys())))}")
    return SPECIES_TO_PDG[key]

def detect_pdg_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("pid", "pdgId", "pdgid", "pdg_id", "pdg"):
        if col in df.columns:
            return col
    return None

def infer_charge_sign_from_species(spec_pdg: Optional[int]) -> Optional[int]:
    if spec_pdg is None:
        return None
    return -1 if spec_pdg < 0 else +1

def parse_int_set(s: str) -> Set[int]:
    toks = re.split(r"[,\s]+", s.strip())
    return {int(t) for t in toks if t}

def event_triplets(inputdir: Path, patterns: List[str]) -> List[Tuple[str, Path, Path, Path]]:
    out = []
    for pat in patterns:
        for p in sorted(inputdir.glob(f"{pat}-particles.csv")):
            event_id = p.name.replace("-particles.csv", "")
            t = inputdir / f"{event_id}-truth.csv"
            h = inputdir / f"{event_id}-hits.csv"
            if t.exists() and h.exists():
                out.append((event_id, p, t, h))
    return out

def pixel_hit_counts_per_particle(truth_df: pd.DataFrame, hits_df: pd.DataFrame, pixel_vols: Set[int]) -> pd.Series:
    # Merge truth with hits to get volumes per (particle,hit)
    m = truth_df[["particle_id", "hit_id"]].merge(
        hits_df[["hit_id", "volume_id"]], on="hit_id", how="inner"
    )
    m = m[m["volume_id"].isin(pixel_vols)]
    if m.empty:
        return pd.Series(dtype=int)
    cnt = m.groupby("particle_id").size()
    cnt.name = "pixel_hits"
    return cnt

def select_particle(particles: pd.DataFrame,
                    truth: pd.DataFrame,
                    hits: pd.DataFrame,
                    pt_min: float,
                    pt_max: float,
                    vtx_max_cm: float,
                    species_pdg: Optional[int],
                    pixel_vols: Set[int],
                    min_pixel_hits: int) -> pd.Series:
    """Return a single row (Series) for the chosen particle satisfying all cuts.
       NOTE: particles (vx,vy,vz) and hits (x,y,z) have been converted to cm already.
    """
    particles = particles.copy()
    particles["pt"] = np.sqrt(particles["px"]**2 + particles["py"]**2)  # momenta unchanged
    particles["rt"] = np.sqrt(particles["vx"]**2 + particles["vy"]**2)  # cm

    # Attach per-particle pixel hit counts
    pix_counts = pixel_hit_counts_per_particle(truth, hits, pixel_vols)
    particles = particles.merge(pix_counts, how="left", left_on="particle_id", right_index=True)
    particles["pixel_hits"] = particles["pixel_hits"].fillna(0).astype(int)

    pdg_col = detect_pdg_column(particles)
    charge_needed = infer_charge_sign_from_species(species_pdg)

    mask = (
        (particles["pt"] > pt_min) &
        (particles["pt"] < pt_max) &
        (particles["rt"] <= vtx_max_cm) &
        (particles["q"] != 0) &
        (particles["pixel_hits"] >= min_pixel_hits)
    )
    if pdg_col and species_pdg is not None:
        mask &= (particles[pdg_col] == species_pdg)
    elif species_pdg is not None:
        mask &= (np.sign(particles["q"]) == charge_needed)

    cands = particles[mask]
    if cands.empty:
        raise RuntimeError("No particle matches the selection (including pixel hit requirement).")
    # Prefer most hits overall, then most pixel hits, then highest pt
    cands = cands.sort_values(["nhits", "pixel_hits", "pt"], ascending=[False, False, False])
    return cands.iloc[0]

def write_skim(event_id: str,
               outdir: Path,
               particles_row: pd.Series,
               truth_df: pd.DataFrame,
               hits_df: pd.DataFrame) -> Tuple[Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    pid = int(particles_row["particle_id"])

    p_out = outdir / f"{event_id}-particles.csv"
    t_out = outdir / f"{event_id}-truth.csv"
    h_out = outdir / f"{event_id}-hits.csv"

    particles_row.to_frame().T.to_csv(p_out, index=False)

    truth_sel = truth_df[truth_df["particle_id"] == pid].copy()
    truth_sel.to_csv(t_out, index=False)

    hit_ids = set(truth_sel["hit_id"].tolist())
    hits_sel = hits_df[hits_df["hit_id"].isin(hit_ids)].copy()
    hits_sel.to_csv(h_out, index=False)

    return p_out, t_out, h_out

def run_binary_builder(builder_path: Path,
                       indir: Path,
                       outdir: Path,
                       modules_csv: Optional[Path],
                       pattern: Optional[str],
                       extra_args: List[str]) -> None:

    cmd = [sys.executable, str(builder_path),
           "--indir", str(indir),
           "--outdir", str(outdir)]
    if pattern is not None:
        cmd += ["--pattern", pattern]
    if modules_csv:
        cmd += ["--modules", str(modules_csv)]
    cmd += extra_args
    print(f"[binaries] running: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)

def main():
    ap = argparse.ArgumentParser(description="Skim TrackML events to a single particle and its hits; optionally build binaries. Converts CSV spatial coords mm -> cm on read.")
    ap.add_argument("--inputdir", default="data/trackml/", help="Directory containing TrackML CSVs.")
    ap.add_argument("--outdir", default=None, help="Output directory. Default: data/trackml/single_{species}/")
    ap.add_argument("--pattern", default="", help="Comma-separated list or glob(s) like 'event000001005,event000001006' or 'event0000010*'. If empty, processes all 'event*-particles.csv'.")
    ap.add_argument("--species", default=None,
                    help="One of: electron, anti-electron (positron), muon, anti-muon, pion, anti-pion.")
    ap.add_argument("--pt-min", type=float, default=10.0)
    ap.add_argument("--pt-max", type=float, default=100.0)
    ap.add_argument("--vtx-max", type=float, default=0.05, help="Transverse vertex cut on sqrt(vx^2+vy^2), in cm (CSV values converted mm->cm on read).")
    ap.add_argument("--pixel-volumes", default="7,8,9", help="Comma/space-separated list of volume_ids considered 'pixel'.")
    ap.add_argument("--min-pixel-hits", type=int, default=3, help="Require at least this many hits in the pixel volumes.")
    ap.add_argument("--make-binaries", action="store_true", help="If set, invoke your binary builder on the skimmed files.")
    ap.add_argument("--binary-builder", default="trackml-hits-particles-binaries.py",
                    help="Path to your existing binary builder script.")
    ap.add_argument("--modules", default=None, help="Path to detectors.csv (required by your binary builder).")
    ap.add_argument("--builder-extra", default="", help="Extra args for the binary builder (e.g. '--volumes pixel strip_short').")
    args = ap.parse_args()

    inputdir = Path(args.inputdir).resolve()
    if not inputdir.exists():
        ap.error(f"Input dir does not exist: {inputdir}")

    species_pdg = parse_species(args.species) if args.species else None
    species_name = (args.species or "any").replace(" ", "_").replace("-", "")
    outdir = Path(args.inputdir + f"/single_{species_name}/")
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    pixel_vols = parse_int_set(args.pixel_volumes)

    if args.pattern.strip():
        raw_patterns = [s.strip() for s in args.pattern.split(",") if s.strip()]
    else:
        raw_patterns = ["event*"]

    events = event_triplets(inputdir, raw_patterns)
    if not events:
        ap.error(f"No events found under {inputdir} matching {raw_patterns}")

    print(f"[skim] Input: {inputdir}")
    print(f"[skim] Output: {outdir}")
    print(f"[skim] Events to process: {len(events)}")
    print(f"[skim] Pixel volumes: {sorted(pixel_vols)} (min hits: {args.min_pixel_hits})")
    print(f"[skim] NOTE: Converting CSV spatial coords mm -> cm on read; outputs will be in cm.")

    selected_summary = []
    for event_id, p_csv, t_csv, h_csv in events:
        particles = pd.read_csv(p_csv)
        truth = pd.read_csv(t_csv)
        hits = pd.read_csv(h_csv)

        # ---- convert spatial coordinates mm -> cm on read ----
        for c in ("vx","vy","vz"):
            if c in particles.columns:
                particles[c] = particles[c] * MM_TO_CM
        for c in ("x","y","z"):
            if c in hits.columns:
                hits[c] = hits[c] * MM_TO_CM
        # ------------------------------------------------------

        try:
            row = select_particle(particles, truth, hits,
                                  args.pt_min, args.pt_max, args.vtx_max,
                                  species_pdg, pixel_vols, args.min_pixel_hits)
        except Exception as e:
            print(f"[skip] {event_id}: {e}")
            continue

        p_out, t_out, h_out = write_skim(event_id, outdir, row, truth, hits)

        sel = {
            "event": event_id,
            "particle_id": int(row["particle_id"]),
            "pt": float(row["pt"]),
            "q": int(row["q"]),
            "vx": float(row["vx"]),
            "vy": float(row["vy"]),
            "vz": float(row["vz"]),
            "rt": float(np.sqrt(row["vx"]**2 + row["vy"]**2)),
            "nhits": int(row["nhits"]),
            "pixel_hits": int(row.get("pixel_hits", 0)),
            "particles_csv": str(p_out),
            "truth_csv": str(t_out),
            "hits_csv": str(h_out),
            "units": "cm"
        }
        pdg_col = detect_pdg_column(particles)
        if pdg_col and pdg_col in row:
            sel["pdgId"] = int(row[pdg_col])
        selected_summary.append(sel)

        print(f"[ok] {event_id}: particle_id={sel['particle_id']} pt={sel['pt']:.2f} q={sel['q']} nhits={sel['nhits']} pixel_hits={sel['pixel_hits']}")

    # Optionally build binaries once at the end
    if args.make_binaries:
        builder = Path(args.binary_builder).resolve()
        if not builder.exists():
            ap.error(f"Binary builder not found: {builder}")
        modules = Path(args.modules).resolve() if args.modules else None
        if args.modules and not modules.exists():
            ap.error(f"Modules file not found: {modules}")
        pattern = args.pattern if args.pattern else None
        extra_args = [tok for tok in args.builder_extra.split() if tok]
        run_binary_builder(builder, outdir, outdir, modules, pattern, extra_args)

    if selected_summary:
        manifest = pd.DataFrame(selected_summary)
        manifest_path = outdir / "skim_manifest.csv"
        manifest.to_csv(manifest_path, index=False)
        print(f"[skim] Wrote manifest: {manifest_path}")
    else:
        print("[skim] Nothing selected. Check your cuts/species and pixel-hit requirement.]")

if __name__ == "__main__":
    main()
