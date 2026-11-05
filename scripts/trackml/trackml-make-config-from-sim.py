import argparse, csv, json
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--template", default="", help="Optional base JSON to merge into")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = {}
    if args.template and Path(args.template).exists():
        cfg = json.loads(Path(args.template).read_text())
    if "producer" not in cfg:
        cfg["producer"] = {}
    if "deviceAlgo" not in cfg:
        cfg["deviceAlgo"] = {}

    cfg["producer"]["doStats"] = True

    rows = []
    with open(args.csv, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)

    rows.sort(key=lambda r: int(r["layerPairId"]))

    layerPairs = []
    cellMinz, cellMaxz, cellMaxr, cellPhiCuts = [], [], [], []
    for row in rows:
        lp = int(row["layerPairId"])
        inner = lp // 100
        outer = lp % 100
        layerPairs.extend([inner, outer])
        cellMinz.append(float(row["minZ"]))
        cellMaxz.append(float(row["maxZ"]))
        cellMaxr.append(float(row["maxR"]))
        cellPhiCuts.append(int(row["phiCut"]))

    cfg["deviceAlgo"]["layerPairs"] = layerPairs
    cfg["deviceAlgo"]["cellMinz"] = cellMinz
    cfg["deviceAlgo"]["cellMaxz"] = cellMaxz
    cfg["deviceAlgo"]["cellMaxr"] = cellMaxr
    cfg["deviceAlgo"]["cellPhiCuts"] = cellPhiCuts

    def pick_conservative(key, cast):
        vals = [cast(r[key]) for r in rows if r.get(key,"")!=""]
        if not vals: return None
        if key in ("z0Cut",):
            return max(vals)
        if key in ("ptCut",):
            return min(vals)
        if key in ("minYsizeB1","minYsizeB2"):
            return min(vals)
        if key in ("maxDYsize12","maxDYsize","maxDYPred"):
            return min(vals)
        return vals[0]

    dev = cfg["deviceAlgo"]
    for key, cast in [
        ("z0Cut", float),
        ("ptCut", float),
        ("minYsizeB1", int),
        ("minYsizeB2", int),
        ("maxDYsize12", int),
        ("maxDYsize", int),
        ("maxDYPred", int),
    ]:
        val = pick_conservative(key, cast)
        if val is not None:
            dev[key] = val

    prod = cfg["producer"]
    prod.setdefault("maxNumberOfDoublets", 5_000_000)
    prod.setdefault("maxNumberOfTuples",   500_000)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cfg, indent=2))
    print(f"[config] wrote {out} with {len(rows)} layer-pair rows.")

if __name__ == "__main__":
    main()