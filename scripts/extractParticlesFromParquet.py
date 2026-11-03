import os

import numpy as np
import vector
import argparse

from datasets import load_dataset

parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

parser.add_argument('-o', '--outFileName', default='')
parser.add_argument('--outDir', default='.', help='Directory for output files (default: current directory).')
parser.add_argument('-m', '--isMuonsOnly', action='store_true')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-e', '--eventsForDebug', type=int, default=1)
parser.add_argument('-n', '--numberOfEvents', type=int, default=10)

args = parser.parse_args()
outDir = args.outDir or '.'
os.makedirs(outDir, exist_ok=True)

muonModifier = ''
if args.isMuonsOnly: muonModifier = 'MuonsOnly'

# Output filename
partFileName = os.path.join(outDir, 'particles' + args.outFileName + muonModifier + '.txt')

# Since the particles information is appended to the partFileName, it is better
# to remove it before executing the code
if os.path.exists(partFileName): os.remove(partFileName)

# Use this in streaming mode so that it doesn't directly download all files into local.
# It will make available some datasets, but all can be iterated via a for loop, as shown below
# More details in https://huggingface.co/datasets/OpenDataDetector/ColliderML-Release-1,
# https://huggingface.co/docs/datasets/v1.10.1/package_reference/main_classes.html#datasets.IterableDataset
# and https://huggingface.co/docs/datasets/stream

particles_iter_dataset = load_dataset(
    "OpenDataDetector/ColliderML-Release-1",
    data_dir="data/ttbar_pu200_particles",
    split="train",
    streaming=True,
    columns=["px","py","pz","vx","vy","vz","energy","mass","charge","particle_id","pdg_id"],
)

numberOfEvents = 0

for particles in particles_iter_dataset:

    ppx = np.asarray(particles["px"], dtype=np.float32)
    ppy = np.asarray(particles["py"], dtype=np.float32)
    ppz = np.asarray(particles["pz"], dtype=np.float32)
    pvx = np.asarray(particles["vx"], dtype=np.float32)
    pvy = np.asarray(particles["vy"], dtype=np.float32)
    pvz = np.asarray(particles["vz"], dtype=np.float32)
    pe = np.asarray(particles["energy"], dtype=np.float32)
    pm = np.asarray(particles["mass"], dtype=np.float32)
    pc = np.asarray(particles["charge"], dtype=np.int16)
    pid = np.asarray(particles["particle_id"], dtype=np.uint32)
    ppid = np.asarray(particles["pdg_id"], dtype=np.int32)

    if args.debug and numberOfEvents >= args.eventsForDebug:
        print(f"Created file {partFileName} WHILE DEBUGGING")
        break

    # Build kinematics from the original event, then combine all selection cuts, similar to what is done in CMS.
    p4 = vector.arr({"px": ppx, "py": ppy, "pz": ppz, "E": pe})

    keep = (
        (np.abs(pvz) < 300.0)
        & (np.hypot(pvx, pvy) < 25.0)
    )

    if args.isMuonsOnly:
        keep &= np.abs(ppid) == 13

    keep &= (p4.pt > 0.9)
    keep &= (np.abs(p4.eta) < 3.2)

    ppx, ppy, ppz, pvx, pvy, pvz, pe, pm, pc, pid, ppid = (
        values[keep]
        for values in (ppx, ppy, ppz, pvx, pvy, pvz, pe, pm, pc, pid, ppid)
    )

    p4_afterCuts = vector.arr({"px": ppx, "py": ppy, "pz": ppz, "E": pe})

    # Check to see if the lists exist. Even if 0 particles, count the event to match with hits file
    if ppx.size == 0:
        with open(partFileName, 'a') as fout:
            fout.write("particles:0\n")
        numberOfEvents += 1
        continue

    order = np.argsort(p4_afterCuts.pt, kind="stable")

    with open(partFileName, "a") as fout:
        fout.write(f"particles:{len(order)}\n")
        for idx in order:
            fout.write(
                f"{float(pvx[idx] / 10.0)},{float(pvy[idx] / 10.0)},{float(pvz[idx] / 10.0)},"
                f"{float(ppx[idx])},{float(ppy[idx])},{float(ppz[idx])},{float(pe[idx])},"
                f"{float(p4_afterCuts.pt[idx])},{float(p4_afterCuts.eta[idx])},"
                f"{float(p4_afterCuts.phi[idx])},{float(pm[idx])},"
                f"{int(pc[idx])},{int(ppid[idx])},{int(pid[idx])}\n"
            )

    if args.debug: print(numberOfEvents)

    numberOfEvents += 1

    if numberOfEvents == args.numberOfEvents: break

print(f"Created file {partFileName} with {numberOfEvents} events")