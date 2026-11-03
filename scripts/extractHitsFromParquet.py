import os

import numpy as np
import math
import argparse

from datasets import load_dataset

import matplotlib.pyplot as plt

def writeToFile(inputData,barrelPos,barrelThres,endcapPos,endcapThres,modules):

    getInitialModulesLen = len(modules)

    hxAux, hyAux, hzAux, hrAux, hpAux, hpartAux = inputData[1:]

    # Only keep hits that are in the current part of detector being checked
    keep = (
        (np.abs(hzAux) < endcapPos[-1] + endcapThres)
        & (hrAux < barrelPos[-1] + barrelThres)
        & (hrAux > barrelPos[0] - barrelThres)
    )

    hxAux = hxAux[keep]
    hyAux = hyAux[keep]
    hzAux = hzAux[keep]
    hrAux = hrAux[keep]
    hpAux = hpAux[keep]
    hpartAux = hpartAux[keep]

    # Initialize the layersID array for all the pixel hits
    saveLayersIDarray = getInitialModulesLen*np.ones(len(hzAux)) - np.ones(len(hzAux))

    abs_z = np.abs(hzAux)
    abs_r = np.abs(hrAux)
    endcap_boundary = endcapPos[0] - endcapThres

    barrel_region = abs_z < endcap_boundary
    pos_endcap_region = hzAux > endcap_boundary
    neg_endcap_region = hzAux < -endcap_boundary

    # Assign layer IDs and count hits per layer.
    for r_idx, radius in enumerate(barrelPos):
        in_layer = (
            barrel_region
            & (abs_r > radius - barrelThres)
            & (abs_r < radius + barrelThres)
        )
        saveLayersIDarray += r_idx * in_layer
        modules.append(np.count_nonzero(in_layer))

    # Positive endcap layers: append all of these counts first.
    for z_idx, z_layer in enumerate(endcapPos):
        in_layer = (
            pos_endcap_region
            & (hzAux > z_layer - endcapThres)
            & (hzAux < z_layer + endcapThres)
        )
        layer_id = z_idx + len(barrelPos)
        saveLayersIDarray += layer_id * in_layer
        modules.append(np.count_nonzero(in_layer))

    # Negative endcap layers: append these after all positive-endcap counts.
    for z_idx, z_layer in enumerate(endcapPos):
        in_layer = (
            neg_endcap_region
            & (hzAux > -(z_layer + endcapThres))
            & (hzAux < -(z_layer - endcapThres))
        )
        layer_id = z_idx + len(barrelPos) + len(endcapPos)
        saveLayersIDarray += layer_id * in_layer
        modules.append(np.count_nonzero(in_layer))

    # Convert hit counts per layer to cumulative counts.
    for i in range(getInitialModulesLen, len(modules)):
        modules[i] += modules[i - 1]

    # Order input data by the LayersID
    layer_ids = saveLayersIDarray.astype("i")
    order = np.argsort(layer_ids, kind="stable")

    return (
        hxAux[order] / 10.0,
        hyAux[order] / 10.0,
        hzAux[order] / 10.0,
        hrAux[order] / 10.0,
        np.arctan2(hyAux[order], hxAux[order]),
        layer_ids[order],
        hpartAux[order],
        modules,
    )

parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

parser.add_argument('-o', '--outFileName', default='')
parser.add_argument('-p', '--isPhase2', action='store_true')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-e', '--eventsForDebug', type=int, default=1)
parser.add_argument('-n', '--numberOfEvents', type=int, default=10)
parser.add_argument('-m', '--makeHitPosPlot', action='store_true')
parser.add_argument('-s', '--shortStrips', action='store_true')
parser.add_argument('-l', '--longStrips', action='store_true')
parser.add_argument('--skipHeaders', action='store_true')
parser.add_argument('--skipModules', action='store_true')
parser.add_argument('-a', '--forAutograph', action='store_true')

args = parser.parse_args()

def phi2short(x: float) -> int:
    p2i = (2**15) / math.pi  # 32768 / pi
    val = round(x * p2i)

    # simular wraparound de int16
    if val > 32767:
        val -= 65536
    elif val < -32768:
        val += 65536

    return val

phaseModifier = 'Phase1'
if args.isPhase2: phaseModifier = 'Phase2'

autographModifier = ''
if args.forAutograph: autographModifier = 'OneEventForAutograph'

# Autograph should only modify the hits file, but adding modifier to maps file if needed
hitsFileName = 'hitsWithParticleID' + args.outFileName + phaseModifier + autographModifier + '.txt'
mapFileName = 'mapHitsToParticles' + args.outFileName + phaseModifier + autographModifier + '.txt'

if os.path.exists(hitsFileName): os.system('rm ' + hitsFileName)
if os.path.exists(mapFileName): os.system('rm ' + mapFileName)

# These define the position of the layers in the open data detector geometry
# taken from https://iopscience.iop.org/article/10.1088/1742-6596/2438/1/012110/pdf
# Thresholds are set to get hits in all of a given layer, but not in other layers
colliderMLPixelBarrel = [34.0,70.0,116.0,172.0]
colliderMLPixelEndcap = [620.0,730.0,830.0]
if args.isPhase2: colliderMLPixelEndcap = [620.0,730.0,830.0,980.0,1120.0,1320.0,1520.0]
colliderMLPixelBarrelThreshold = 15.0
colliderMLPixelEndcapThreshold = 50.0

coliderMLShortStripsBarrel = [250.0,350.0,500.0,650.0]
coliderMLShortStripsEndcap = [1300.0,1580.0,1820.0,2200.0,2580.0,2980.0]
coliderMLShortStripsBarrelThreshold = 60.0
coliderMLShortStripsEndcapThreshold = 80.0

coliderMLLongStripsBarrel = [820.0,1020.0]
coliderMLLongStripsEndcap = [1300.0,1580.0,1900.0,2250.0,2600.0,3000.0]
coliderMLLongStripsBarrelThreshold = 60.0
coliderMLLongStripsEndcapThreshold = 100.0

# Use this in streaming mode so that it doesn't directly download all files into local.
# It will make available some datasets, but all can be iterated via a for loop, as shown below
# More details in https://huggingface.co/datasets/OpenDataDetector/ColliderML-Release-1,
# https://huggingface.co/docs/datasets/v1.10.1/package_reference/main_classes.html#datasets.IterableDataset
# and https://huggingface.co/docs/datasets/stream

x_key, y_key, z_key = ("x", "y", "z")
if args.forAutograph:
    x_key, y_key, z_key = ("true_x", "true_y", "true_z")

tracker_hits_iter_dataset = load_dataset(
    "OpenDataDetector/ColliderML-Release-1",
    data_dir="data/ttbar_pu200_tracker_hits",
    split="train",
    streaming=True,
    columns=[x_key,y_key,z_key,"particle_id"],
)

numberOfEvents = 0

for tracker_hits in tracker_hits_iter_dataset:

    # Need to set this for autograph; if needed, change the number of events to have more data
    # when using autograph
    if args.forAutograph:
        args.skipHeaders = True
        args.skipModules = True
        args.numberOfEvents = 1

    hx = np.asarray(tracker_hits[x_key], dtype=np.float32)
    hy = np.asarray(tracker_hits[y_key], dtype=np.float32)
    hz = np.asarray(tracker_hits[z_key], dtype=np.float32)
    hn = np.asarray((len(hx)), dtype=np.uint32)
    hr = np.sqrt(hx**2 + hy**2)
    hp = np.arctan2(hy,hx)
    hpart = np.asarray(tracker_hits["particle_id"], dtype=np.uint64)

    inputData = [hn,hx,hy,hz,hr,hp,hpart]

    if args.makeHitPosPlot and numberOfEvents == 0:

        plt.scatter(hz, hr, s=0.5)
        plt.title("ColliderML data")
        plt.xlabel("Z [mm]")
        plt.ylabel("R [mm]")
        plt.xlim(-1700,1700)
        plt.ylim(0, 200)
        if args.shortStrips:
            plt.xlim(-3100,3100)
            plt.ylim(0, 750)
        if args.longStrips:
            plt.xlim(-3100,3100)
            plt.ylim(0, 1200)
        plt.savefig("colliderMLGeometry.pdf")

    if args.debug:
        if numberOfEvents >= args.eventsForDebug: 
            print(f"Created files {hitsFileName} and {mapFileName} WHILE DEBUGGING")
            exit()

    saveX = []
    saveY = []
    saveZ = []
    saveR = []
    saveP = []
    saveLayersID = []
    saveParticleID = []

    saveModules = [0]

    extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(inputData,colliderMLPixelBarrel,colliderMLPixelBarrelThreshold,colliderMLPixelEndcap,colliderMLPixelEndcapThreshold,saveModules)

    saveX.extend(extractedX)
    saveY.extend(extractedY)
    saveZ.extend(extractedZ)
    saveR.extend(extractedR)
    saveP.extend(extractedP)
    saveLayersID.extend(extractedLayersID)
    saveParticleID.extend(extractedParticleID)

    if args.shortStrips:

        extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(inputData,coliderMLShortStripsBarrel,coliderMLShortStripsBarrelThreshold,coliderMLShortStripsEndcap,coliderMLShortStripsEndcapThreshold,saveModules)

        saveX.extend(extractedX)
        saveY.extend(extractedY)
        saveZ.extend(extractedZ)
        saveR.extend(extractedR)
        saveP.extend(extractedP)
        saveLayersID.extend(extractedLayersID)
        saveParticleID.extend(extractedParticleID)

    if args.longStrips:

        extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(inputData,coliderMLLongStripsBarrel,coliderMLLongStripsBarrelThreshold,coliderMLLongStripsEndcap,coliderMLLongStripsEndcapThreshold,saveModules)

        saveX.extend(extractedX)
        saveY.extend(extractedY)
        saveZ.extend(extractedZ)
        saveR.extend(extractedR)
        saveP.extend(extractedP)
        saveLayersID.extend(extractedLayersID)
        saveParticleID.extend(extractedParticleID)

    # Information per event is appended to file hits.txt as:
    # events:number_of_events (only added once)
    # hits:total_number_of_hits
    # modules:total_number_of_layers
    # ph,ph,dummy_error,dummy_error,xg,yg,zg,rg,phiShort,ph,ph,ph,layerID,partID
    # 0,nHits_layer(0),nHits_layer(1)+prev,...,nHitsLayer(N-1)+prev,nHitsLayer(N)+prev
    with open(hitsFileName, 'a') as fout:
        if not args.skipHeaders:
            if numberOfEvents == 0:
                if not args.debug: fout.write(f"events:{args.numberOfEvents}\n")
                else: fout.write(f"events:{args.eventsForDebug}\n")
        if not args.skipHeaders: fout.write(f"hits:{len(saveX)}\n")
        if not args.skipHeaders: fout.write(f"module:{len(saveModules)-1}\n")
        for j in range(len(saveX)):
            fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]:.5f},{saveY[j]:.5f},{saveZ[j]:.5f},{saveR[j]:.5f},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveModules[1]},{saveParticleID[j]}\n")
        if not args.skipModules:
            writeHelper = ""
            for m in saveModules: writeHelper = writeHelper + str(m) + ","
            writeHelper = writeHelper[:-1]
            fout.write(f"{writeHelper}\n")

    # Information per event is appended to file maps.txr as:
    # events:number_of_events (only added once)
    # hits:total_number_of_hits
    # partID (should match the same partID/order as in hits.txt)
    with open(mapFileName, 'a') as fout:
        if not args.skipHeaders:
            if numberOfEvents == 0:
                if not args.debug: fout.write(f"events:{args.numberOfEvents}\n")
                else: fout.write(f"events:{args.eventsForDebug}\n")
        if not args.skipHeaders: fout.write(f"hits:{len(saveParticleID)}\n")
        for j in range(len(saveParticleID)):
            fout.write(f"{saveParticleID[j]}\n")

    if args.debug: print(numberOfEvents)

    numberOfEvents += 1

    if numberOfEvents == args.numberOfEvents: break

print(f"Created files {hitsFileName} and {mapFileName} with {numberOfEvents} events")