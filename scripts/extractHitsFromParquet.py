import os

import numpy as np
import polars as pl
import math
import argparse

from colliderml.core import load_tables, collect_tables
from colliderml.polars import explode_tracker_hits

import matplotlib.pyplot as plt

pl.Config.set_tbl_rows(12)
pl.Config.set_tbl_cols(20)

def writeToFile(i,n,inputData,barrelPos,barrelThres,endcapPos,endcapThres,modules):

    getInitialModulesLen = len(modules)

    hxAux = inputData[1][n:(n + inputData[0][i])]
    hyAux = inputData[2][n:(n + inputData[0][i])]
    hzAux = inputData[3][n:(n + inputData[0][i])]
    hrAux = inputData[4][n:(n + inputData[0][i])]
    hpAux = inputData[5][n:(n + inputData[0][i])]
    hpartAux = inputData[6][n:(n + inputData[0][i])]

    # Only save hits that are in the current part of detector being checked
    zCheck = (np.abs(hzAux) < ((endcapPos[-1] + endcapThres)*np.ones(len(hzAux))))

    hxAux = hxAux[zCheck]
    hyAux = hyAux[zCheck]
    hzAux = hzAux[zCheck]
    hrAux = hrAux[zCheck]
    hpAux = hpAux[zCheck]
    hpartAux = hpartAux[zCheck]

    xyCheck = (hrAux < ((barrelPos[-1] + barrelThres)*np.ones(len(hrAux))))
    
    hxAux = hxAux[xyCheck]
    hyAux = hyAux[xyCheck]
    hzAux = hzAux[xyCheck]
    hrAux = hrAux[xyCheck]
    hpAux = hpAux[xyCheck]
    hpartAux = hpartAux[xyCheck]

    xyCheck = (hrAux > ((barrelPos[0] - barrelThres)*np.ones(len(hrAux))))
    
    hxAux = hxAux[xyCheck]
    hyAux = hyAux[xyCheck]
    hzAux = hzAux[xyCheck]
    hrAux = hrAux[xyCheck]
    hpAux = hpAux[xyCheck]
    hpartAux = hpartAux[xyCheck]

    # Initialize the layersID array for all the pixel hits
    saveLayersIDarray = getInitialModulesLen*np.ones(len(hzAux)) - np.ones(len(hzAux))

    # This block saves layer ID in barrel, i.e., |z| < closestEndcapDisk
    barrelCheck = (np.abs(hzAux) < ((endcapPos[0] - endcapThres)*np.ones(len(hzAux))))

    # Do a loop over all of the radii for the barrel layers
    for rIdx in range(len(barrelPos)):
        rHit = np.abs(hrAux)
        rBarrelLayer = barrelPos[rIdx]*np.ones(len(hxAux))
        layerBarrelCheck = barrelCheck * ((rHit > (rBarrelLayer - barrelThres)))
        layerBarrelCheck = layerBarrelCheck * ((rHit < (rBarrelLayer + barrelThres)))
        saveLayersIDarray = saveLayersIDarray + (rIdx*layerBarrelCheck)

    # This block saves layer ID in the positive endcap, i.e., z > closestEndcapDisk
    posEndcapCheck = (hzAux > ((endcapPos[0] - endcapThres)*np.ones(len(hzAux))))

    # Do a loop over all of the longitudes for the endcap layers
    for zIdx in range(len(endcapPos)):
        zHit = hzAux
        zEndcapLayer = endcapPos[zIdx]*np.ones(len(hxAux))
        layerPosEndcapCheck = posEndcapCheck * ((zHit > (zEndcapLayer - endcapThres)))
        layerPosEndcapCheck = layerPosEndcapCheck * ((zHit < (zEndcapLayer + endcapThres)))
        saveLayersIDarray = saveLayersIDarray + ((zIdx + len(barrelPos))*layerPosEndcapCheck)

    # This block saves layer ID in the negative endcap, i.e., z < -closestEndcapDisk
    negEndcapCheck = (hzAux < (-(endcapPos[0] - endcapThres)*np.ones(len(hzAux))))

    # Do a loop over all of the longitudes for the endcap layers; the actual longitudinal
    # value has a signal change below to only consider z < 0
    for zIdx in range(len(endcapPos)):
        zHit = hzAux
        zEndcapLayer = endcapPos[zIdx]*np.ones(len(hxAux))
        layerNegEndcapCheck = negEndcapCheck *((zHit < -(zEndcapLayer - endcapThres)))
        layerNegEndcapCheck = layerNegEndcapCheck *((zHit > -(zEndcapLayer + endcapThres)))
        saveLayersIDarray = saveLayersIDarray + ((zIdx + len(barrelPos) + len(endcapPos))*layerNegEndcapCheck)

    # Check the amount of hits per layer of the pixel barrel, i.e., |z| < closestEndcapDisk
    barrelCheck = (np.abs(hzAux) < ((endcapPos[0] - endcapThres)*np.ones(len(hzAux))))

    barrelHxAux = hxAux[barrelCheck]
    barrelHrAux = hrAux[barrelCheck]

    for r in barrelPos:
        rHit = np.abs(barrelHrAux)
        rBarrelLayer = r*np.ones(len(barrelHxAux))
        layerBarrelCheck = ((rHit > (rBarrelLayer - barrelThres)))
        datasetBarrelLayers = rHit[layerBarrelCheck]
        rBarrelLayer = r*np.ones(len(datasetBarrelLayers))
        layerBarrelCheck = ((datasetBarrelLayers < (rBarrelLayer + barrelThres)))
        # if detector == 2: print(len(datasetBarrelLayers[layerBarrelCheck]))
        modules.append(len(datasetBarrelLayers[layerBarrelCheck]))

    # Check the amount of hits per layer of the pixel positive endcap, i.e., z > closestEndcapDisk
    posEndcapCheck = ((hzAux) > ((endcapPos[0] - endcapThres)*np.ones(len(hzAux))))
    
    posEndcapHxAux = hxAux[posEndcapCheck]
    posEndcapHzAux = hzAux[posEndcapCheck]

    for z in endcapPos:
        zHit = posEndcapHzAux
        zEndcapLayer = z*np.ones(len(posEndcapHxAux))
        posEndcapLayerCheck = (zHit > (zEndcapLayer - endcapThres))
        datasetPosEndcapLayers = posEndcapHzAux[posEndcapLayerCheck]
        zEndcapLayer = z*np.ones(len(datasetPosEndcapLayers))
        posEndcapLayerCheck = (datasetPosEndcapLayers < (zEndcapLayer + endcapThres))
        # if detector == 2: print(len(datasetPosEndcapLayers[posEndcapLayerCheck]))
        modules.append(len(datasetPosEndcapLayers[posEndcapLayerCheck]))

    # Check the amount of hits per layer of the pixel negative endcap, i.e., z < -closestEndcapDisk
    negEndcapCheck = ((hzAux) < (-(endcapPos[0] - endcapThres)*np.ones(len(hzAux))))

    negEndcapHxAux = hxAux[negEndcapCheck]
    negEndcapHzAux = hzAux[negEndcapCheck]

    for z in endcapPos:
        zHit = negEndcapHzAux
        zEndcapLayer = z*np.ones(len(negEndcapHxAux))
        negEndcapLayerCheck = (zHit < (-(zEndcapLayer - endcapThres)))
        datasetNegEndcapLayers = zHit[negEndcapLayerCheck]
        zEndcapLayer = z*np.ones(len(datasetNegEndcapLayers))
        negEndcapLayerCheck = (datasetNegEndcapLayers > (-(zEndcapLayer + endcapThres)))
        # if detector == 2: print(len(datasetNegEndcapLayers[negEndcapLayerCheck]))
        modules.append(len(datasetNegEndcapLayers[negEndcapLayerCheck]))

    # Adds the values of hits per module cumulativelly to mimic the CMS hits input
    for i in range(getInitialModulesLen,len(modules)):
        modules[i] = modules[i] + modules[i-1]

    # To apply the layers ID ordering to all of the information to be saved, a tuple has to be
    # created so that the ordering acts in the same way over all columns. Afterwards, the tuple
    # can be split. Index 5 in x[5] means the layer ID column
    tuples = list(zip((hxAux/10.0).tolist(),(hyAux/10.0).tolist(),(hzAux/10.0).tolist(),(hrAux/10.0).tolist(),(np.arctan2(hyAux,hxAux)).tolist(),(saveLayersIDarray.astype('i')).tolist(),(hpartAux).tolist()))
    tuples_sorted = sorted(tuples, key = lambda x : x[5])

    saveX,saveY,saveZ,saveR,saveP,saveLayersID,saveParticleID = zip(*tuples_sorted)

    return saveX,saveY,saveZ,saveR,saveP,saveLayersID,saveParticleID,modules

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

hitsFileName = 'hitsWithParticleID' + args.outFileName + phaseModifier + '.txt'
mapFileName = 'mapHitsToParticles' + args.outFileName + phaseModifier + '.txt'

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

home_directory = os.path.expanduser("~")

cfg = {
    "dataset_id": "CERN/ColliderML-Release-1",
    "channels": "ttbar",
    "pileup": "pu200",
    "objects": ["tracker_hits"],
    "split": "train",
    "lazy": True,
    "max_events": args.numberOfEvents,
    "data_dir": home_directory + "/.cache/colliderml"
}

tables = load_tables(cfg)
frames = collect_tables(tables)

tracker_evt = frames["tracker_hits"]

tracker_mult = tracker_evt.select(
    pl.col("x").list.len().alias("n_tracker_hits"),
)

hits_flat = explode_tracker_hits(tracker_evt)

hn = (tracker_mult["n_tracker_hits"]).to_numpy().astype(int)
hx = hits_flat["x"].to_numpy().astype(float)
hy = hits_flat["y"].to_numpy().astype(float)
hz = hits_flat["z"].to_numpy().astype(float)
hr = np.sqrt(hx**2 + hy**2)
hp = np.arctan2(hy,hx)
hpart = hits_flat["particle_id"].to_numpy().astype(int)

inputData = [hn,hx,hy,hz,hr,hp,hpart]

counter = {"n": 0}
nAux = 0

if args.makeHitPosPlot:

    plt.scatter(hz, hr, s=0.5)
    plt.title("ColliderML data")
    plt.xlabel("Z [mm]")
    plt.ylabel("R [mm]")
    plt.savefig("colliderMLGeometry.pdf")

for i in range(len(hn)):

    if args.debug:
        if counter["n"] >= args.eventsForDebug: 
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

    extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(i,nAux,inputData,colliderMLPixelBarrel,colliderMLPixelBarrelThreshold,colliderMLPixelEndcap,colliderMLPixelEndcapThreshold,saveModules)

    saveX.extend(list(extractedX))
    saveY.extend(list(extractedY))
    saveZ.extend(list(extractedZ))
    saveR.extend(list(extractedR))
    saveP.extend(list(extractedP))
    saveLayersID.extend(list(extractedLayersID))
    saveParticleID.extend(list(extractedParticleID))

    if args.shortStrips:

        extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(i,nAux,inputData,coliderMLShortStripsBarrel,coliderMLShortStripsBarrelThreshold,coliderMLShortStripsEndcap,coliderMLShortStripsEndcapThreshold,saveModules)

        saveX.extend(list(extractedX))
        saveY.extend(list(extractedY))
        saveZ.extend(list(extractedZ))
        saveR.extend(list(extractedR))
        saveP.extend(list(extractedP))
        saveLayersID.extend(list(extractedLayersID))
        saveParticleID.extend(list(extractedParticleID))

    if args.longStrips:

        extractedX,extractedY,extractedZ,extractedR,extractedP,extractedLayersID,extractedParticleID,saveModules = writeToFile(i,nAux,inputData,coliderMLLongStripsBarrel,coliderMLLongStripsBarrelThreshold,coliderMLLongStripsEndcap,coliderMLLongStripsEndcapThreshold,saveModules)

        saveX.extend(list(extractedX))
        saveY.extend(list(extractedY))
        saveZ.extend(list(extractedZ))
        saveR.extend(list(extractedR))
        saveP.extend(list(extractedP))
        saveLayersID.extend(list(extractedLayersID))
        saveParticleID.extend(list(extractedParticleID))

    nAux += inputData[0][i]

    # Information per event is appended to file hits.txt as:
    # hits:total_number_of_hits
    # modules:total_number_of_layers
    # ph,ph,dummy_error,dummy_error,xg,yg,zg,rg,phiShort,ph,ph,ph,layerID,partID
    # 0,nHits_layer(0),nHits_layer(1)+prev,...,nHitsLayer(N-1)+prev,nHitsLayer(N)+prev
    with open(hitsFileName, 'a') as fout:
        if counter["n"] == 0: fout.write(f"events:{len(hn)}\n")
        fout.write(f"hits:{len(saveX)}\n")
        fout.write(f"module:{len(saveModules)-1}\n")
        for j in range(len(saveX)):
            fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveModules[1]},{saveParticleID[j]}\n")
        writeHelper = ""
        for m in saveModules: writeHelper = writeHelper + str(m) + ","
        # 
        writeHelper = writeHelper[:-1]
        fout.write(f"{writeHelper}\n")

    with open(mapFileName, 'a') as fout:
        if counter["n"] == 0: fout.write(f"events:{len(hn)}\n")
        fout.write(f"hits:{len(saveParticleID)}\n")
        for j in range(len(saveParticleID)):
            fout.write(f"{saveParticleID[j]}\n")

    if args.debug: print(counter["n"])

    counter["n"] += 1

print(f"Created files {hitsFileName} and {mapFileName} with {args.numberOfEvents} events")