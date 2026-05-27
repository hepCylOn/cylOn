import os

import numpy as np
import polars as pl
import math
import argparse

from colliderml.core import load_tables, collect_tables
from colliderml.polars import explode_tracker_hits

pl.Config.set_tbl_rows(12)
pl.Config.set_tbl_cols(20)

parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

parser.add_argument('-o', '--outFileName', default='')
parser.add_argument('-p', '--isPhase2', action='store_true')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-e', '--eventsForDebug', type=int, default=1)
parser.add_argument('-n', '--numberOfEvents', type=int, default=10)

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

saveX = []
saveY = []
saveZ = []
saveR = []
saveP = []
saveM = []
saveID = []
saveLayersID = []
saveParticleID = []

hits_flat = explode_tracker_hits(tracker_evt)

hn = (tracker_mult["n_tracker_hits"]).to_numpy().astype(int)
hx = hits_flat["x"].to_numpy().astype(float)
hy = hits_flat["y"].to_numpy().astype(float)
hz = hits_flat["z"].to_numpy().astype(float)
hr = np.sqrt(hx**2 + hy**2)
hp = np.arctan2(hy,hx)
hpart = hits_flat["particle_id"].to_numpy().astype(int)

counter = {"n": 0}
nAux = 0

for i in range(len(hn)):

    if args.debug:
        if counter["n"] >= args.eventsForDebug: 
            print(f"Created files {hitsFileName} and {mapFileName} WHILE DEBUGGING")
            exit()

    saveModules = [0]

    hxAux = hx[nAux:(nAux + hn[i])]
    hyAux = hy[nAux:(nAux + hn[i])]
    hzAux = hz[nAux:(nAux + hn[i])]
    hrAux = hr[nAux:(nAux + hn[i])]
    hpAux = hp[nAux:(nAux + hn[i])]
    hpartAux = hpart[nAux:(nAux + hn[i])]

    nAux += hn[i]

    # Only save hits that are in the pixel detector, i.e., |z| < 1600 and r < 200
    zCheck = (np.abs(hzAux) < ((colliderMLPixelEndcap[-1] + 80)*np.ones(len(hzAux))))

    hxAux = hxAux[zCheck]
    hyAux = hyAux[zCheck]
    hzAux = hzAux[zCheck]
    hrAux = hrAux[zCheck]
    hpAux = hpAux[zCheck]
    hpartAux = hpartAux[zCheck]

    xyCheck = (hrAux < (200.0*np.ones(len(hrAux))))
    
    hxAux = hxAux[xyCheck]
    hyAux = hyAux[xyCheck]
    hzAux = hzAux[xyCheck]
    hrAux = hrAux[xyCheck]
    hpAux = hpAux[xyCheck]
    hpartAux = hpartAux[xyCheck]

    # Initialize the layersID array for all the pixel hits
    saveLayersIDarray = np.zeros(len(hzAux))

    # This block saves layer ID in barrel, i.e., |z| < 550
    barrelCheck = (np.abs(hzAux) < (550*np.ones(len(hzAux))))

    # Do a loop over all of the radii for the barrel layers
    for rIdx in range(len(colliderMLPixelBarrel)):
        rHit = np.abs(hrAux)
        rBarrelLayer = colliderMLPixelBarrel[rIdx]*np.ones(len(hxAux))
        layerBarrelCheck = barrelCheck * ((rHit > (rBarrelLayer - colliderMLPixelBarrelThreshold)))
        layerBarrelCheck = layerBarrelCheck * ((rHit < (rBarrelLayer + colliderMLPixelBarrelThreshold)))
        saveLayersIDarray = saveLayersIDarray + (rIdx*layerBarrelCheck)

    # This block saves layer ID in the positive endcap, i.e., z > 550
    posEndcapCheck = (hzAux > (550*np.ones(len(hzAux))))

    # Do a loop over all of the longitudes for the endcap layers
    for zIdx in range(len(colliderMLPixelEndcap)):
        zHit = hzAux
        zEndcapLayer = colliderMLPixelEndcap[zIdx]*np.ones(len(hxAux))
        layerPosEndcapCheck = posEndcapCheck * ((zHit > (zEndcapLayer - colliderMLPixelEndcapThreshold)))
        layerPosEndcapCheck = layerPosEndcapCheck * ((zHit < (zEndcapLayer + colliderMLPixelEndcapThreshold)))
        saveLayersIDarray = saveLayersIDarray + ((zIdx + len(colliderMLPixelBarrel))*layerPosEndcapCheck)

    # This block saves layer ID in the negative endcap, i.e., z < -550
    negEndcapCheck = (hzAux < (-550*np.ones(len(hzAux))))

    # Do a loop over all of the longitudes for the endcap layers; the actual longitudinal
    # value has a signal change below to only consider z < 0
    for zIdx in range(len(colliderMLPixelEndcap)):
        zHit = hzAux
        zEndcapLayer = colliderMLPixelEndcap[zIdx]*np.ones(len(hxAux))
        layerNegEndcapCheck = negEndcapCheck *((zHit < -(zEndcapLayer - colliderMLPixelEndcapThreshold)))
        layerNegEndcapCheck = layerNegEndcapCheck *((zHit > -(zEndcapLayer + colliderMLPixelEndcapThreshold)))
        saveLayersIDarray = saveLayersIDarray + ((zIdx + len(colliderMLPixelBarrel) + len(colliderMLPixelEndcap))*layerNegEndcapCheck)

    # Check the amount of hits per layer of the pixel barrel, i.e., |z| < 550
    barrelCheck = (np.abs(hzAux) < (550*np.ones(len(hzAux))))

    barrelHxAux = hxAux[barrelCheck]
    barrelHrAux = hrAux[barrelCheck]

    for r in colliderMLPixelBarrel:
        rHit = np.abs(barrelHrAux)
        rBarrelLayer = r*np.ones(len(barrelHxAux))
        layerBarrelCheck = ((rHit > (rBarrelLayer - colliderMLPixelBarrelThreshold)))
        datasetBarrelLayers = rHit[layerBarrelCheck]
        rHit = np.abs(datasetBarrelLayers)
        rBarrelLayer = r*np.ones(len(datasetBarrelLayers))
        layerBarrelCheck = ((rHit < (rBarrelLayer + colliderMLPixelBarrelThreshold)))
        saveModules.append(len(datasetBarrelLayers[layerBarrelCheck]))

    # Check the amount of hits per layer of the pixel positive endcap, i.e., z > 550
    posEndcapCheck = ((hzAux) > (550*np.ones(len(hzAux))))
    
    posEndcapHxAux = hxAux[posEndcapCheck]
    posEndcapHzAux = hzAux[posEndcapCheck]

    for z in colliderMLPixelEndcap:
        zHit = posEndcapHzAux
        zEndcapLayer = z*np.ones(len(posEndcapHxAux))
        posEndcapLayerCheck = (zHit > (zEndcapLayer - colliderMLPixelEndcapThreshold))
        datasetPosEndcapLayers = posEndcapHzAux[posEndcapLayerCheck]
        zHit = datasetPosEndcapLayers
        zEndcapLayer = z*np.ones(len(datasetPosEndcapLayers))
        posEndcapLayerCheck = (zHit < (zEndcapLayer + colliderMLPixelEndcapThreshold))
        saveModules.append(len(datasetPosEndcapLayers[posEndcapLayerCheck]))

    # Check the amount of hits per layer of the pixel negative endcap, i.e., z < -550
    negEndcapCheck = ((hzAux) < (-550*np.ones(len(hzAux))))

    negEndcapHxAux = hxAux[negEndcapCheck]
    negEndcapHzAux = hzAux[negEndcapCheck]

    for z in colliderMLPixelEndcap:
        zHit = negEndcapHzAux
        zEndcapLayer = z*np.ones(len(negEndcapHxAux))
        negEndcapLayerCheck = (zHit < (-(zEndcapLayer - colliderMLPixelEndcapThreshold)))
        datasetNegEndcapLayers = zHit[negEndcapLayerCheck]
        zHit = datasetNegEndcapLayers
        zEndcapLayer = z*np.ones(len(datasetNegEndcapLayers))
        negEndcapLayerCheck = (zHit > (-(zEndcapLayer + colliderMLPixelEndcapThreshold)))
        saveModules.append(len(datasetNegEndcapLayers[negEndcapLayerCheck]))

    # Adds the values of hits per module cumulativelly to mimic the CMS hits input
    for i in range(1,len(saveModules)):
        saveModules[i] = saveModules[i] + saveModules[i-1]

    # To apply the layers ID ordering to all of the information to be saved, a tuple has to be
    # created so that the ordering acts in the same way over all columns. Afterwards, the tuple
    # can be split. Index 5 in x[5] means the layer ID column
    tuples = list(zip((hxAux/10.0).tolist(),(hyAux/10.0).tolist(),(hzAux/10.0).tolist(),(hrAux/10.0).tolist(),(np.arctan2(hyAux,hxAux)).tolist(),(saveLayersIDarray.astype('i')).tolist(),(hpartAux).tolist()))
    tuples_sorted = sorted(tuples, key = lambda x : x[5])

    saveX,saveY,saveZ,saveR,saveP,saveLayersID,saveParticleID = zip(*tuples_sorted)

    saveX = list(saveX)
    saveY = list(saveY)
    saveZ = list(saveZ)
    saveR = list(saveR)
    saveP = list(saveP)
    saveLayersID = list(saveLayersID)
    saveParticleID = list(saveParticleID)

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
            # fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveModules[1]}\n")
            fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveModules[1]},{saveParticleID[j]}\n") # Used to extract particles with more layers and for muons
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