import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import math
import os
import pickle
import struct

def phi2short(x: float) -> int:
    p2i = (2**15) / math.pi  # 32768 / pi
    val = round(x * p2i)

    # simular wraparound de int16
    if val > 32767:
        val -= 65536
    elif val < -32768:
        val += 65536

    return val

saveX = []
saveY = []
saveZ = []
saveR = []
saveP = []
saveM = []
saveID = []
saveLayersID = []
saveParticleID = []

# These define the position of the layers in the open data detector geometry
# taken from https://iopscience.iop.org/article/10.1088/1742-6596/2438/1/012110/pdf
# Thresholds are set to get hits in all of a given layer, but not in other layers
colliderMLPixelBarrel = [34.0,70.0,116.0,172.0]
# colliderMLPixelEndcap = [620.0,730.0,830.0,980.0,1120.0,1320.0,1520.0]
colliderMLPixelEndcap = [620.0,730.0,830.0]
colliderMLPixelBarrelThreshold = 15.0
colliderMLPixelEndcapThreshold = 50.0

# Added counter to run over a given amount of events for debugging
counter = {"n": 0}
debug = False
# debug = True

# Since the hits information is appended to the hits.txt file, it is better
# to remove it before executing the code
if os.path.exists('hitsWithoutParticleId.txt'): os.system('rm hitsWithoutParticleId.txt')

# Function that accesses the information inside of the .h5 file
def showContent(name):
    if debug:
        if counter["n"] >= 1: return True
    global f
    obj = f[name]
    saveModules = [0]
    if isinstance(obj, h5py.Dataset):
        time3 = time.time()
        # Gets names of variables if they are lower or upper case; might be useful
        # when using this code for other files or other file types
        saveXstr = ''
        saveYstr = ''
        saveZstr = ''
        saveParticleIDstr = ''
        if obj.dtype.names:
            for field in obj.dtype.names:
                if field.lower() == "x":
                    saveXstr = field
                if field.lower() == "y":
                    saveYstr = field
                if field.lower() == "z":
                    saveZstr = field
                if field.lower() == "particle_id":
                    saveParticleIDstr = field

        time4 = time.time()

        # Only save hits that are in the pixel detector, i.e., |z| < 1600 and r < 200
        zCheck = (np.abs(obj[()][saveZstr]) < ((colliderMLPixelEndcap[-1] + 80)*np.ones(len(obj[()][saveZstr]))))
        dataset = obj[()][zCheck]
        xyCheck = ((np.sqrt((dataset[saveXstr]*dataset[saveXstr]) + (dataset[saveYstr]*dataset[saveYstr]))) < (200.0*np.ones(len(dataset[saveXstr]))))
        dataset = dataset[xyCheck]

        # Initialize the layersID array for all the pixel hits
        saveLayersIDarray = np.zeros(len(dataset[saveZstr]))

        # This block saves layer ID in barrel, i.e., |z| < 550
        barrelCheck = (np.abs(dataset[saveZstr]) < (550*np.ones(len(dataset[saveZstr]))))
        # Do a loop over all of the radii for the barrel layers
        for rIdx in range(len(colliderMLPixelBarrel)):
            rHit = np.abs(np.sqrt((dataset[saveXstr]*dataset[saveXstr]) + (dataset[saveYstr]*dataset[saveYstr])))
            rBarrelLayer = colliderMLPixelBarrel[rIdx]*np.ones(len(dataset[saveXstr]))
            layerBarrelCheck = barrelCheck * ((rHit > (rBarrelLayer - colliderMLPixelBarrelThreshold)))
            layerBarrelCheck = layerBarrelCheck * ((rHit < (rBarrelLayer + colliderMLPixelBarrelThreshold)))
            saveLayersIDarray = saveLayersIDarray + (rIdx*layerBarrelCheck)

        # This block saves layer ID in the positive endcap, i.e., z > 550
        posEndcapCheck = (dataset[saveZstr] > (550*np.ones(len(dataset[saveZstr]))))
        # Do a loop over all of the longitudes for the endcap layers
        for zIdx in range(len(colliderMLPixelEndcap)):
            zHit = dataset[saveZstr]
            zEndcapLayer = colliderMLPixelEndcap[zIdx]*np.ones(len(dataset[saveXstr]))
            layerPosEndcapCheck = posEndcapCheck * ((zHit > (zEndcapLayer - colliderMLPixelEndcapThreshold)))
            layerPosEndcapCheck = layerPosEndcapCheck * ((zHit < (zEndcapLayer + colliderMLPixelEndcapThreshold)))
            saveLayersIDarray = saveLayersIDarray + ((zIdx + len(colliderMLPixelBarrel))*layerPosEndcapCheck)

        # This block saves layer ID in the negative endcap, i.e., z < -550
        negEndcapCheck = (dataset[saveZstr] < (-550*np.ones(len(dataset[saveZstr]))))
        # Do a loop over all of the longitudes for the endcap layers; the actual longitudinal
        # value has a signal change below to only consider z < 0
        for zIdx in range(len(colliderMLPixelEndcap)):
            zHit = dataset[saveZstr]
            zEndcapLayer = colliderMLPixelEndcap[zIdx]*np.ones(len(dataset[saveXstr]))
            layerNegEndcapCheck = negEndcapCheck *((zHit < -(zEndcapLayer - colliderMLPixelEndcapThreshold)))
            layerNegEndcapCheck = layerNegEndcapCheck *((zHit > -(zEndcapLayer + colliderMLPixelEndcapThreshold)))
            saveLayersIDarray = saveLayersIDarray + ((zIdx + len(colliderMLPixelBarrel) + len(colliderMLPixelEndcap))*layerNegEndcapCheck)

        # Check the amount of hits per layer of the pixel barrel, i.e., |z| < 550
        barrelCheck = (np.abs(dataset[saveZstr]) < (550*np.ones(len(dataset[saveZstr]))))
        datasetBarrelLayers = dataset[barrelCheck]
        for r in colliderMLPixelBarrel:
            rHit = np.abs(np.sqrt((datasetBarrelLayers[saveXstr]*datasetBarrelLayers[saveXstr]) + (datasetBarrelLayers[saveYstr]*datasetBarrelLayers[saveYstr])))
            rBarrelLayer = r*np.ones(len(datasetBarrelLayers[saveXstr]))

            layerBarrelCheck = ((rHit > (rBarrelLayer - colliderMLPixelBarrelThreshold)))
            datasetBarrelLayers = datasetBarrelLayers[layerBarrelCheck]

            rHit = np.abs(np.sqrt((datasetBarrelLayers[saveXstr]*datasetBarrelLayers[saveXstr]) + (datasetBarrelLayers[saveYstr]*datasetBarrelLayers[saveYstr])))
            rBarrelLayer = r*np.ones(len(datasetBarrelLayers[saveXstr]))

            layerBarrelCheck = ((rHit < (rBarrelLayer + colliderMLPixelBarrelThreshold)))
            saveModules.append(len(datasetBarrelLayers[layerBarrelCheck]))

        # Check the amount of hits per layer of the pixel positive endcap, i.e., z > 550
        posEndcapCheck = ((dataset[saveZstr]) > (550*np.ones(len(dataset[saveZstr]))))
        datasetPosEndcapLayers = dataset[posEndcapCheck]
        for z in colliderMLPixelEndcap:
            zHit = datasetPosEndcapLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(datasetPosEndcapLayers[saveXstr]))

            posEndcapLayerCheck = (zHit > (zEndcapLayer - colliderMLPixelEndcapThreshold))
            datasetPosEndcapLayers = datasetPosEndcapLayers[posEndcapLayerCheck]

            zHit = datasetPosEndcapLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(datasetPosEndcapLayers[saveXstr]))

            posEndcapLayerCheck = (zHit < (zEndcapLayer + colliderMLPixelEndcapThreshold))
            saveModules.append(len(datasetPosEndcapLayers[posEndcapLayerCheck]))

        # Check the amount of hits per layer of the pixel negative endcap, i.e., z < -550
        negEndcapCheck = ((dataset[saveZstr]) < (-550*np.ones(len(dataset[saveZstr]))))
        datasetNegEndcapLayers = dataset[negEndcapCheck]
        for z in colliderMLPixelEndcap:
            zHit = datasetNegEndcapLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(datasetNegEndcapLayers[saveXstr]))

            negEndcapLayerCheck = (zHit < (-(zEndcapLayer - colliderMLPixelEndcapThreshold)))
            datasetNegEndcapLayers = datasetNegEndcapLayers[negEndcapLayerCheck]

            zHit = datasetNegEndcapLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(datasetNegEndcapLayers[saveXstr]))

            negEndcapLayerCheck = (zHit > (-(zEndcapLayer + colliderMLPixelEndcapThreshold)))
            saveModules.append(len(datasetNegEndcapLayers[negEndcapLayerCheck]))

        # Adds the values of hits per module cumulativelly to mimic the CMS hits input
        for i in range(1,len(saveModules)):
            saveModules[i] = saveModules[i] + saveModules[i-1]

        time5 = time.time()

        # To apply the layers ID ordering to all of the information to be saved, a tuple has to be
        # created so that the ordering acts in the same way over all columns. Afterwards, the tuple
        # can be split. Index 5 in x[5] means the layer ID column
        tuples = list(zip((dataset[saveXstr]/10.0).tolist(),(dataset[saveYstr]/10.0).tolist(),(dataset[saveZstr]/10.0).tolist(),(np.sqrt((dataset[saveXstr]*dataset[saveXstr]) + (dataset[saveYstr]*dataset[saveYstr]))/10.0).tolist(),(np.arctan2(dataset[saveYstr],dataset[saveXstr])).tolist(),(saveLayersIDarray.astype('i')).tolist(),(dataset[saveParticleIDstr]).tolist()))
        tuples_sorted = sorted(tuples, key = lambda x : x[5])
        saveX,saveY,saveZ,saveR,saveP,saveLayersID,saveParticleID = zip(*tuples_sorted)

        time6 = time.time()

        saveX = list(saveX)
        saveY = list(saveY)
        saveZ = list(saveZ)
        saveR = list(saveR)
        saveP = list(saveP)
        saveLayersID = list(saveLayersID)
        saveParticleID = list(saveParticleID)

        time7 = time.time()

        # Information per event is appended to file hits.txt as:
        # hits:total_number_of_hits
        # ph,ph,dummy_error,dummy_error,xg,yg,zg,rg,phiShort,ph,ph,ph,layerID,partID
        # modules:total_number_of_layers
        # 0,nHits_layer(0),nHits_layer(1)+prev,...,nHitsLayer(N-1)+prev,nHitsLayer(N)+prev
        with open('hitsWithoutParticleId.txt', 'a') as fout:
            if counter["n"] == 0: fout.write(f"events:1000\n")
            fout.write(f"hits:{len(saveX)}\n")
            fout.write(f"module:{len(saveModules)-1}\n")
            for j in range(len(saveX)):
                #fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveParticleID[j]}\n")
                fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveModules[1]}\n")
            writeHelper = ""
            for m in saveModules: writeHelper = writeHelper + str(m) + ","
            
            writeHelper = writeHelper[:-1]
            fout.write(f"{writeHelper}\n")

        time8 = time.time()
        # Checks execution time in distinct blocks:
        # - To get names (very fast -> 0.2 ms)
        # - To get data and layers (slow -> 160 ms)
        # - To sort by layer ID (fast -> 38 ms)
        # - To convert to lists (very fast -> 3 ms)
        # - To write to file (very slow -> 247 ms)
        # Still process about 2.5 ev/s which is manageable
        print("=====================================")
        print(f"Checking times: to get object names {time4 - time3} -- to get data and layers {time5 - time4} -- to sort by layer ID {time6 - time5} -- to convert to lists {time7 - time6} -- to write to file {time8 - time7}")

        counter["n"] += 1
        print(counter["n"])
    
    return None

f = h5py.File('/home/breno/data/full_pileup_pilot/ttbar/v2/reco/tracker_hits/events0-999.h5', 'r')

fKeys = list(f.keys())

print(f.visit(showContent))

