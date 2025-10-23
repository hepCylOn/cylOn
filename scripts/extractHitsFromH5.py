import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import math

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

colliderMLPixelBarrel = [34.0,70.0,116.0,172.0]
# colliderMLPixelEndcap = [650.0,780.0,900.0,1000.0,1200.0,1400.0,1600.0]
colliderMLPixelEndcap = [620.0,730.0,830.0]
colliderMLPixelBarrelThreshold = 15.0
colliderMLPixelEndcapThreshold = 50.0

counter = {"n": 0}
debug = False

def showContent(name):
    print("=============================================")
    print("Started function!!")
    if debug:
        if counter["n"] >= 2: return True
    global f
    obj = f[name]
    saveModules = [0]
    if isinstance(obj, h5py.Dataset):
        time3 = time.time()
        saveXstr = ''
        saveYstr = ''
        saveZstr = ''
        saveIDstr = ''
        saveParticleIDstr = ''
        # print(name)
        if obj.dtype.names:
            for field in obj.dtype.names:
                print(field)
                if field.lower() == "x":
                    saveXstr = field
                if field.lower() == "y":
                    saveYstr = field
                if field.lower() == "z":
                    saveZstr = field
                if field.lower() == "surface_id":
                    saveIDstr = field
                if field.lower() == "particle_id":
                    saveParticleIDstr = field
        time4 = time.time()
        mXY = ~((np.abs(obj[()][saveZstr]) > (950*np.ones(len(obj[()][saveZstr])))))
        dsetXY = obj[()][mXY]
        mXY = ~((np.abs(np.sqrt((dsetXY[saveXstr]*dsetXY[saveXstr]) + (dsetXY[saveYstr]*dsetXY[saveYstr]))) > (200.0*np.ones(len(dsetXY[saveXstr])))))
        dsetXY = dsetXY[mXY]
        # dsetXY = obj[()]
        saveLayersIDarray = np.zeros(len(dsetXY[saveZstr]))

        # Save layer ID in barrel
        mIDBarrel = ~((np.abs(dsetXY[saveZstr]) > (550*np.ones(len(dsetXY[saveZstr])))))
        for ri in range(len(colliderMLPixelBarrel)):
            rXYBarrel = np.abs(np.sqrt((dsetXY[saveXstr]*dsetXY[saveXstr]) + (dsetXY[saveYstr]*dsetXY[saveYstr])))
            rBarrelLayer = colliderMLPixelBarrel[ri]*np.ones(len(dsetXY[saveXstr]))
            mIDLayerBarrel = mIDBarrel *(~((rXYBarrel < (rBarrelLayer - colliderMLPixelBarrelThreshold))))
            mIDLayerBarrel = mIDLayerBarrel *(~((rXYBarrel > (rBarrelLayer + colliderMLPixelBarrelThreshold))))
            saveLayersIDarray = saveLayersIDarray + (ri*mIDLayerBarrel)

        # Save layer ID in positive endcap
        mIDEndcap = ~((np.abs(dsetXY[saveZstr]) < (550*np.ones(len(dsetXY[saveZstr])))))
        for zi in range(len(colliderMLPixelEndcap)):
            zZEndcap = dsetXY[saveZstr]
            zEndcapLayer = colliderMLPixelEndcap[zi]*np.ones(len(dsetXY[saveXstr]))
            mIDLayerEndcap = mIDEndcap *(~((zZEndcap < (zEndcapLayer - colliderMLPixelEndcapThreshold))))
            mIDLayerEndcap = mIDLayerEndcap *(~((zZEndcap > (zEndcapLayer + colliderMLPixelEndcapThreshold))))
            saveLayersIDarray = saveLayersIDarray + ((zi + 4)*mIDLayerEndcap)

        # Save layer ID in negative endcap
        mIDEndcap = ~((np.abs(dsetXY[saveZstr]) < (550*np.ones(len(dsetXY[saveZstr])))))
        for zi in range(len(colliderMLPixelEndcap)):
            zZEndcap = dsetXY[saveZstr]
            zEndcapLayer = colliderMLPixelEndcap[zi]*np.ones(len(dsetXY[saveXstr]))
            mIDLayerEndcap = mIDEndcap *(~((zZEndcap > -(zEndcapLayer - colliderMLPixelEndcapThreshold))))
            mIDLayerEndcap = mIDLayerEndcap *(~((zZEndcap < -(zEndcapLayer + colliderMLPixelEndcapThreshold))))
            saveLayersIDarray = saveLayersIDarray + ((zi + 7)*mIDLayerEndcap)

        # Check the amount of hits per layer of the pixel barrel
        mXY = ~((np.abs(dsetXY[saveZstr]) > (550*np.ones(len(dsetXY[saveZstr])))))
        dCheckLayers = dsetXY[mXY]
        for r in colliderMLPixelBarrel:
            rXYBarrel = np.abs(np.sqrt((dCheckLayers[saveXstr]*dCheckLayers[saveXstr]) + (dCheckLayers[saveYstr]*dCheckLayers[saveYstr])))
            rBarrelLayer = r*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~((rXYBarrel < (rBarrelLayer - colliderMLPixelBarrelThreshold)))
            dCheckLayers = dCheckLayers[mXY]

            rXYBarrel = np.abs(np.sqrt((dCheckLayers[saveXstr]*dCheckLayers[saveXstr]) + (dCheckLayers[saveYstr]*dCheckLayers[saveYstr])))
            rBarrelLayer = r*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~((rXYBarrel > (rBarrelLayer + colliderMLPixelBarrelThreshold)))
            saveModules.append(len(dCheckLayers[mXY]))

        # Check the amount of hits per layer of the pixel positive endcap
        mXY = ~(((dsetXY[saveZstr]) < (550*np.ones(len(dsetXY[saveZstr])))))
        dCheckLayers = dsetXY[mXY]
        for z in colliderMLPixelEndcap:
            zZEndcap = dCheckLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~(zZEndcap < (zEndcapLayer - colliderMLPixelEndcapThreshold))
            dCheckLayers = dCheckLayers[mXY]

            zZEndcap = dCheckLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~(zZEndcap > (zEndcapLayer + colliderMLPixelEndcapThreshold))
            saveModules.append(len(dCheckLayers[mXY]))

        # Check the amount of hits per layer of the pixel negative endcap
        mXY = ~(((dsetXY[saveZstr]) > (-550*np.ones(len(dsetXY[saveZstr])))))
        dCheckLayers = dsetXY[mXY]
        for z in colliderMLPixelEndcap:
            zZEndcap = dCheckLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~(zZEndcap > (-(zEndcapLayer - colliderMLPixelEndcapThreshold)))
            dCheckLayers = dCheckLayers[mXY]

            zZEndcap = dCheckLayers[saveZstr]
            zEndcapLayer = z*np.ones(len(dCheckLayers[saveXstr]))

            mXY = ~(zZEndcap < (-(zEndcapLayer + colliderMLPixelEndcapThreshold)))
            saveModules.append(len(dCheckLayers[mXY]))

        for i in range(1,len(saveModules)):
            saveModules[i] = saveModules[i] + saveModules[i-1]

        time5 = time.time()

        tuples = list(zip((dsetXY[saveXstr]).tolist(),(dsetXY[saveYstr]).tolist(),(dsetXY[saveZstr]).tolist(),(np.sqrt((dsetXY[saveXstr]*dsetXY[saveXstr]) + (dsetXY[saveYstr]*dsetXY[saveYstr]))).tolist(),(np.arctan2(dsetXY[saveYstr],dsetXY[saveXstr])).tolist(),(saveLayersIDarray.astype('i')).tolist(),(dsetXY[saveParticleIDstr]).tolist()))
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

        with open('hits.txt', 'a') as fout:
            fout.write(f"hits:{len(saveX)}\n")
            for j in range(len(saveX)):
                fout.write(f"1.0,1.0,0.01,0.01,{saveX[j]},{saveY[j]},{saveZ[j]},{saveR[j]},{phi2short(saveP[j])},2,3,4,{saveLayersID[j]},{saveParticleID[j]}\n")
            fout.write(f"module:{len(saveModules)-1}\n")
            writeHelper = ""
            for m in saveModules: writeHelper = writeHelper + str(m) + ","
            writeHelper = writeHelper[:-1]
            fout.write(f"{writeHelper}\n")

        time8 = time.time()

        print(f"Checking times: to get object names {time4 - time3} -- to get layers {time5 - time4} -- to sort by layer ID {time6 - time5} -- to convert to lists {time7 - time6} -- to write to file {time8 - time7}")

        counter["n"] += 1
        print(counter["n"])
    
    return None

f = h5py.File('/data/user/borzari/cmssw/pixeltrack-standalone/data/hits.h5', 'r')

fKeys = list(f.keys())

print(f.visit(showContent))