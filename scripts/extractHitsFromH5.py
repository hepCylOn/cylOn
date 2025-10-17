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
debug = True

def showContent(name):
    if debug:
        if counter["n"] >= 2: return True
    obj = f[name]
    saveModules = [0]
    if isinstance(obj, h5py.Dataset):
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

        saveX.extend([(dsetXY[saveXstr]).tolist()])
        saveY.extend([(dsetXY[saveYstr]).tolist()])
        saveZ.extend([(dsetXY[saveZstr]).tolist()])
        saveR.extend([(np.sqrt((dsetXY[saveXstr]*dsetXY[saveXstr]) + (dsetXY[saveYstr]*dsetXY[saveYstr]))).tolist()])
        saveP.extend([(np.arctan2(dsetXY[saveYstr],dsetXY[saveXstr])).tolist()])
        saveM.extend([saveModules])
        saveID.extend([(dsetXY[saveIDstr]).tolist()])
        saveParticleID.extend([(dsetXY[saveParticleIDstr]).tolist()])
        saveLayersIDarray = saveLayersIDarray.astype('i')
        # plt.hist((dsetXY[saveParticleIDstr]),bins=500)
        # plt.show()
        saveLayersID.extend([saveLayersIDarray.tolist()])
        counter["n"] += 1
    
    return None

f = h5py.File('data/full_pileup_pilot/ttbar/v2/reco/tracker_hits/events0-999.h5', 'r')

fKeys = list(f.keys())

print(f.visit(showContent))

plt.scatter(saveZ[0],saveR[0])
plt.show()

# Sorting hits by layerID
for i in range(len(saveX)):
    tuples = list(zip(saveX[i],saveY[i],saveZ[i],saveR[i],saveP[i],saveLayersID[i],saveParticleID[i]))
    tuples_sorted = sorted(tuples, key = lambda x : x[5])
    saveX[i],saveY[i],saveZ[i],saveR[i],saveP[i],saveLayersID[i],saveParticleID[i] = zip(*tuples_sorted)

begtime = time.time()

with open('hitsPixel.txt', 'w') as f:
    for i in range(len(saveX)):
        f.write(f"hits:{len(saveX[i])}\n")
        for j in range(len(saveX[i])):
            # f.write(f"1.0,1.0,0.01,0.01,{saveX[i][j]},{saveY[i][j]},{saveZ[i][j]},{saveR[i][j]},{phi2short(saveP[i][j])},2,3,4,{saveID[i][j]}\n")
            f.write(f"1.0,1.0,0.01,0.01,{saveX[i][j]},{saveY[i][j]},{saveZ[i][j]},{saveR[i][j]},{phi2short(saveP[i][j])},2,3,4,{saveLayersID[i][j]},{saveParticleID[i][j]}\n")
        f.write(f"module:{len(saveM[0])-1}\n")
        writeHelper = ""
        for m in saveM[i]: writeHelper = writeHelper + str(m) + ","
        writeHelper = writeHelper[:-1]
        f.write(f"{writeHelper}\n")

endtime = time.time()

tottime = endtime - begtime
print(tottime)
