import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import math
import vector

counter = {"n": 0}
debug = False

def showContent(name):
    print("=============================================")
    print("Started function!!")
    time1 = time.time()
    if debug:
        if counter["n"] >= 2: return True
    global f
    time2 = time.time()
    obj = f[name]
    if isinstance(obj, h5py.Dataset):
        time3 = time.time()
        savePXstr = ''
        savePYstr = ''
        savePZstr = ''
        saveVXstr = ''
        saveVYstr = ''
        saveVZstr = ''
        saveEstr = ''
        saveMstr = ''
        saveQstr = ''
        savePDGIDstr = ''
        savePartIDstr = ''
        if obj.dtype.names:
            for field in obj.dtype.names:
                if field.lower() == "px":
                    savePXstr = field
                if field.lower() == "py":
                    savePYstr = field
                if field.lower() == "pz":
                    savePZstr = field
                if field.lower() == "vx":
                    saveVXstr = field
                if field.lower() == "vy":
                    saveVYstr = field
                if field.lower() == "vz":
                    saveVZstr = field
                if field.lower() == "energy":
                    saveEstr = field
                if field.lower() == "mass":
                    saveMstr = field
                if field.lower() == "charge":
                    saveQstr = field
                if field.lower() == "particle_id":
                    savePartIDstr = field
                if field.lower() == "pdg_id":
                    savePDGIDstr = field
        dsetXY = obj[()]
        time4 = time.time()

        mXY = ~(abs(dsetXY[saveVZstr]) > (300.0*np.ones(len(dsetXY[saveVZstr]))))
        dsetXY = dsetXY[mXY]

        mXY = ~(np.sqrt((dsetXY[saveVXstr]*dsetXY[saveVXstr])+(dsetXY[saveVYstr]*dsetXY[saveVYstr])) > (25.0*np.ones(len(dsetXY[saveVXstr]))))
        dsetXY = dsetXY[mXY]

        p4_beforeCut = vector.arr({"px": dsetXY[savePXstr], "py": dsetXY[savePYstr], "pz": dsetXY[savePZstr], "E": dsetXY[saveEstr]})

        mXY = ~(p4_beforeCut.pt < (0.9*np.ones(len(dsetXY[savePXstr]))))
        dsetXY = dsetXY[mXY]

        p4_afterCut = vector.arr({"px": dsetXY[savePXstr], "py": dsetXY[savePYstr], "pz": dsetXY[savePZstr], "E": dsetXY[saveEstr]})

        time5 = time.time()

        tuples = list(zip((dsetXY[savePXstr]).tolist(),(dsetXY[savePYstr]).tolist(),(dsetXY[savePZstr]).tolist(),(dsetXY[saveEstr]).tolist(),(p4_afterCut.pt).tolist(),(p4_afterCut.eta).tolist(),(p4_afterCut.phi).tolist(),(dsetXY[saveMstr]).tolist(),(dsetXY[saveQstr]).tolist(),(dsetXY[savePDGIDstr]).tolist(),(dsetXY[savePartIDstr]).tolist(),(dsetXY[saveVXstr]).tolist(),(dsetXY[saveVYstr]).tolist(),(dsetXY[saveVZstr]).tolist()))
        tuples_sorted = sorted(tuples, key = lambda x : x[4])
        savePX,savePY,savePZ,saveE,savePT,saveEta,savePhi,saveM,saveQ,savePDGID,savePartID,saveVX,saveVY,saveVZ = zip(*tuples_sorted)

        time6 = time.time()

        savePX = list(savePX)
        savePY = list(savePY)
        savePZ = list(savePZ)
        saveE = list(saveE)
        savePT = list(savePT)
        saveEta = list(saveEta)
        savePhi = list(savePhi)
        saveM = list(saveM)
        saveQ = list(saveQ)
        savePDGID = list(savePDGID)
        savePartID = list(savePartID)
        saveVX = list(saveVX)
        saveVY = list(saveVY)
        saveVZ = list(saveVZ)

        time7 = time.time()

        with open('particles.txt', 'a') as fout:
            fout.write(f"particles:{len(savePX)}\n")
            for j in range(len(savePX)):
                fout.write(f"{saveVX[j]},{saveVY[j]},{saveVZ[j]},{savePX[j]},{savePY[j]},{savePZ[j]},{saveE[j]},{savePT[j]},{saveEta[j]},{savePhi[j]},{saveM[j]},{saveQ[j]},{savePDGID[j]},{savePartID[j]}\n")
            fout.close()

        time8 = time.time()

        print(f"Checking times: to get object names {time4 - time3} -- to perform cuts {time5 - time4} -- to sort by pt {time6 - time5} -- to convert to lists {time7 - time6} -- to write to file {time8 - time7}")

        counter["n"] += 1
        print(counter["n"])
    
    return None

f = h5py.File('/data/user/borzari/cmssw/pixeltrack-standalone/data/particles.h5', 'r')

fKeys = list(f.keys())

print(f.visit(showContent))