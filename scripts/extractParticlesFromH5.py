import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import math
import vector

savePX = []
savePY = []
savePZ = []
saveE = []
savePT = []
saveEta = []
savePhi = []
saveM = []
saveVX = []
saveVY = []
saveVZ = []
saveQ = []
savePDGID = []
savePartID = []

counter = {"n": 0}
debug = True

def showContent(name):
    if debug:
        if counter["n"] >= 2: return True
    obj = f[name]
    if isinstance(obj, h5py.Dataset):
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
                # print(field)
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

        p4_beforeCut = vector.arr({"px": dsetXY[savePXstr], "py": dsetXY[savePYstr], "pz": dsetXY[savePZstr], "E": dsetXY[saveEstr]})

        mXY = ~(p4_beforeCut.pt < (1.0*np.ones(len(dsetXY[savePXstr]))))
        dsetXY = dsetXY[mXY]

        p4_afterCut = vector.arr({"px": dsetXY[savePXstr], "py": dsetXY[savePYstr], "pz": dsetXY[savePZstr], "E": dsetXY[saveEstr]})

        plt.hist((dsetXY[savePartIDstr]),bins=500)
        plt.show()

        savePX.extend([(dsetXY[savePXstr]).tolist()])
        savePY.extend([(dsetXY[savePYstr]).tolist()])
        savePZ.extend([(dsetXY[savePZstr]).tolist()])
        saveE.extend([(dsetXY[saveEstr]).tolist()])
        savePT.extend([(p4_afterCut.pt).tolist()])
        saveEta.extend([(p4_afterCut.eta).tolist()])
        savePhi.extend([(p4_afterCut.phi).tolist()])
        saveM.extend([(dsetXY[saveMstr]).tolist()])
        saveQ.extend([(dsetXY[saveQstr]).tolist()])
        savePDGID.extend([(dsetXY[savePDGIDstr]).tolist()])
        savePartID.extend([(dsetXY[savePartIDstr]).tolist()])
        saveVX.extend([(dsetXY[saveVXstr]).tolist()])
        saveVY.extend([(dsetXY[saveVYstr]).tolist()])
        saveVZ.extend([(dsetXY[saveVZstr]).tolist()])

        counter["n"] += 1
    
    return None

f = h5py.File('data/full_pileup_pilot/ttbar/v2/truth/particles/events0-999.h5', 'r')

fKeys = list(f.keys())

print(f.visit(showContent))

# Sorting particles by pt
for i in range(len(savePX)):
    tuples = list(zip(savePX[i],savePY[i],savePZ[i],saveE[i],savePT[i],saveEta[i],savePhi[i],saveM[i],saveQ[i],savePDGID[i],savePartID[i],saveVX[i],saveVY[i],saveVZ[i]))
    tuples_sorted = sorted(tuples, key = lambda x : x[4])
    savePX[i],savePY[i],savePZ[i],saveE[i],savePT[i],saveEta[i],savePhi[i],saveM[i],saveQ[i],savePDGID[i],savePartID[i],saveVX[i],saveVY[i],saveVZ[i] = zip(*tuples_sorted)

begtime = time.time()

with open('particles.txt', 'w') as f:
    for i in range(len(savePX)):
        f.write(f"particles:{len(savePX[i])}\n")
        for j in range(len(savePX[i])):
            # print(f"{saveVX[i][j]},{saveVY[i][j]},{saveVZ[i][j]},{savePX[i][j]},{savePY[i][j]},{savePZ[i][j]},{saveE[i][j]},{savePT[i][j]},{saveEta[i][j]},{savePhi[i][j]},{saveM[i][j]},{saveQ[i][j]},{savePDGID[i][j]},{savePartID[i][j]}\n")
            f.write(f"{saveVX[i][j]},{saveVY[i][j]},{saveVZ[i][j]},{savePX[i][j]},{savePY[i][j]},{savePZ[i][j]},{saveE[i][j]},{savePT[i][j]},{saveEta[i][j]},{savePhi[i][j]},{saveM[i][j]},{saveQ[i][j]},{savePDGID[i][j]},{savePartID[i][j]}\n")

endtime = time.time()

tottime = endtime - begtime
print(tottime)
