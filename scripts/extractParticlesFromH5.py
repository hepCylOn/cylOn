import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import vector

# Added counter to run over a given amount of events for debugging
counter = {"n": 0}
debug = False

# Output file name
outFile = 'particles'

# Add option to only save muons (goood for debugging)
muonsOnly = False
if muonsOnly: outFile = outFile + "MuonsOnly"

# Since the particles information is appended to the outFile, it is better
# to remove it before executing the code
if os.path.exists(outFile + '.txt'): os.system('rm ' + outFile + '.txt')

# Function that accesses the information inside of the .h5 file
def showContent(name):
    print("=============================================")
    print("Started function!!")
    if debug:
        if counter["n"] >= 1: return True
    global f
    obj = f[name]
    if isinstance(obj, h5py.Dataset):
        time3 = time.time()
        # Gets names of variables if they are lower or upper case; might be useful
        # when using this code for other files or other file types
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

        time4 = time.time()

        # Apply some cuts to mimic tracking studies at CMS
        # Only save particles with zVertex < 30 cm
        vzCheck = (abs(obj[()][saveVZstr]) < (300.0*np.ones(len(obj[()][saveVZstr]))))
        dataset =  obj[()][vzCheck]

        # Only save particles with rVertex < 2.5 cm
        vrCheck = (np.sqrt((dataset[saveVXstr]*dataset[saveVXstr])+(dataset[saveVYstr]*dataset[saveVYstr])) < (25.0*np.ones(len(dataset[saveVXstr]))))
        dataset = dataset[vrCheck]

        if muonsOnly:
            # Only save muons
            pdgidCheck = (abs(dataset[savePDGIDstr]) == (13*np.ones(len(dataset[savePDGIDstr]))))
            dataset = dataset[pdgidCheck]

        p4_beforeCut = vector.arr({"px": dataset[savePXstr], "py": dataset[savePYstr], "pz": dataset[savePZstr], "E": dataset[saveEstr]})

        # Only save particles with pT > 0.9 GeV
        pTCheck = (p4_beforeCut.pt > (0.9*np.ones(len(dataset[savePXstr]))))
        dataset = dataset[pTCheck]

        p4_afterPtCut = vector.arr({"px": dataset[savePXstr], "py": dataset[savePYstr], "pz": dataset[savePZstr], "E": dataset[saveEstr]})

        # Only save particles that usually pass through layer 0, i.e., |eta| < 3.2
        etaCheck = (abs(p4_afterPtCut.eta) < (3.2*np.ones(len(dataset[savePXstr]))))
        dataset = dataset[etaCheck]

        p4_afterCuts = vector.arr({"px": dataset[savePXstr], "py": dataset[savePYstr], "pz": dataset[savePZstr], "E": dataset[saveEstr]})

        time5 = time.time()

        # Convert columns of dataset to be more readable
        datasetPxCol = (dataset[savePXstr]).tolist()
        datasetPyCol = (dataset[savePYstr]).tolist()
        datasetPzCol = (dataset[savePZstr]).tolist()
        datasetECol = (dataset[saveEstr]).tolist()
        datasetPTCol = (p4_afterCuts.pt).tolist()
        datasetEtaCol = (p4_afterCuts.eta).tolist()
        datasetPhiCol = (p4_afterCuts.phi).tolist()
        datasetMCol = (dataset[saveMstr]).tolist()
        datasetQCol = (dataset[saveQstr]).tolist()
        datasetPDGIDCol = (dataset[savePDGIDstr]).tolist()
        datasetPartIDCol = (dataset[savePartIDstr]).tolist()
        datasetVxCol = (dataset[saveVXstr]).tolist()
        datasetVyCol = (dataset[saveVYstr]).tolist()
        datasetVzCol = (dataset[saveVZstr]).tolist()

        # Check to see if the lists exist. Even if 0 particles, count the event to match with hits file
        if len(datasetPxCol) == 0:
            with open(outFile + '.txt', 'a') as fout:
                fout.write(f"particles:0\n")
            counter["n"] += 1
            print(counter["n"])
            return None

        # To apply the particles' pT ordering to all of the information to be saved, a tuple has to be
        # created so that the ordering acts in the same way over all columns. Afterwards, the tuple
        # can be split. Index 4 in x[4] means the pT column
        tuples = list(zip(datasetPxCol,datasetPyCol,datasetPzCol,datasetECol,datasetPTCol,datasetEtaCol,datasetPhiCol,datasetMCol,datasetQCol,datasetPDGIDCol,datasetPartIDCol,datasetVxCol,datasetVyCol,datasetVzCol))
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

        # Information per event is appended to outFile as:
        # particles:total_number_of_particles
        # vx,vy,vz,px,py,pz,E,pT,eta,phi,m,q,pdgID,partID
        with open(outFile + '.txt', 'a') as fout:
            fout.write(f"particles:{len(savePX)}\n")
            for j in range(len(savePX)):
                fout.write(f"{saveVX[j]},{saveVY[j]},{saveVZ[j]},{savePX[j]},{savePY[j]},{savePZ[j]},{saveE[j]},{savePT[j]},{saveEta[j]},{savePhi[j]},{saveM[j]},{saveQ[j]},{savePDGID[j]},{savePartID[j]}\n")
            fout.close()

        time8 = time.time()
        # Checks execution time in distinct blocks:
        # - To get names (very fast -> 0.2 ms)
        # - To perform cuts (slow -> 52 ms)
        # - To sort by pT (very fast -> 1 ms)
        # - To convert to lists (very fast -> 0.06 ms)
        # - To write to file (slow -> 69 ms)
        # Still process about 10 ev/s which is manageable
        print(f"Checking times: to get object names {time4 - time3} -- to perform cuts {time5 - time4} -- to sort by pt {time6 - time5} -- to convert to lists {time7 - time6} -- to write to file {time8 - time7}")

        counter["n"] += 1
        print(counter["n"])
    
    return None

f = h5py.File('/home/breno/data/full_pileup_pilot/ttbar/v2/truth/particles/events0-999.h5', 'r')

fKeys = list(f.keys())

f.visit(showContent)