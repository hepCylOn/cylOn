import os

import numpy as np
import polars as pl
import vector
import argparse

from colliderml.core import load_tables, collect_tables
from colliderml.polars import explode_particles

pl.Config.set_tbl_rows(12)
pl.Config.set_tbl_cols(20)

parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

parser.add_argument('-o', '--outFileName', default='')
parser.add_argument('-m', '--isMuonsOnly', action='store_true')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-e', '--eventsForDebug', type=int, default=1)
parser.add_argument('-n', '--numberOfEvents', type=int, default=10)

args = parser.parse_args()

muonModifier = ''
if args.isMuonsOnly: muonModifier = 'MuonsOnly'

# Output filename
partFileName = 'particles' + args.outFileName + muonModifier + '.txt'

# Since the particles information is appended to the partFileName, it is better
# to remove it before executing the code
if os.path.exists(partFileName): os.system('rm ' + partFileName)

home_directory = os.path.expanduser("~")

cfg = {
    "dataset_id": "CERN/ColliderML-Release-1",
    "channels": "ttbar",
    "pileup": "pu200",
    "objects": ["particles"],
    "split": "train",
    "lazy": True,
    "max_events": args.numberOfEvents,
    "data_dir": home_directory + "/.cache/colliderml"
}

tables = load_tables(cfg)
frames = collect_tables(tables)

particles_evt = frames["particles"]

particles_mult = particles_evt.select(
    pl.col("particle_id").list.len().alias("n_particles"),
)
particles_flat = explode_particles(particles_evt)

pn = (particles_mult["n_particles"]).to_numpy().astype(int)
ppx = particles_flat["px"].to_numpy().astype(float)
ppy = particles_flat["py"].to_numpy().astype(float)
ppz = particles_flat["pz"].to_numpy().astype(float)
pvx = particles_flat["vx"].to_numpy().astype(float)
pvy = particles_flat["vy"].to_numpy().astype(float)
pvz = particles_flat["vz"].to_numpy().astype(float)
pe = particles_flat["energy"].to_numpy().astype(float)
pm = particles_flat["mass"].to_numpy().astype(float)
pc = particles_flat["charge"].to_numpy().astype(int)
pid = particles_flat["particle_id"].to_numpy().astype(int)
ppid = particles_flat["pdg_id"].to_numpy().astype(int)

counter = {"n": 0}
nAux = 0

for i in range(len(pn)):

    if args.debug:
        if counter["n"] >= args.eventsForDebug: 
            print(f"Created file {partFileName} WHILE DEBUGGING")
            exit()

    ppxAux = ppx[nAux:(nAux + pn[i])]
    ppyAux = ppy[nAux:(nAux + pn[i])]
    ppzAux = ppz[nAux:(nAux + pn[i])]
    pvxAux = pvx[nAux:(nAux + pn[i])]
    pvyAux = pvy[nAux:(nAux + pn[i])]
    pvzAux = pvz[nAux:(nAux + pn[i])]
    peAux = pe[nAux:(nAux + pn[i])]
    pmAux = pm[nAux:(nAux + pn[i])]
    pcAux = pc[nAux:(nAux + pn[i])]
    pidAux = pid[nAux:(nAux + pn[i])]
    ppidAux = ppid[nAux:(nAux + pn[i])]

    nAux += pn[i]

    # Apply some cuts to mimic tracking studies at CMS
    # Only save particles with zVertex < 30 cm
    vzCheck = (abs(pvzAux) < (300.0*np.ones(len(pvzAux))))

    ppxAux = ppxAux[vzCheck]
    ppyAux = ppyAux[vzCheck]
    ppzAux = ppzAux[vzCheck]
    pvxAux = pvxAux[vzCheck]
    pvyAux = pvyAux[vzCheck]
    pvzAux = pvzAux[vzCheck]
    peAux = peAux[vzCheck]
    pmAux = pmAux[vzCheck]
    pcAux = pcAux[vzCheck]
    pidAux = pidAux[vzCheck]
    ppidAux = ppidAux[vzCheck]

    # Only save particles with rVertex < 2.5 cm
    vrCheck = (np.sqrt((pvxAux*pvxAux)+(pvyAux*pvyAux)) < (25.0*np.ones(len(pvxAux))))

    ppxAux = ppxAux[vrCheck]
    ppyAux = ppyAux[vrCheck]
    ppzAux = ppzAux[vrCheck]
    pvxAux = pvxAux[vrCheck]
    pvyAux = pvyAux[vrCheck]
    pvzAux = pvzAux[vrCheck]
    peAux = peAux[vrCheck]
    pmAux = pmAux[vrCheck]
    pcAux = pcAux[vrCheck]
    pidAux = pidAux[vrCheck]
    ppidAux = ppidAux[vrCheck]

    if args.isMuonsOnly:
        # Only save muons
        pdgidCheck = (abs(ppidAux) == (13*np.ones(len(ppidAux))))

        ppxAux = ppxAux[pdgidCheck]
        ppyAux = ppyAux[pdgidCheck]
        ppzAux = ppzAux[pdgidCheck]
        pvxAux = pvxAux[pdgidCheck]
        pvyAux = pvyAux[pdgidCheck]
        pvzAux = pvzAux[pdgidCheck]
        peAux = peAux[pdgidCheck]
        pmAux = pmAux[pdgidCheck]
        pcAux = pcAux[pdgidCheck]
        pidAux = pidAux[pdgidCheck]
        ppidAux = ppidAux[pdgidCheck]

    p4_beforeCut = vector.arr({"px": ppxAux, "py": ppyAux, "pz": ppzAux, "E": peAux})

    # Only save particles with pT > 0.9 GeV
    pTCheck = (p4_beforeCut.pt > (0.9*np.ones(len(ppxAux))))

    ppxAux = ppxAux[pTCheck]
    ppyAux = ppyAux[pTCheck]
    ppzAux = ppzAux[pTCheck]
    pvxAux = pvxAux[pTCheck]
    pvyAux = pvyAux[pTCheck]
    pvzAux = pvzAux[pTCheck]
    peAux = peAux[pTCheck]
    pmAux = pmAux[pTCheck]
    pcAux = pcAux[pTCheck]
    pidAux = pidAux[pTCheck]
    ppidAux = ppidAux[pTCheck]

    p4_afterPtCut = vector.arr({"px": ppxAux, "py": ppyAux, "pz": ppzAux, "E": peAux})

    # Only save particles that usually pass through layer 0, i.e., |eta| < 3.2
    etaCheck = (abs(p4_afterPtCut.eta) < (3.2*np.ones(len(ppxAux))))

    ppxAux = ppxAux[etaCheck]
    ppyAux = ppyAux[etaCheck]
    ppzAux = ppzAux[etaCheck]
    pvxAux = pvxAux[etaCheck]
    pvyAux = pvyAux[etaCheck]
    pvzAux = pvzAux[etaCheck]
    peAux = peAux[etaCheck]
    pmAux = pmAux[etaCheck]
    pcAux = pcAux[etaCheck]
    pidAux = pidAux[etaCheck]
    ppidAux = ppidAux[etaCheck]

    p4_afterCuts = vector.arr({"px": ppxAux, "py": ppyAux, "pz": ppzAux, "E": peAux})
    
    # Convert columns of dataset to be more readable
    datasetPxCol = (ppxAux).tolist()
    datasetPyCol = (ppyAux).tolist()
    datasetPzCol = (ppzAux).tolist()
    datasetECol = (peAux).tolist()
    datasetPTCol = (p4_afterCuts.pt).tolist()
    datasetEtaCol = (p4_afterCuts.eta).tolist()
    datasetPhiCol = (p4_afterCuts.phi).tolist()
    datasetMCol = (pmAux).tolist()
    datasetQCol = (pcAux).tolist()
    datasetPDGIDCol = (ppidAux).tolist()
    datasetPartIDCol = (pidAux).tolist()
    datasetVxCol = (pvxAux/10.0).tolist()
    datasetVyCol = (pvyAux/10.0).tolist()
    datasetVzCol = (pvzAux/10.0).tolist()

    # Check to see if the lists exist. Even if 0 particles, count the event to match with hits file
    if len(datasetPxCol) == 0:
        with open(partFileName, 'a') as fout:
            fout.write(f"particles:0\n")
        counter["n"] += 1
        continue

    # To apply the particles' pT ordering to all of the information to be saved, a tuple has to be
    # created so that the ordering acts in the same way over all columns. Afterwards, the tuple
    # can be split. Index 4 in x[4] means the pT column
    tuples = list(zip(datasetPxCol,datasetPyCol,datasetPzCol,datasetECol,datasetPTCol,datasetEtaCol,datasetPhiCol,datasetMCol,datasetQCol,datasetPDGIDCol,datasetPartIDCol,datasetVxCol,datasetVyCol,datasetVzCol))
    tuples_sorted = sorted(tuples, key = lambda x : x[4])
    savePX,savePY,savePZ,saveE,savePT,saveEta,savePhi,saveM,saveQ,savePDGID,savePartID,saveVX,saveVY,saveVZ = zip(*tuples_sorted)
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

    # Information per event is appended to partFileName as:
    # particles:total_number_of_particles
    # vx,vy,vz,px,py,pz,E,pT,eta,phi,m,q,pdgID,partID
    with open(partFileName, 'a') as fout:
        fout.write(f"particles:{len(savePX)}\n")
        for j in range(len(savePX)):
            fout.write(f"{saveVX[j]},{saveVY[j]},{saveVZ[j]},{savePX[j]},{savePY[j]},{savePZ[j]},{saveE[j]},{savePT[j]},{saveEta[j]},{savePhi[j]},{saveM[j]},{saveQ[j]},{savePDGID[j]},{savePartID[j]}\n")
        fout.close()

    if args.debug: print(counter["n"])

    counter["n"] += 1

print(f"Created files {partFileName} with {args.numberOfEvents} events")