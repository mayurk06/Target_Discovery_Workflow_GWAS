import pandas as pd
import numpy as np
import re

masterOutPath = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Results/masterLoci.csv"

data = pd.read_csv(masterOutPath, header=0)

cur = data.iloc[229]
print(cur)
print()

print(data.info())
curBP = cur["bp"]
curChr = cur["Chr"]
curSNP = cur["SNP"]
curGene = cur["GENE1"]
print(curBP, curChr, curSNP)
# Prepare data for hyprcoloc
exit()

# From the masterLoci csv, get rsID, pos, chr and 
# then search through the raw original tables for each and
# get the betas and se for all the rsIDs within 500kbp+/-

#^^ instead of basing it off of the bp, use the gwas
# data to get list of rsIDs, then find that subset from
# QTL data

# T1, T2, ...
# rsID val val ...

gwasPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"
eQTL1Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/GTEx.v8.hg38.Kidney_Cortex.eQTL.txt"
eQTL2Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/NephQTL_Glomeruli.txt"
eQTL3Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/NephQTL_Tubule.txt"

# GTEx - !!!20GB file!!! ~5 mins load time
eQTL1 = pd.read_csv(eQTL1Path, sep="\t", header=0)
gwasData = pd.read_csv(gwasPath, sep="\t", header=0)

# Do not use this method
# eQTL1["bp"] = eQTL1["MarkerID"].str.extract(r'(_[0-9]+)', expand=False)

print(eQTL1.head(5))

window = 500000
def exportSNPList(chrIn, SNP):
    # chrom = int(chrIn[3:])
    chrom = chrIn

    # indexSNP = gwasData.loc[gwasData["SNP"] == SNP] # used to get bp position (might not need)
    # print("indexSNP:",indexSNP)

    curPos = curBP
    # print("CURPOS:",curPos)

    endPos = curPos + window
    startPos = curPos - window

    locusData = gwasData.loc[(gwasData["Chr"] == chrIn) & (gwasData["Pos"] >= startPos) & (gwasData["Pos"] <= endPos)]
    
    listSNP = list(locusData["SNP"])
    print(listSNP[0:10])

    return listSNP

# curLoci = exportSNPList(curChr, curSNP)

# eQTL_Loci = eQTL1.loc[eQTL1['SNP'].isin(curLoci)]

# print(eQTL_Loci.head(10))
