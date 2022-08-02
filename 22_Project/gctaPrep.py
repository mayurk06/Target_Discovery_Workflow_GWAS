import pandas as pd
import numpy as np
import sys

"""
This script will process 
"""

pathGWAS = "/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"

gwasData = pd.read_csv(pathGWAS, sep="\t", header=0)

toReplace = {"rs77924615": 1.19e-348,
            "rs1617634": 2.62E-325,
            "rs1617984": 2.93E-331,
            "rs1719245": 2.74E-330,
            "rs1719246": 1.63E-370,
            "rs1145097": 9.24E-369,
            "rs1153862": 6.80E-372,
            "rs1153860": 4.18E-374,
            "rs35861938": 7.07E-374,
            "rs2467865": 2.51E-373,
            "rs2461700": 4.43E-374,
            "rs2453533": 9.19E-375,
            "rs2467862": 1.39E-374,
            "rs1145093": 3.68E-374,
            "rs1145089": 7.49E-339,
            "rs1049518": 1.37E-374,
            "rs1145086": 5.80E-373,
            "rs1153857": 6.58E-372,
            "rs1153855": 7.32E-373,
            "rs2486274": 1.90E-374,
            "rs1145084": 1.67E-375,
            "rs1145080": 3.22E-375,
            "rs1145077": 1.47E-375,
            "rs1145076": 2.68E-375,
            "rs1426932": 5.21E-374,
            "rs1346267": 4.00E-374,
            "rs2467853": 5.35E-370,
            "rs2486288": 1.07E-373,
            "rs2433601": 1.56E-373
}

for SNP in toReplace:
    gwasData.loc[gwasData["SNP"] == SNP, "PVAL"] = toReplace[SNP]

window = 500000
# print(gwasData.head(10))

def getSNPList(chrIn, SNP, df):
    chrom = int(chrIn[3:])

    # for each value in the index files, collect all the SNPs that are 
    # 500K +/- index position
    # output single file with just the rsIDs
    indexSNP = gwasData.loc[gwasData["SNP"] == SNP] #filter this to just curChr or above
    # print("indexSNP:",indexSNP)

    curPos = int(indexSNP["Pos"])
    # print("CURPOS:",curPos)

    endPos = curPos + window
    startPos = curPos - window

    locusData = gwasData.loc[(gwasData["Chr"] == chrIn) & (gwasData["Pos"] >= startPos) & (gwasData["Pos"] <= endPos)]
    
    listSNP = list(locusData["SNP"])

    # print(listSNP[0:10])

    with open(path+"/"+SNP+".loci", mode='w+') as outFile:
        outFile.write('\n'.join(listSNP))


getSNPList(str(curChr), str(curSNP), str(outPath))