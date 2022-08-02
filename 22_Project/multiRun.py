import pandas as pd
import numpy as np
import multiprocessing
import sys

from functools import partial

 #TODO: need to add the replacement SNPs from the byChrom file. Look into making this 
 # a class for ease of use.
 #TODO: pass in the windows as a parameter in the parallel function
 #TODO: add option to pass in window size from command line (and then also SNP/file for one single analysis)

pathIn = sys.argv[1]
outPath = sys.argv[2]

def prune(dataCut, chrm):
    chromData = pd.DataFrame(dataCut[dataCut["Chr"] == chrm])

    for index, row in chromData.iterrows():
        # Get current row values
        window = 500000
        curPVal = row["PVAL"]
        curPos = row["Pos"] 
        curChr = row["Chr"]
        start = curPos - window
        end = curPos + window

        if curPVal <= chromData.loc[(chromData["Pos"] >= start) & (chromData["Pos"] <= end) & (chromData["Chr"] == curChr), "PVAL"].min():
            chromData.loc[index, "IndexSNP"] = 1

    # chromData.to_csv(str(outPath) + "/chroms/gwas_" + str(window/1000) + "K_" + str(chrm) + ".csv", index=False)
    chromData.to_csv(str(outPath) + "/gwas_" + str(window/1000) + "K_" + str(chrm) + ".csv", index=False)

    sortedHits = pd.DataFrame(chromData[chromData["IndexSNP"] == 1])
    # sortedHits.to_csv(str(outPath) + "/chroms_hits/gwas_" + str(window/1000) + "K_" + str(chrm) + "_HITS.csv", index=False)
    sortedHits.to_csv(str(outPath) + "/gwas_" + str(window/1000) + "K_" + str(chrm) + "_HITS.csv", index=False)
    print("Analyzed " + chrm)

if __name__ == "__main__":
    # path = "/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"
    path = pathIn

    data = pd.read_csv(path, sep="\t", header=0)
    dataCut = pd.DataFrame(data[data["PVAL"] <= 5E-8])
    dataCut.sort_values(by=["PVAL", "Pos"], inplace=True)

    # # Manually replace values that have a p value of "0" in raw data
    # toReplace = {"rs77924615": 1.19e-348,
    #             "rs1617634": 2.62E-325,
    #             "rs1617984": 2.93E-331,
    #             "rs1719245": 2.74E-330,
    #             "rs1719246": 1.63E-370,
    #             "rs1145097": 9.24E-369,
    #             "rs1153862": 6.80E-372,
    #             "rs1153860": 4.18E-374,
    #             "rs35861938": 7.07E-374,
    #             "rs2467865": 2.51E-373,
    #             "rs2461700": 4.43E-374,
    #             "rs2453533": 9.19E-375,
    #             "rs2467862": 1.39E-374,
    #             "rs1145093": 3.68E-374,
    #             "rs1145089": 7.49E-339,
    #             "rs1049518": 1.37E-374,
    #             "rs1145086": 5.80E-373,
    #             "rs1153857": 6.58E-372,
    #             "rs1153855": 7.32E-373,
    #             "rs2486274": 1.90E-374,
    #             "rs1145084": 1.67E-375,
    #             "rs1145080": 3.22E-375,
    #             "rs1145077": 1.47E-375,
    #             "rs1145076": 2.68E-375,
    #             "rs1426932": 5.21E-374,
    #             "rs1346267": 4.00E-374,
    #             "rs2467853": 5.35E-370,
    #             "rs2486288": 1.07E-373,
    #             "rs2433601": 1.56E-373}

    # for SNP in toReplace:
    #     dataCut.loc[dataCut["SNP"] == SNP, "PVAL"] = toReplace[SNP]

    # Assign dummy column with no hits
    dataCut["IndexSNP"] = 0

    # list of chromosomes
    chroms = list(dataCut["Chr"].unique())

    # Currently have to set window manually above

    # !CHECK YOUR COMPUTING UNITS FIRST!
    N_THREADS = 22
    p = multiprocessing.Pool(N_THREADS)

    func = partial(prune, dataCut)

    results = p.map(func, chroms)


