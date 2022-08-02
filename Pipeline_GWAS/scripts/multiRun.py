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
window = int(sys.argv[3])

def prune(dataCut, chrm):
    chromData = pd.DataFrame(dataCut[dataCut["Chr"] == chrm])

    for index, row in chromData.iterrows():
        # Get current row values
        curPVal = row["PVAL"]
        curPos = row["Pos"] 
        curChr = row["Chr"]
        start = curPos - window
        end = curPos + window

        if curPVal <= chromData.loc[(chromData["Pos"] >= start) & (chromData["Pos"] <= end) & (chromData["Chr"] == curChr), "PVAL"].min():
            chromData.loc[index, "IndexSNP"] = 1
    
    windowText = str(window/1000)

    # In case a list of all filtered values per each chromosome is required
    #chromData.to_csv(str(outPath) + "/gwas_" + windowText + "K_" + str(chrm) + ".csv", index=False)

    sortedHits = pd.DataFrame(chromData[chromData["IndexSNP"] == 1])
    sortedHits.to_csv("./gwas_" + windowText + "K_" + str(chrm) + "_Loci.csv", index=False)
    print("Analyzed " + chrm)

if __name__ == "__main__":
    path = pathIn

    data = pd.read_csv(path, sep="\t", header=0)
    dataCut = pd.DataFrame(data[data["PVAL"] <= 5E-8])
    dataCut.sort_values(by=["PVAL", "Pos"], inplace=True)

    # Assign dummy column with no hits
    dataCut["IndexSNP"] = 0

    # list of chromosomes
    chroms = list(dataCut["Chr"].unique())

    N_THREADS = 22
    p = multiprocessing.Pool(N_THREADS)

    func = partial(prune, dataCut)

    results = p.map(func, chroms)


