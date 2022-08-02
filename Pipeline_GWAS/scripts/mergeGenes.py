import pandas as pd
import numpy as np
import sys

# Get list of chromosomes
# rawDataPath = "/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"
rawDataPath = genePath = sys.argv[1]
dataCut = pd.read_csv(rawDataPath, sep="\t", header=0)
# dataCut = pd.DataFrame(data[data["PVAL"] <= 5E-8])
chroms = list(dataCut["Chr"].unique())

# Gene mapping file
# genePath = "/home/khannm06/renal/gene_pos_hg19_coding_gene_release_87.txt"
genePath = sys.argv[2]
geneList = pd.read_csv(genePath, sep="\t", header=0)

# Get input path
targsIn = sys.argv[3]

# windows = [500000, 1000000]
windows = [500000]

for window in windows:
    for chrom in chroms:
        # Convert chromosome string to int
        curChrom = int(chrom[3:])
        try:
            # targsRoot = "/home/khannm06/renal/parallel/chroms_hits/gwas_"+ str(window/1000) +"K_chrCHROM_HITS.csv"
            targsRoot = str(targsIn) + "/gwas_"+ str(window/1000) +"K_chrCHROM_Loci.csv"
            targsPath = targsRoot.replace("CHROM", str(curChrom)) 

            targsData = pd.read_csv(targsPath, header=0)
            targsData.sort_values(by=["Pos"], inplace=True)
            targsData["ClosestGene"] = "NA"
            targsData["ClosestGeneID"] = 0
            targsData["Delta"] = -1

        except:
            print("passed")
            continue

        # Filter and reference gene list by current chromosome and position
        geneListChrom = geneList.loc[(geneList["chrom"] == curChrom)]
        geneListChrom.sort_values(by=["start"])
        
        print(chrom)
        
        for index1, row1 in targsData.iterrows():
            curPos = row1["Pos"]
            nearest = 9E100
            curIndex = None
            curDelta = -1
            start, end = 0, 0

            # Find closest gene
            for index2, row2 in geneListChrom.iterrows():
                refPosStart = row2["start"]
                refPosEnd = row2["end"]

                deltaStart = abs(curPos - refPosStart)
                deltaEnd = abs(curPos - refPosEnd)

                if deltaStart < nearest:
                    nearest = deltaStart
                    curDelta = deltaStart
                    curIndex = index2
                    start = refPosStart
                    end = refPosEnd

                if deltaEnd < nearest:
                    nearest = deltaEnd
                    curDelta = deltaEnd
                    curIndex = index2
                    start = refPosStart
                    end = refPosEnd

            if start <= curPos <= end:
                curDelta = 0

            # Once nearest gene has been selected, pick the largest window gene if there is an isoform.
            if curDelta == 0:
                closestGene = geneListChrom.loc[curIndex, "gene_name"]
                print(closestGene)
                # targsData.loc[index1, "ClosestGene"] = geneListChrom.loc[nearest[1], "gene_name"]
                isoforms = geneListChrom.loc[geneListChrom["gene_name"] == closestGene]
                print(isoforms)

                for index, row in isoforms.iterrows():
                    refPosStart = row2["start"]
                    refPosEnd = row2["end"]
                    isoDelta = abs(refPosStart - refPosEnd)

                    if isoDelta > curDelta:
                        curIndex = index
            
            # targsData.loc[index1, "ClosestGene"] = geneListChrom.loc[nearest, "gene_name"]
            # targsData.loc[index1, "Delta"] = curDelta

            if curIndex:
                targsData.loc[index1, "ClosestGene"] = geneListChrom.loc[curIndex, "gene_name"]
                targsData.loc[index1, "ClosestGeneID"] = geneListChrom.loc[curIndex, "gene_id"]
                targsData.loc[index1, "Delta"] = curDelta

        # outPath = targsPath.replace("chroms_hits", "chroms_hits_genes")
        # outPath = outPath.replace("_HITS", "_HITS_GENES")
        outPath = targsPath.replace("_HITS", "_HITS_GENES")
        targsData.to_csv(outPath, index=False)


