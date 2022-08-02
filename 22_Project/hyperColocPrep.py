import pandas as pd
import numpy as np
import sys
import sqlite3
import csv

# For each column that has a QTL signal grab the flanked list 
# from the appropriate source
# Can just look for up to 3 for now
# Turn column names into list and then count how many QTL sources 

# SNP = sys.argv[1]
# outPath = sys.argv[2]

# For raw data that doesnt have Pos, use the sum stats or GWAS sum stats
lociPath = "/home/khannm06/users/khannam/Renal_Project/Data/GCTA-slct/Index_Loci_GWAS_GCTA.tsv"
eQTLNephPath = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Results/eQTL_neph.csv"
mQTLNephPath = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Results/mQTL_neph.csv"
mQTLBloodPath = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Results/mQTL_blood.csv"

dataSetPaths = ["Index_Loci_GWAS_GCTA", "eQTL_neph", "mQTL_neph", "mQTL_blood"]

gwasPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"
eQTLSumPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/Kidney_eQTL_Meta_S686_Significant.q0.01.txt"
mQTLSumPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/Kidney_eQTL_Meta_S686_Significant.q0.01.txt"

# QTL/GWAS loci datasets
GWASLoci = pd.read_csv(lociPath, sep="\t", header=0)
eQTLLoci = pd.read_csv(eQTLNephPath, sep=",", header=0)
mQTLNephLoci = pd.read_csv(mQTLNephPath, sep=",", header=0)
mQTLBloodLoci = pd.read_csv(mQTLBloodPath, sep=",", header=0)

GWASLoci.drop_duplicates(subset="SNP", inplace=True)
eQTLLoci.drop_duplicates(subset="SNP", inplace=True)
mQTLNephLoci.drop_duplicates(subset="SNP", inplace=True)
mQTLBloodLoci.drop_duplicates(subset="SNP", inplace=True)

gwasData = pd.read_csv(gwasPath, sep="\t", usecols=["SNP", "Chr", "Pos", "BETA", "SE"], header=0)

# eQTL are al seperate, but mQTL in one
eQTLDBName = "eQTL.db"
mQTLDBName = "mQTL.db"

# print(gwasQTL)
# print("GWAS")
# print(gwasQTL2.head())
# print(gwasQTL2.info())
# exit()
# print(GWASLoci.head())
# print(GWASLoci.info())

window = 500000
c=0
# Match based on matches from output lists
for indexG, rowG in GWASLoci.iterrows():
    gwasSNP = rowG["SNP"]

    for indexQ, rowQ in eQTLLoci.iterrows():
        qtlSNP = rowQ["Index_SNP"]

        if qtlSNP == gwasSNP:
            pass
        else:
            continue

        qtlGeneID = rowQ["GeneID"]

        for indexMn, rowMn in mQTLNephLoci.iterrows():
            mqtlNSNP = rowMn["Index_SNP"]

            if mqtlNSNP == gwasSNP:
                pass
            else:
                continue

            mqtlNcpg = rowMn["CpG"]

            for indexMb, rowMb in mQTLBloodLoci.iterrows():
                mqtlBSNP = rowMb["Index_SNP"]

                if mqtlBSNP == gwasSNP:
                    pass
                else:
                    continue

                mqtlBcpg = rowMb["CpG"]

                if gwasSNP == qtlSNP == mqtlNSNP == mqtlBSNP and mqtlNcpg == mqtlBcpg:
                    gwasCols = ["SNP", "BETA", "SE"]
                    IndexChr = "chr"+str(rowG["Chr"])
                    indexPos = int(rowG["bp"])
                    startPos = indexPos - window
                    endPos   = indexPos + window

                    hyprDataMatch = gwasData.loc[(gwasData["Chr"] == IndexChr) & (gwasData["Pos"] >= startPos) & (gwasData["Pos"] <= endPos)]
                    hyprData = hyprDataMatch.rename(columns={"BETA": "BETA_GWAS", "SE": "SE_GWAS"})

                    # matrixDB = sqlite3.connect('mode=memory')
                    matrixDB = sqlite3.connect(':memory:')
                    matrixDBCur = matrixDB.cursor()

                    # hyprDataMatch.to_sql('gwasIndexTemp', con=matrixDB, if_exists="replace")
                    hyprData.to_sql('gwasIndexTemp', con=matrixDB, if_exists="replace")
                    matrixDBCur.executescript("ATTACH DATABASE 'GTEx.db' AS GTEx;")
                    matrixDBCur.executescript("ATTACH DATABASE 'NephQTL_Glomeruli.db' AS eQTL_G;")
                    matrixDBCur.executescript("ATTACH DATABASE 'NephQTL_Tubule.db' AS eQTL_T;")
                    matrixDBCur.executescript("ATTACH DATABASE 'mQTL.db' AS mQTL;")

                    print(hyprData.head())

                    c+=1
                    # if c>3:
                    #     exit()
                    try:
                        # If table names are unique accross all dbs, then dont need the prefix. format
                        joinTableQuery = '''
                        SELECT 
                            gwasIndexTemp.SNP, 
                            gwasIndexTemp.BETA_GWAS, 
                            gwasIndexTemp.SE_GWAS, 
                            GTEx_Cortex.BETA,
                            GTEx_Cortex.SE,
                            Neph_Glomeruli.BETA, 
                            Neph_Glomeruli.SE,
                            Neph_Tubule.BETA,
                            Neph_Tubule.SE,
                            meQTL_Blood.mBETA_B,
                            meQTL_Blood.mSE_B,
                            meQTL_Summary_Stats_Sig.mBETA_N,
                            meQTL_Summary_Stats_Sig.mSE_N

                        FROM gwasIndexTemp 

                        INNER JOIN GTEx_Cortex
                            ON gwasIndexTemp.SNP = GTEx_Cortex.SNP

                        INNER JOIN Neph_Glomeruli
                            ON gwasIndexTemp.SNP = Neph_Glomeruli.SNP
                        
                        INNER JOIN Neph_Tubule
                            ON gwasIndexTemp.SNP = Neph_Tubule.SNP
                        
                        LEFT JOIN meQTL_Blood
                            ON gwasIndexTemp.SNP = meQTL_Blood.mSNP_B AND
                            meQTL_Blood.mCpG_B = '$M'
                        
                        LEFT JOIN meQTL_Summary_Stats_Sig
                            ON gwasIndexTemp.SNP = meQTL_Summary_Stats_Sig.mSNP_N AND
                            meQTL_Summary_Stats_Sig.mCpG_N = '$M'

                        WHERE Neph_Tubule.GeneID = '$G' AND Neph_Glomeruli.GeneID = '$G' AND GTEx_Cortex.GeneID = '$G';'''


                        sqlQuery = joinTableQuery.replace("$G", str(qtlGeneID))
                        sqlQuery = sqlQuery.replace("$M", str(mqtlBcpg))

                        # print(sqlQuery)
                        # exit()
                        csvName = "hyperDataFull/"+str(gwasSNP)+"_"+str(qtlGeneID)+"_"+str(mqtlBcpg)+".csv"

                        db_df = pd.read_sql_query(sqlQuery, matrixDB)
                        db_df.to_csv(csvName, index=False)

                        print("Ran SQL Query")

                    except Exception as e:
                        print(e) 
                        continue

exit()


for index, row in gwasQTL.iterrows():
    IndexSNP = row["Index_SNP"]
    eSNPGeneID = row["GeneID"]
    mCpG = row["CpG"]

    IndexChr = "chr"+str(row["ChrNone"])
    indexPos = int(row["bp"])
    startPos = indexPos - window
    endPos   = indexPos + window

    print(IndexSNP, indexPos, IndexChr, eSNPGeneID, mCpG)

    gwasCols = ["SNP", "BETA", "SE"]
    hyprDataMatch = gwasData.loc[(gwasData["Chr"] == IndexChr) & (gwasData["Pos"] >= startPos) & (gwasData["Pos"] <= endPos)]
    hyprData = hyprDataMatch.rename(columns={"BETA": "BETA_GWAS", "SE": "SE_GWAS"})

    # print(hyprData.head())
    print("222")
    # matrixDB = sqlite3.connect('mode=memory')
    matrixDB = sqlite3.connect(':memory:')
    matrixDBCur = matrixDB.cursor()

    # hyprDataMatch.to_sql('gwasIndexTemp', con=matrixDB, if_exists="replace")
    hyprData.to_sql('gwasIndexTemp', con=matrixDB, if_exists="replace")
    print("333")
    matrixDBCur.executescript("ATTACH DATABASE 'eQTL.db' AS eQTL;")
    matrixDBCur.executescript("ATTACH DATABASE 'mQTL.db' AS mQTL;")
    print("444")
    c+=1
    if c>4:
        exit()
    try:
        # c+=1
        # if c>4:
        #     exit()

        # more verbose query generation: why inner join, when can do conditional select?

#SELECT gwasIndexTemp.SNP, gwasIndexTemp.BETA_GWAS, gwasIndexTemp.SE_GWAS, eQTL.NephQTL_Glomeruli.eSNP_G, eQTL.NephQTL_Glomeruli.eGeneID_G, eQTL.NephQTL_Glomeruli.eBETA_G, eQTL.NephQTL_Glomeruli.eSE_G
        joinTableQuery = '''SELECT gwasIndexTemp.SNP, gwasIndexTemp.BETA_GWAS, gwasIndexTemp.SE_GWAS, eQTL.NephQTL_Glomeruli.eBETA_G, eQTL.NephQTL_Glomeruli.eSE_G
            FROM gwasIndexTemp 
            INNER JOIN eQTL.NephQTL_Glomeruli
                ON gwasIndexTemp.SNP = eQTL.NephQTL_Glomeruli.eSNP_G
            INNER JOIN eQTL.NephQTL_Tubule
                ON gwasIndexTemp.SNP = eQTL.NephQTL_Tubule.eSNP_T
                WHERE eQTL.NephQTL_Tubule.eGeneID_T = '$G' OR eQTL.NephQTL_Glomeruli.eGeneID_G = '$G';'''

        print("Test")


        sqlQuery = joinTableQuery.replace("$G", str(eSNPGeneID))
        sqlQuery = sqlQuery.replace("$M", str(mCpG))
        # print(sqlQuery)
       
        csvName = "hyperData/"+str(IndexSNP)+"_"+str(eSNPGeneID)+"_"+str(mCpG)+".csv"

        db_df = pd.read_sql_query(sqlQuery, matrixDB)
        db_df.to_csv(csvName, index=False)

        # matrixDBCur.execute(sqlQuery)

        # with open(csvName, "w", newline='') as csv_file:   
        #     csv_writer = csv.writer(csv_file)
        #     csv_writer.writerow([i[0] for i in matrixDBCur.description]) # write headers
        #     csv_writer.writerows(matrixDBCur)

        #     for row in matrixDBCur:
        #         print(row)
        #         csv_writer.writerow(row)
    
    except Exception as e: 
        print(e)
        continue

    c+=1


exit()

# SQL multiple INNER JOIN queries:
# FROM gwasIndexTemp 
# INNER JOIN eQTL.NephQTL_Glomeruli
#     ON gwasIndexTemp.SNP = eQTL.NephQTL_Glomeruli.eSNP_G
#     WHERE eQTL.NephQTL_Glomeruli.eGeneID_G = '$G'
# INNER JOIN eQTL.NephQTL_Tubule
#     ON gwasIndexTemp.SNP = eQTL.NepthQTL_Tubule.eSNP_T
#     WHERE eQTL.NepthQTL_Tubule.eGeneID_T = 'G;'''

# for index, row in GWASLoci.iterrows():
#     IndexSNP = row["SNP"]

#     for indexE, rowE in eQTLLoci.iterrows():
#         eSNPIndex = rowE["Index_SNP"]
#         eSNPGeneID = rowE["GeneID"]

#         for indexM, rowM in mQTLNephLoci.iterrows():
#             mSNPIndex = rowM["Index_SNP"]
#             mCpG = rowM["CpG"]

#             # Left join x2 the gwas and qtl pd.dfs and use that list to run through 
#             # the values
#             if IndexSNP == eSNPIndex == mSNPIndex:
#                 IndexChr = "chr" + str(row["Chr"])
#                 indexPos = int(row["bp"])
#                 startPos = indexPos - window
#                 endPos   = indexPos + window

#                 print(IndexSNP, eSNPGeneID, eSNPGeneID, mSNPIndex, mCpG)

#                 gwasCols = ["SNP", "BETA", "SE"]
#                 hyprDataMatch = gwasData.loc[(gwasData["Chr"] == IndexChr) & (gwasData["Pos"] >= startPos) & (gwasData["Pos"] <= endPos)]
#                 hyprData = hyprDataMatch.rename(columns={"BETA": "BETA_GWAS", "SE": "SE_GWAS"})

#                 matrixDB = sqlite3.connect('file:cachedb?mode=memory')
#                 matrixDBCur = matrixDB.cursor()

#                 hyprData.to_sql('gwasIndexTemp', con=matrixDB, if_exists="replace")

#                 matrixDBCur.executescript("ATTACH DATABASE 'eQTL.db' AS eQTL;")
#                 matrixDBCur.executescript("ATTACH DATABASE 'meQTL.db' AS mQTL;")

#                 try:
#                     joinTableQuery = '''SELECT SNP, BETA_GWAS, SE_GWAS 
#                         FROM gwasIndexTemp 
#                         INNER JOIN eQTL.NephQTL_Glomeruli
#                             ON gwasIndexTemp.SNP = eQTL.NephQTL_Glomeruli.eSNP_G
#                         WHERE eQTL.NephQTL_Glomeruli.eGeneID_G = $G;''' #change to proper

#                     matrixDBCur.execute(joinTableQuery.replace("$G", str(eSNPGeneID)))

#                     csvName = str(IndexSNP)+"_"+str(eSNPIndex)+"-"+str(eSNPGeneID)+"_"+str(mSNPIndex)+"_"+str(mCpG)+".csv"
#                     with open(csvName, "w", newline='') as csv_file:   
#                         csv_writer = csv.writer(csv_file)
#                         csv_writer.writerow([i[0] for i in matrixDBCur.description]) # write headers
#                         csv_writer.writerows(matrixDBCur)

#                     rows = c.fetchall()
#                     for row in rows:
#                         # do your stuff
#                         csvWriter.writerow(row)

#                     print(result)
                
#                 except: 
#                     continue

#                 c+=1
#                 if c>1:
#                     exit()




    # set equal to values that are associated w/ index SNP and then get gene to 
    # search by SNP + GeneID

    # For each gwas index if index snp exists in all the qtl sources, 
    # then collect SNPs for each of the qtl sources, loop through every trait 


