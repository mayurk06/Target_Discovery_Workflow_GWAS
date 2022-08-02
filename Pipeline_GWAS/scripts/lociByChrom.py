import pandas as pd
import numpy as np
import sys

# path = "/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas"
path = sys.argv[1]
windows = [500000, 1000000]
# windows = [500000]
# windows = [sys.argv[1]]

data = pd.read_csv(path, sep="\t", header=0)

# Manually replace values that have a p value of "0" in raw data
# toReplace = {"rs77924615": 1.19e-348,
# 			"rs1617634": 2.62E-325,
# 			"rs1617984": 2.93E-331,
# 			"rs1719245": 2.74E-330,
# 			"rs1719246": 1.63E-370,
# 			"rs1145097": 9.24E-369,
# 			"rs1153862": 6.80E-372,
# 			"rs1153860": 4.18E-374,
# 			"rs35861938": 7.07E-374,
# 			"rs2467865": 2.51E-373,
# 			"rs2461700": 4.43E-374,
# 			"rs2453533": 9.19E-375,
# 			"rs2467862": 1.39E-374,
# 			"rs1145093": 3.68E-374,
# 			"rs1145089": 7.49E-339,
# 			"rs1049518": 1.37E-374,
# 			"rs1145086": 5.80E-373,
# 			"rs1153857": 6.58E-372,
# 			"rs1153855": 7.32E-373,
# 			"rs2486274": 1.90E-374,
# 			"rs1145084": 1.67E-375,
# 			"rs1145080": 3.22E-375,
# 			"rs1145077": 1.47E-375,
# 			"rs1145076": 2.68E-375,
# 			"rs1426932": 5.21E-374,
# 			"rs1346267": 4.00E-374,
# 			"rs2467853": 5.35E-370,
# 			"rs2486288": 1.07E-373,
# 			"rs2433601": 1.56E-373
# }

# for SNP in toReplace:
# 	data.loc[data["SNP"] == SNP, "PVAL"] = toReplace[SNP]

# Filter and sort SNPs with p-val <= 5-E8
dataCut = pd.DataFrame(data[data["PVAL"] <= 5E-8])
dataCut.sort_values(by=["PVAL", "Pos"], inplace=True)

# Assign dummy column with no hits
dataCut["Targets"] = 0

for window in windows:
	testD = dataCut.loc[:]
	# testD = testT.head(5000)
	for index, row in testD.iterrows():
		# Get current row values
		curPVal = row["PVAL"]
		curPos = row["Pos"] 
		curChr = row["Chr"]
		start = curPos - window
		end = curPos + window

		if curPVal <= testD.loc[(testD["Pos"] >= start) & (testD["Pos"] <= end) & (testD["Chr"] == curChr), "PVAL"].min():
			testD.loc[index, "Targets"] = 1

		# testD.to_csv("./gwas_" + str(window/1000) + "K_"+ chrom.upper() + "_.csv", index=False)
		testD.to_csv(sys.argv[2] + "/gwas_" + str(window/1000) + "K.csv", index=False)