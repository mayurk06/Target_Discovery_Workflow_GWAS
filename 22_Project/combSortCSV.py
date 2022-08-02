import pandas as pd

path = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/IndexSNPs_Nearest_Gene/comb_unsorted.csv"

data = pd.read_csv(path, header=0)
# Locus***ClosestGene_Chr-Pos_rsID (e.g. Locus00_PRKCZ_chr1-2049599_rs3107148)
data["ChrNum"] = data["Chr"].str[3:]

data["ChrNum"] = pd.to_numeric(data["ChrNum"])

data.sort_values(by=["ChrNum", "Pos"], inplace=True)

data.insert(0, 'LocusID', range(1, len(data)+1))
# data['LocusID'] = data['LocusID'] + 1
# data.reindex(copy=False)
# data['LocusID'] = data['LocusID'].apply(lambda x: "Locus"+'{0:0>3}'.format(x)+"_"
#                                                     +str(data.loc[x-1, "ClosestGene"])+"_"
#                                                     +str(data.loc[x-1, "Chr"])+"-"
#                                                     +str(data.loc[x-1, "Pos"])+"_"
#                                                     +str(data.loc[x-1, "SNP"]))

data['LocusID'] = data['LocusID'].apply(lambda x: "Locus"+'{0:0>3}'.format(x))
data['LocusID'] = data['LocusID']+"_"+data["ClosestGene"].astype(str)+"_"+data["Chr"].astype(str)+"-"+data["Pos"].astype(str)+"_"+data["SNP"].astype(str)

data.to_csv("/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/IndexSNPs_Nearest_Gene/GWAS_500K_Loci.csv", index=False)

# Add gcta results
gctaCounts = "/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/Loci/gctaCountsByLoci.csv"

counts = pd.read_csv(gctaCounts)

mergedData = data.merge(counts, how="left",on="SNP")

mergedData.to_csv("/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/IndexSNPs_Nearest_Gene/GWAS_500K_Loci.csv", index=False)

