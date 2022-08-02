import sqlite3
import pandas as pd

dbName = "NephQTL_Glomeruli.db"
connection = sqlite3.connect(dbName)

eQTL1Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/GTEx.v8.hg38.Kidney_Cortex.eQTL.txt"
eQTL2Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/NephQTL_Glomeruli.txt"
eQTL3Path = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/NephQTL_Tubule.txt"

mQTLKidneySumPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/Kidney_meQTL_S443_Significant.q0.01.txt"
mQTLBloodSumPath = "/hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/Kidney_mqtl_FinalSum.txt"

pathSelect = eQTL2Path
tableName = "Neph_Glomeruli"
primaryKey = "SNP"

# load the data into a Pandas DataFrame
QTLdf = pd.read_csv(pathSelect, sep="\t", header=0)
QTLdf.set_index(primaryKey)
# multi = df.set_index([primaryKey, 'GeneID']) # multiindex but idt this is working

print(QTLdf.info())
# QTLdf.columns = QTLdf.columns.to_flat_index()

# write the data to a sqlite table

QTLdf.to_sql(tableName, con=connection, if_exists='replace', index=True)

connection.close()
