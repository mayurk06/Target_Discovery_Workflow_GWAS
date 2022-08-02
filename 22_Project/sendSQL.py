import sqlite3

# conn = sqlite3.connect('GTEx.db')
conn = sqlite3.connect('mQTL.db')
cur = conn.cursor()

# query = "CREATE INDEX eIndex_eGeneID_G ON NephQTL_Glomeruli(eGeneID_G);"
query = "CREATE INDEX mId_SNP_N ON meQTL_Summary_Stats_Sig(mSNP_N);"

# query = "CREATE INDEX eId_SNP_X on GTEx_Cortex(SNP);"

cur.executescript(query)

print("DONE mqtlSNP")
