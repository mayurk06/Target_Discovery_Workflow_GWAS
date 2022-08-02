import pandas as pd
import os

QTL_csv_folder = "/home/khannm06/users/khannam/Renal_Project/Scripts/mQTLSignals"

filename_root = "INDEX_meQTL_Blood_SNP"

CSV_filenames = os.listdir(QTL_csv_folder)

print(CSV_filenames)