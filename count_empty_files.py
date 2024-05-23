import csv
import glob
import pandas as pd

files = glob.glob("results/successful_mutation/*")

print(len(files))

count_non_empty = 0

for i in range(len(files)):
    peptide_data = pd.read_csv(files[i])
    # print(len(peptide_data))

    if len(peptide_data) > 0:
        count_non_empty += 1

print(count_non_empty)
print(1 - count_non_empty/len(files))