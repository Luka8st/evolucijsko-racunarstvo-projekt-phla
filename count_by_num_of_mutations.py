import csv
import glob
import pandas as pd

files = glob.glob("results/mutation/*")
filenames = []

counters = [0,0,0,0]

for i in range(len(files)):
    peptide_data = pd.read_csv(files[i])
    num = peptide_data['mutation_AA_number']
    preds = peptide_data['y_pred'].tolist()[1:]

    # print(len(peptide_data))
    # print("preds:", preds)

    # print(preds.index(1) + 1, num[preds.index(1)+1])

    counters[num[preds.index(1)+1] - 1] += 1


print(counters)