import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt 

def attention_filename_for_hla_and_peptide(hla, peptide):
    return f"results/attention\\{hla.replace(':', '_').replace('*', '_')}_{peptide}_attention.csv"

aa_array = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

contribs = np.zeros((20, 14))
contribs_pos = np.zeros((20, 14))
contribs_neg = np.zeros((20, 14))

pairs_data = pd.read_csv("1696_iedb_random_negative_pairs.csv")
# print(pairs_data)
for i in range(len(pairs_data)):
    pep = pairs_data['peptide'][i]
    hla = pairs_data['hla'][i]

    if os.path.isfile(f"./{attention_filename_for_hla_and_peptide(hla, pep)}"):
        # print(f"{attention_filename_for_hla_and_peptide(hla, pep)} exists")
        data = pd.read_csv(attention_filename_for_hla_and_peptide(hla, pep))

        aminoacids = [data.columns[i][0] for i in range(1, len(data.columns))]

        for i in range(1, len(data.columns)):
            # print("pos:", i, "aa:", aminoacids[i-1], "contribution:", data.iloc[-1, i])
            idx = aa_array.index(aminoacids[i-1])

            contribs[idx, i-1] += data.iloc[-1, i]

transphla_pairs = pd.read_csv("results/predict_results.csv")

# print(transphla_pairs)
# print(transphla_pairs['y_pred'])
for i in range(len(transphla_pairs)):
    pep = transphla_pairs['peptide'][i]
    hla = transphla_pairs['HLA'][i]
    pred = transphla_pairs['y_pred'][i]

    if pred == 1:
        # print("pred=1")
        if os.path.isfile(f"./{attention_filename_for_hla_and_peptide(hla, pep)}"):
            # print(f"{attention_filename_for_hla_and_peptide(hla, pep)} exists")
            data = pd.read_csv(attention_filename_for_hla_and_peptide(hla, pep))

            aminoacids = [data.columns[i][0] for i in range(1, len(data.columns))]

            for i in range(1, len(data.columns)):
                # print("pos:", i, "aa:", aminoacids[i-1], "contribution:", data.iloc[-1, i])
                idx = aa_array.index(aminoacids[i-1])

                contribs_pos[idx, i-1] += data.iloc[-1, i]

    else:
        # print("pred=0")
        if os.path.isfile(f"./{attention_filename_for_hla_and_peptide(hla, pep)}"):
            # print(f"{attention_filename_for_hla_and_peptide(hla, pep)} exists")
            data = pd.read_csv(attention_filename_for_hla_and_peptide(hla, pep))

            aminoacids = [data.columns[i][0] for i in range(1, len(data.columns))]

            for i in range(1, len(data.columns)):
                # print("pos:", i, "aa:", aminoacids[i-1], "contribution:", data.iloc[-1, i])
                idx = aa_array.index(aminoacids[i-1])

                contribs_neg[idx, i-1] += data.iloc[-1, i]


fig, axs = plt.subplots(1, 3, figsize=(15, 15))

# First subplot
axs[0].imshow(contribs)
axs[0].set_title("All samples")
axs[0].set_xticks(range(14))
axs[0].set_xticklabels(range(1, 15))
axs[0].set_yticks(range(len(aa_array)))
axs[0].set_yticklabels(aa_array)
fig.colorbar(axs[0].imshow(contribs), ax=axs[0])

# Second subplot
axs[1].imshow(contribs_pos[:-1])
axs[1].set_title("Samples TransPHLA classified \nas positive")
axs[1].set_xticks(range(14))
axs[1].set_xticklabels(range(1, 15))
axs[1].set_yticks(range(len(aa_array)))
axs[1].set_yticklabels(aa_array)
fig.colorbar(axs[1].imshow(contribs_pos[:-1]), ax=axs[1])

# Third subplot
axs[2].imshow(contribs_neg)
axs[2].set_title("Samples TransPHLA classified \nas negative")
axs[2].set_xticks(range(14))
axs[2].set_xticklabels(range(1, 15))
axs[2].set_yticks(range(len(aa_array)))
axs[2].set_yticklabels(aa_array)
fig.colorbar(axs[2].imshow(contribs_neg), ax=axs[2])

# Adjust layout to prevent overlap
plt.subplots_adjust(hspace=0.4)
# Display the plots
plt.show()
