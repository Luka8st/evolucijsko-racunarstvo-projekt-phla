import glob
import pandas as pd
import shlex, subprocess
import csv

def prediciton(peptide, allele, length):
    command = 'curl --data "method=recommended&sequence_text=' + peptide + '&allele=' + allele + '&length=' + str(length) + '" http://tools-cluster-interface.iedb.org/tools_api/mhci/'
    args = shlex.split(command)
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    output_list = output[0].decode('utf8').split('\t')
    
    print(output_list, len(output_list))
    # Provera dužine liste pre pristupa elementu
    if len(output_list) == 19:
        consesus_percentile = output_list[-1]
        score = output_list[-2]
    else:
        consesus_percentile = '' # ili neka druga podrazumevana vrednost
        score = ''
    

    return consesus_percentile, score

def peptide_to_string(peptide, index):
    return f"%3Epeptide{index+1}%0A{peptide}"

files = glob.glob("results/mutation/*")
filenames = []

for file in files:
    # print("file:",file.replace("results/mutation\\", "").replace("_mutation.csv", ""))
    filenames.append(str(file.replace("results/mutation\\", "").replace("_mutation.csv", "")))
    hla, peptide = file.replace("results/mutation\\", "").replace("_mutation.csv", "").rsplit('_', 1)
    hla = hla.replace('_', '*', 1).replace('_', ':')
    # print(hla, peptide)

print(len(filenames))

files = glob.glob("results/mutation/*")


for i in range(len(files)):
    peptide_data = pd.read_csv(files[i])
    alleles = peptide_data['HLA'].tolist()
    num_of_mutations = peptide_data['mutation_AA_number'].tolist()
    peptides = peptide_data['mutation_peptide'].tolist()
    
    mutated_peptides = []
    percentiles = []
    numMuts= []
    scores = []
    generated = []


    for j in range(len(alleles)):
        print("j=", j)
        peptide = peptides[j]
        allele = alleles[j]
        num = num_of_mutations[j]
        length = len(peptide)
        consesus_percentile, score = prediciton(peptide, allele, length)

        if score != '' and float(score) >= 0.5:
            print("found!!!!")
            mutated_peptides.append(peptide)
            percentiles.append(consesus_percentile)
            scores.append(score)
            numMuts.append(num_of_mutations)
            generated.append((peptide, num, consesus_percentile, score))

    reduced_name = filenames[i]
    csv_file = f'results/successful_mutation/{reduced_name}.csv'

    # Pisanje podataka u CSV datoteku
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Pisanje zaglavlja (opciono)
        writer.writerow(['mutated_peptide', 'number_of_mutations', 'percentile', 'score'])
        # Pisanje podataka
        writer.writerows(generated)