import csv

true_negs = []
hlas = []
hla_sequences = []
peptides = []
y_preds = []
y_probs = []

with open('results/predict_results.csv', mode ='r') as file: 
    csvFile = csv.DictReader(file)
    for lines in csvFile:
	    # print(lines['y_pred']=='0')
        if lines['y_pred'] == '0':
            print(lines)
            hlas.append(lines['HLA'])
            hla_sequences.append(lines['HLA_sequence'])
            peptides.append(lines['peptide'])
            y_preds.append(lines['y_pred'])
            y_probs.append(lines['y_prob'])
            true_negs.append((lines['HLA'], lines['HLA_sequence'], lines['peptide'], lines['y_pred'], lines['y_prob']))
      
csv_file = 'true_negatives.csv'

# Pisanje podataka u CSV datoteku
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Pisanje zaglavlja (opciono)
    writer.writerow(['HLA', 'HLA_sequence', 'peptide', 'y_pred', 'y_prob'])
    # Pisanje podataka
    writer.writerows(true_negs)