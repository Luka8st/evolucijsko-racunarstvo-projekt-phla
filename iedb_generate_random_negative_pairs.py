import requests
import json
import pandas as pd
from io import StringIO
import numpy as np
import csv

base_uri='https://query-api.iedb.org'

# funciton to print the CURL command given a request
def print_curl_cmd(req):
    url = req.url
    print("curl -X 'GET' '" + url + "'")

def generate_10_random_negative_peptides_for_a_hla(hla, hla_sequence, peptide_length):
    # hla = hlas[index]
    # hla_sequence = hla_sequences[index]

    pairs = []
    generated = 0
    max_iterations = 30

    # result = requests.get("https://query-api.iedb.org/mhc_search?mhc_restriction=eq." + hla + "&mhc_class=eq.I&qualitative_measure=eq.Negative&host_organism_name=like.%human%&linear_sequence_length=gte.8&linear_sequence_length=lte.14&linear_sequence_length=eq." + str(peptide_length) + "&select=elution_id,elution_iri,structure_description,mhc_restriction")
    result = requests.get("https://query-api.iedb.org/mhc_search?mhc_restriction=eq." + hla + "&mhc_class=eq.I&qualitative_measure=eq.Negative&host_organism_name=like.%human%&linear_sequence_length=eq." + str(peptide_length) + "&select=elution_id,elution_iri,structure_description,mhc_restriction")
    print_curl_cmd(result)
    df = pd.json_normalize(result.json())
    print("len(df)", len(df))

    elution_ids = []
        
    it = 0
    while generated < min(10, len(df)):
        if it > max_iterations:
            break
        it += 1
            
        r = np.random.randint(len(df))
        peptide = df.loc[r, 'structure_description']
        elution_id = df.loc[r, 'elution_id']

        if elution_id in elution_ids or len(set(peptide).difference(set('ARNDCQEGHILKMFPSTWYV'))) != 0:
            continue

        pairs.append((peptide, hla, hla_sequence))
        elution_ids.append(elution_id)

        generated += 1

        # while True:
        #     print(i, ". iteration")
        #     u = np.random.randint(len(hlas))
        #     hla = hlas[u]
        #     hla_sequence = hla_sequences[u]

        #     result = requests.get("https://query-api.iedb.org/mhc_search?mhc_restriction=eq." + hla + "&mhc_class=eq.I&qualitative_measure=eq.Negative&host_organism_name=like.%human%&linear_sequence_length=gte.8&linear_sequence_length=lte.14&select=elution_id,elution_iri,structure_description,mhc_restriction")
        #     print_curl_cmd(result)
        #     df = pd.json_normalize(result.json())

        #     if len(df) > 0:
        #         break

    return pairs

def find_hla_index(hlas, hla):
    for i in range(len(hlas)):
        if hlas[i] == hla:
            return i
    return -1

hla_peptide_length_dict = {
    'HLA-A*01:01': (8, 14),
    'HLA-A*02:01': (8, 14),
    'HLA-A*02:02': (9, 10),
    'HLA-A*02:03': (9, 11),
    'HLA-A*02:04': (9, 11),
    'HLA-A*02:05': (9, 11),
    'HLA-A*02:06': (9, 10),
    'HLA-A*02:07': (9, 11),
    'HLA-A*02:11': (9, 9),
    'HLA-A*02:12': (9, 9),
    'HLA-A*02:16': (9, 9),
    'HLA-A*02:17': (9, 10),
    'HLA-A*02:19': (9, 9),
    'HLA-A*02:20': (9, 9),
    'HLA-A*02:50': (9, 9),
    'HLA-A*03:01': (8, 13),
    'HLA-A*11:01': (8, 13),
    'HLA-A*23:01': (9, 11),
    'HLA-A*24:02': (8, 14),
    'HLA-A*24:03': (9, 9),
    'HLA-A*24:06': (9, 11),
    'HLA-A*24:13': (9, 9),
    'HLA-A*25:01': (9, 9),
    'HLA-A*26:01': (9, 10),
    'HLA-A*26:02': (9, 9),
    'HLA-A*26:03': (9, 9),
    'HLA-A*29:02': (8, 13),
    'HLA-A*30:01': (9, 10),
    'HLA-A*30:02': (9, 10),
    'HLA-A*31:01': (9, 13),
    'HLA-A*32:01': (9, 11),
    'HLA-A*32:07': (9, 9),
    'HLA-A*32:15': (9, 9),
    'HLA-A*33:01': (9, 10),
    'HLA-A*66:01': (9, 9),
    'HLA-A*68:01': (9, 12),
    'HLA-A*68:02': (9, 14),
    'HLA-A*68:23': (9, 9),
    'HLA-A*69:01': (9, 10),
    'HLA-A*80:01': (9, 9),
    'HLA-B*07:02': (8, 14),
    'HLA-B*08:01': (8, 14),
    'HLA-B*13:02': (8, 10),
    'HLA-B*14:01': (9, 9),
    'HLA-B*14:02': (8, 10),
    'HLA-B*15:01': (8, 14),
    'HLA-B*15:02': (9, 9),
    'HLA-B*15:03': (9, 9),
    'HLA-B*15:09': (9, 9),
    'HLA-B*15:11': (9, 9),
    'HLA-B*15:17': (9, 9),
    'HLA-B*15:18': (9, 9),
    'HLA-B*15:42': (9, 9),
    'HLA-B*18:01': (8, 14),
    'HLA-B*18:03': (8, 9),
    'HLA-B*27:01': (9, 13),
    'HLA-B*27:02': (9, 13),
    'HLA-B*27:03': (9, 12),
    'HLA-B*27:04': (9, 11),
    'HLA-B*27:05': (8, 14),
    'HLA-B*27:06': (9, 11),
    'HLA-B*27:07': (9, 12),
    'HLA-B*27:08': (9, 13),
    'HLA-B*27:09': (8, 14),
    'HLA-B*27:20': (9, 9),
    'HLA-B*35:01': (8, 14),
    'HLA-B*35:03': (9, 11),
    'HLA-B*35:08': (9, 11),
    'HLA-B*37:01': (8, 11),
    'HLA-B*38:01': (9, 9),
    'HLA-B*39:01': (8, 11),
    'HLA-B*39:06': (9, 9),
    'HLA-B*39:24': (8, 9),
    'HLA-B*40:01': (8, 12),
    'HLA-B*40:02': (8, 12),
    'HLA-B*41:01': (9, 10),
    'HLA-B*44:02': (8, 13),
    'HLA-B*44:03': (8, 12),
    'HLA-B*44:27': (10, 10),
    'HLA-B*45:01': (9, 11),
    'HLA-B*45:06': (9, 9),
    'HLA-B*46:01': (8, 11),
    'HLA-B*48:01': (9, 9),
    'HLA-B*49:01': (8, 11),
    'HLA-B*50:01': (9, 10),
    'HLA-B*51:01': (8, 13),
    'HLA-B*51:08': (8, 9),
    'HLA-B*52:01': (8, 9),
    'HLA-B*53:01': (9, 10),
    'HLA-B*54:01': (8, 11),
    'HLA-B*56:01': (9, 11),
    'HLA-B*57:01': (8, 14),
    'HLA-B*57:03': (8, 13),
    'HLA-B*58:01': (8, 13),
    'HLA-B*73:01': (9, 9),
    'HLA-B*83:01': (9, 9),
    'HLA-C*01:02': (8, 12),
    'HLA-C*02:02': (8, 11),
    'HLA-C*03:03': (8, 11),
    'HLA-C*03:04': (8, 11),
    'HLA-C*04:01': (8, 14),
    'HLA-C*05:01': (8, 14),
    'HLA-C*06:02': (8, 14),
    'HLA-C*07:01': (8, 12),
    'HLA-C*07:02': (8, 11),
    'HLA-C*07:04': (8, 10),
    'HLA-C*08:02': (8, 11),
    'HLA-C*12:03': (8, 9),
    'HLA-C*14:02': (8, 10),
    'HLA-C*15:02': (8, 9),
    'HLA-C*16:01': (8, 11),
    'HLA-C*17:01': (8, 9)
}

hlas = []
hla_sequences = []
with open('common_hla_sequence.csv', mode ='r') as file: 
	csvFile = csv.DictReader(file)
	for lines in csvFile:
			hlas.append(lines['HLA'])
			hla_sequences.append(lines['HLA_sequence'])

generated = []
for hla in hla_peptide_length_dict.keys():
    index = find_hla_index(hlas, hla)
    hla_seq = hla_sequences[index]
    # print(index)
    for length in range(hla_peptide_length_dict[hla][0], hla_peptide_length_dict[hla][1]+1):
        
        res = generate_10_random_negative_peptides_for_a_hla(hla, hla_seq, length)
        print(f"for HLA {hla}, peptide_length {length}, {len(res)} pairs were generated")
        generated.extend(res)
        print(f"total generated {len(generated)}")

csv_file = str(len(generated)) + '_iedb_random_negative_pairs.csv'

# Pisanje podataka u CSV datoteku
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Pisanje zaglavlja (opciono)
    writer.writerow(['peptide', 'hla', 'hla_sequence'])
    # Pisanje podataka
    writer.writerows(generated)