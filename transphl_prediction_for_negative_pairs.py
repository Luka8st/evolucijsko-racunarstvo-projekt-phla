import csv
import pandas as pd
from model import *
from mutation import *

threshold = 0.5
output_dir = './results/'
mut_savepath = output_dir + '/mutation/'
attn_savepath = output_dir + '/attention/'
fig_savepath = output_dir + '/figures/'
batch_size = 1024

def find_hla_index(hlas, hla):
    for i in range(len(hlas)):
        if hlas[i] == hla:
            return i
    return -1

available_lengths_for_hlas_dict = {
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
peptides = []

with open('1696_iedb_random_negative_pairs.csv', mode='r') as file:
    csvFile = csv.DictReader(file)
    for lines in csvFile:
        # hlas.append(lines['hla'].replace('*', '_').replace(':', '_'))
        hlas.append(lines['HLA'])
        hla_sequences.append(lines['HLA_sequence'])
        peptides.append(lines['peptide'])

# print(hlas)
# print(hlas, hla_sequences, peptides)

predict_data = pd.DataFrame([hlas, hla_sequences, peptides], index = ['HLA', 'HLA_sequence', 'peptide']).T

predict_data, predict_pep_inputs, predict_hla_inputs, predict_loader = read_predict_data(predict_data, batch_size)

use_cuda = True
device = torch.device("cuda")
print('device:', device)

model_file = 'pHLAIformer.pkl'

model_eval = Transformer().to(device)
# u state dictionaryju su obično pohranjeni weights i biases
model_eval.load_state_dict(torch.load(model_file, map_location='cpu'), strict = True)

state_dict_keys = model_eval.state_dict().keys()

model_eval.eval()
y_pred, y_prob, attns = eval_step(model_eval, predict_loader, threshold, True)

predict_data['y_pred'], predict_data['y_prob'] = y_pred, y_prob
predict_data = predict_data.round({'y_prob': 4})

predict_data.to_csv(output_dir + '/predict_results.csv', index = False)






"""iterator=1
for hla, pep in zip(predict_data.HLA, predict_data.peptide):
    print(iterator)
    iterator += 1

    print(hla, pep)"""
    #pHLA_attns_draw_save(predict_data, attns, hla, pep, attn_savepath, fig_savepath)


# print("attns_new", attns)
attns_cp = [tensor.cpu().numpy().copy() for tensor in attns]
predict_data_cp = predict_data.copy()
for idx in range(predict_data.shape[0]):
     print("idx=", idx)
    #  if idx < 72:
    #       continue
     #print("idx =", idx, len(predict_data))
     #attns = attns_cp.copy()

     peptide = predict_data.iloc[idx].peptide
     hla = predict_data.iloc[idx].HLA


     pHLA_attns_draw_save(predict_data, attns, hla, peptide, attn_savepath, fig_savepath)



    #predict_data = predict_data_cp
     mut_peptides_df = pHLA_mutation_peptides(predict_data, attns, hla = hla, peptide = peptide)

     print("mut_peptides_df:", mut_peptides_df)

     index = find_hla_index(hlas, hla)
     hla_seq = hla_sequences[index]

    # for i in range(len(mut_peptides_df)):
    #     mut_peptides_df[i, 'HLA_sequence'] = hla_seq
     mut_peptides_df['HLA_sequence'] = hla_seq

     mut_data, _, _, mut_loader = read_predict_data(mut_peptides_df, batch_size)
    
    # print("mut_data:", mut_data)

     model_eval = Transformer().to(device)
     model_eval.load_state_dict(torch.load(model_file, map_location='cpu'), strict = True)

     model_eval.eval()
     y_pred, y_prob, attns = eval_step(model_eval, mut_loader, threshold, True)
    
     mut_data['y_pred'], mut_data['y_prob'] = y_pred, y_prob
     mut_data = mut_data.round({'y_prob': 4})

     print("mut_data:", mut_data)
     print("path:", mut_savepath + '{}_{}_mutation.csv'.format(hla.replace('*', '_').replace(':', '_'), peptide))
     mut_data.to_csv(mut_savepath + '{}_{}_mutation.csv'.format(hla.replace('*', '_').replace(':', '_'), peptide), index = False)

     print("------------------")
     print('********** {} | {} → # Mutation peptides = {}'.format(hla, peptide, mut_data.shape[0]-1))
     print("=======================")  
    
     mut_peptides_IEDBfmt = ' '.join(mut_data.mutation_peptide)
     print('If you want to use IEDB tools to predict IC50, please use these format: \n {}'.format(mut_peptides_IEDBfmt))

