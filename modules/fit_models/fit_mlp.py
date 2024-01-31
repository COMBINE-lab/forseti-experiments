import os
import argparse
import pandas as pd
import numpy as np
from scipy import interpolate
from collections import Counter
import scipy.sparse as sp
import math
import random
import pickle
import regex
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score, mean_squared_error

"""
This script trains a MLP classifier to predict the binding affinity of a 30 bases sequence w.r.t. the obtained priming window.

The training data is obtained from the upstream 30 bases of the alignments of R1 in our paired-end read alignment results from the genome. The training data has the same orientation as its read2, so, it should contain a polyA site.
As reads might come from polyA tails, which are not included in the genome, we would expect that some of the input data come from intergenic regions and do not correspond to a polyA site. This is saying our labels are noisy, so we would not expect to obtain perfect accuracy.

We train a MLP classifier by reading the training data from each SRR and call `partial_fit` to update the model. 

"""

parser = argparse.ArgumentParser(
    description="Count and summarize non-zero length reads.")
parser.add_argument("parent_dir", type=str,
                    help="the parent dir having all fragment length outputs")
parser.add_argument("PE_sheet", type=str,
                    help="path to the spreadsheet with GSR,SRR list")
parser.add_argument("random_seed", type=int, default=1,
                    help="random seed for selecting training data")
parser.add_argument("snr_len", type=int, default=6,
                    help="random seed for selecting training data")
parser.add_argument("outdir", type=str,
                    help="where to save the output")
args = parser.parse_args()
# args = parser.parse_args(["/fs/nexus-projects/sc_frag_len/nextflow/workflow_output",
#                           "/fs/nexus-projects/sc_frag_len/nextflow/input_files/sample_url_sheet.csv",
#                           "1",
#                           "6",
#                           "/fs/nexus-projects/sc_frag_len/nextflow/test/mlp"])
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

encoder = OneHotEncoder(categories=[['A', 'C', 'G', 'T', 'N']] * 30, handle_unknown='ignore')

parent_dir = args.parent_dir
# 1----- Get the PE datasets spreadsheet
PE_sheet = pd.read_csv(args.PE_sheet)

# 2----- loop through GSE(s), combine all tlen from its SRR
# check if we have all datasets processed
missing_files = []

for (GSE, group_gse_lst) in PE_sheet.groupby('GSE'):
    SRR_lst = group_gse_lst['SRR']
    for SRR in SRR_lst:
        polya_path = os.path.join(
            parent_dir, "process_data", "frag_len_dist", GSE, SRR, "priming_site_seqs", "polya_seq.txt")

        if os.path.exists(polya_path):
            check_file = os.path.getsize(polya_path)
            if (check_file == 0):
                missing_files.append(f"{GSE}-{SRR}")
                error_occur = True
        else:
            missing_files.append(f"{GSE}-{SRR}")
            error_occur = True

if missing_files:
    raise ValueError(f"Please re-run the previous step, the output of following dataset(s) is either missing or empty: {missing_files}")

polya_mlp = MLPClassifier(hidden_layer_sizes=(100,), max_iter=500, alpha=1e-4,
                    solver='adam', verbose=10, random_state=args.random_seed, shuffle=True,
                    learning_rate_init=.001)

held_out_polya = sp.csr_matrix((0, 150))
held_out_polya_bg = sp.csr_matrix((0, 150)) 

random.seed(args.random_seed)
train_list = PE_sheet.loc[PE_sheet['GSE'].isin(random.sample(PE_sheet['GSE'].unique().tolist(), 8))]
num_srr = {GSE: group_gse_lst.shape[0] for (GSE, group_gse_lst) in train_list.groupby('GSE')}
# initialize the sparse matrix 
training_batch = sp.csr_matrix((0, 150))
held_out_polya_bg = sp.csr_matrix((0, 150)) 

# we group by GSE, then loop through SRRs to get the polya seq
for i, (GSE, group_gse_lst) in enumerate(train_list.groupby('GSE')):
    SRR_lst = group_gse_lst['SRR']
    num_held_out = 1000//len(SRR_lst)
    
    for SRR in SRR_lst:
        # The 30 bases are defined according to R2. So, they should always be polyA

        polya_path = os.path.join(
            parent_dir, "process_data", "frag_len_dist", GSE, SRR, "priming_site_seqs", "polya_seq.txt")
        bg_path = os.path.join(
            parent_dir, "process_data", "frag_len_dist", GSE, SRR, "priming_site_seqs", "polya_bg_seq.txt")

        # get polya and bg seq
        # we want to make sure the polya seq has a polyA six mer with at most one mismatch
        polya_batch = encoder.fit_transform([list(x.strip()) for x in open(polya_path).readlines() if len(x.strip()) == 30 and regex.search("A(" + 'A' * (args.snr_len-2)+ "){s<=1}A", x.strip())])
        # polya_batch = encoder.fit_transform([list(x.strip()) for x in open(polya_path).readlines() if len(x.strip()) == 30])
        
        # if we do not have enough data, then skip this SRR
        if polya_batch.shape[0] < num_held_out * 2:
            print("Not enough training examples for", GSE, SRR)
            continue
        
        bg_batch = encoder.fit_transform([list(x.strip()) for x in open(bg_path).readlines() if len(x.strip()) == 30])
        
        # if we do not have enough data, then skip this SRR
        if bg_batch.shape[0] < num_held_out * 2:
            print("Not enough background examples for", GSE, SRR)
            continue
        
        # num_held_out = min(math.ceil(polya_batch.shape[0] * 0.1), 20)
        num_train_fg = polya_batch.shape[0]
        num_train_bg = min(num_train_fg, bg_batch.shape[0])

        polya_batch = polya_batch[np.random.randint(0, polya_batch.shape[0], size = num_train_fg), : ]
        bg_batch = bg_batch[np.random.randint(0, bg_batch.shape[0], size = num_train_bg), : ]

        # append held out data
        polya_held_out = polya_batch[:num_held_out,]
        bg_held_out = bg_batch[:num_held_out,]
        
        held_out_polya = sp.vstack([held_out_polya, polya_held_out])
        held_out_polya_bg = sp.vstack([held_out_polya_bg, bg_held_out])
        
        fg_shuf_id = np.arange((num_train_fg+num_train_bg-num_held_out*2))
        np.random.shuffle(fg_shuf_id)

        # build training data
        train_polya = sp.vstack([
            polya_batch[num_held_out:],
            bg_batch[num_held_out:]
        ])[fg_shuf_id,:]
        
        del polya_batch
        del bg_batch
        
        label_polya = np.hstack([
            np.ones(num_train_fg - num_held_out), 
            np.zeros(num_train_bg - num_held_out)
        ])[fg_shuf_id]
        
        # fit the model using the data from this SRR
        polya_mlp.partial_fit(
            train_polya,
            label_polya,
            classes=[0,1]
        )

        # get the accuracy on the held out data
        print("Accuracy:", 
            accuracy_score(
                np.hstack(
                    [np.ones(num_held_out), 
                    np.zeros(num_held_out)
                    ]
                ), 
                polya_mlp.predict(
                    sp.vstack([
                        polya_held_out, 
                        bg_held_out]
                    )
                )
            )
        )


# finally, we want to get the overall accuracy on all held out data
predictions = polya_mlp.predict(sp.vstack([held_out_polya, held_out_polya_bg]) )
print("Accuracy:", accuracy_score(np.hstack([np.ones(held_out_polya.shape[0]), np.zeros(held_out_polya_bg.shape[0])]), predictions))

# we a;sp also want to test the accuracy on test data
# we read in the test data
test_list = PE_sheet.loc[~PE_sheet['GSE'].isin(random.sample(PE_sheet['GSE'].unique().tolist(), 8))]
test_acc_list = []
# we group by GSE, then loop through SRRs to get the polya seq
for i, (GSE, group_gse_lst) in enumerate(test_list.groupby('GSE')):
    SRR_lst = group_gse_lst['SRR']
    for SRR in SRR_lst:
        # The 30 bases are defined according to R2. So, they should always be polyA

        polya_path = os.path.join(
            parent_dir, "process_data", "frag_len_dist", GSE, SRR, "priming_site_seqs", "polya_seq.txt")
        bg_path = os.path.join(
            parent_dir, "process_data", "frag_len_dist", GSE, SRR, "priming_site_seqs", "polya_bg_seq.txt")

        # get polya and bg seq
        # we want to make sure the polya seq has a polyA six mer with at most one mismatch
        polya_batch = encoder.fit_transform([list(x.strip()) for x in open(polya_path).readlines() if len(x.strip()) == 30 and regex.search("A(" + 'A' * (args.snr_len-2)+ "){s<=1}A", x.strip())])
        # polya_batch = encoder.fit_transform([list(x.strip()) for x in open(polya_path).readlines() if len(x.strip()) == 30])
        
        # if we do not have enough data, then skip this SRR
        if polya_batch.shape[0] < num_held_out * 2:
            print("Not enough training examples for", GSE, SRR)
            continue
        
        bg_batch = encoder.fit_transform([list(x.strip()) for x in open(bg_path).readlines() if len(x.strip()) == 30])
        
        # if we do not have enough data, then skip this SRR
        if bg_batch.shape[0] < num_held_out * 2:
            print("Not enough background examples for", GSE, SRR)
            continue
        
        # build training data
        train_polya = sp.vstack([
            polya_batch[random.sample(range(polya_batch.shape[0]), min(polya_batch.shape[0], bg_batch.shape[0])), : ],
            bg_batch[random.sample(range(bg_batch.shape[0]), min(polya_batch.shape[0], bg_batch.shape[0])), : ]
        ])
        label_polya = np.hstack([
            np.ones(min(polya_batch.shape[0], bg_batch.shape[0])), 
            np.zeros(min(polya_batch.shape[0], bg_batch.shape[0]))
        ])
        
        
        del polya_batch
        del bg_batch
        
        # get the accuracy on the held out data
        acc_scores = accuracy_score(
            label_polya, 
            polya_mlp.predict(
                train_polya
            )
        )
        test_acc_list.append(acc_scores)
        
        print("Accuracy:", acc_scores)



# priming/binding affinity 
# later use

# mlp.predict_proba(sp.vstack([held_out_polya, held_out_polya_bg]) )
# mlp.predict_proba(held_out_polya)[:,1]
# mlp.predict_proba(held_out_polya_bg)[:,1]

"""
    ------------                               |
        AAAAAAAA
        ------------                           |

"""

# For classification

model_pkl_path = os.path.join(outdir, 'mlp_pa6with1mm.pkl')
with open(model_pkl_path, 'wb') as file_model:
    pickle.dump(polya_mlp, file_model)

