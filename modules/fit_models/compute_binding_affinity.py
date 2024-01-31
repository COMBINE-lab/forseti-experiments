# In this script, we will use the trained MLP classifier to predict the binding affinity of the 30 mers of each transcript.

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import math
import pickle
from scipy import interpolate
from Bio import SeqIO
from sklearn.preprocessing import OneHotEncoder

parser = argparse.ArgumentParser(
    description="calculate the probability of a read arising from each of its references.")
# genome or txome. Currently need txome to find the polyA/T sites
# parser.add_argument("spliced_reads_genome_bam", type=str, help="genome-based BAM file containing the classified spliced reads")
# parser.add_argument("unspliced_reads_genome_bam", type=str, help="genome-based BAM file containing the classified spliced reads")
parser.add_argument("mlp_model_file", type=str,
                    help="trained MLP model for predicting the priming site in a pickle file")
parser.add_argument("spliceu_fasta_file", type=str,
                    help="spliceu reference sequence fasta file")
parser.add_argument("pa_min_len", type=int,help="polyA minimum length")
parser.add_argument("max_frag_len", type=int, default=500, help="polyA minimum length")
parser.add_argument("num_threads", type=int, default=2, help="number of threads")
parser.add_argument("outdir", type=str, help="output directory")

# args = parser.parse_args()
args = parser.parse_args([
    "/fs/nexus-projects/sc_frag_len/nextflow/test/mlp/mlp.pkl",
    "/fs/nexus-projects/sc_frag_len/nextflow/data/spliceu_ref/refdata-gex-GRCh38-2020-A/spliceu.fa",
    "7",
    "500",
    "16",
    "/fs/nexus-projects/sc_frag_len/nextflow/test/mlp"
])

# define a function that takes a read alignment, and the polyA intervals of the reference, and return the maximum probability of the read arising from the reference
#
# the read alignment is a pysam.AlignedSegment object
# the alignment tells us the position of the read on the reference
# and if the read is sense or antisense w.r.t. the reference
#
# The polyA and polyT are two dictionary of lists, where the keys are the transcript ids and the values are the polyA/T sites on the transcript

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:(i + ksize)]
        kmers.append(kmer)

    return kmers

encoder = OneHotEncoder(categories=[['A', 'C', 'G', 'T', 'N']] * 30, handle_unknown='ignore')

# read in the mlp model
with open(args.mlp_model_file, 'rb') as file_model:
    mlp = pickle.load(file_model)

# read in the spliceu fasta file
txome = {record.id: record.seq for record in SeqIO.parse(
    args.spliceu_fasta_file, "fasta")}

# calculate binding affinity
def calculate_binding_affinity(batch, mlp, encoder, batch_size=10_000):
    # we compute the binding affinity of each 30-mer of each transcript
    # if we have short transcript, we compute directly
    # if we have long transcript, we compute in batches
    out_list = []
    for tname1, seq1 in batch:
        if len(seq1) < 30:
            out_list.append((tname1, np.array([])))
        else:
        # elif len(seq) <= batch_size:
            # try:
            out_list.append((tname1, mlp.predict_proba(encoder.fit_transform([list(x) for x in build_kmers(seq1, 30)]))[:,1] * 0.75)) 
            # except IndexError as e:
            #     raise IndexError(f"tname: {tname} returned error: {e}")
            
        # else:
        #     # we compute in batches
        #     # we first get the number of batches
        #     num_batches = math.ceil(len(seq) / batch_size ) 
            
        #     # we use the first batch to initialize the array
        #     ba = mlp.predict_proba(encoder.fit_transform([list(x) for x in build_kmers(seq[0:batch_size], 30)]))[:,1] * 0.75

        #     for i in range(1, num_batches):
                
        #         batch_probs = mlp.predict_proba(encoder.fit_transform([list(x) for x in build_kmers(seq[(i*batch_size-30+1):min((i+1)*batch_size, len(seq))], 30)]))[:,1] * 0.75
        #         ba = np.hstack((ba,batch_probs))
        #     out_list.append((tname, ba))
    return out_list
    
# Use joblib's Parallel to parallelize predictions
# reads = spliced_reads_txome_bam.fetch(until_eof=True)

batch_size = math.ceil(len(txome) / args.num_threads)
shuf_keys = list(txome.keys())
np.random.shuffle(shuf_keys)

batches = []
for batch_start in range(0, len(shuf_keys), batch_size):
    batches.append(list(zip(shuf_keys[batch_start:min(batch_start+batch_size, len(shuf_keys))], [txome[k] for k in shuf_keys[batch_start:min(batch_start+batch_size, len(shuf_keys))]])))

# joblib version
from joblib import Parallel, delayed

# Use joblib's Parallel to parallelize predictions
tx_binding_affinity_batched = Parallel(n_jobs=args.num_threads)(
    delayed(calculate_binding_affinity)(batch, mlp, encoder) for batch in batches
)

tx_binding_affinity_dict = {}
for batch in tx_binding_affinity_batched:
    for (tname, ba) in batch:
        tx_binding_affinity_dict[tname] = ba


model_pkl_path = os.path.join(args.outdir, 'tx_binding_affinity_dict.pkl')
with open(model_pkl_path, 'wb') as file_model:
    pickle.dump(tx_binding_affinity_dict, file_model)

