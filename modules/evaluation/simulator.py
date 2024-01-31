import os
import argparse
import numpy as np
import pandas as pd
import pyranges as pr
import random
import pickle
import math
from scipy import interpolate
import pyroe
import time
from joblib import Parallel, delayed
from datetime import datetime
import gzip
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord

from iss.error_models import kde

parser = argparse.ArgumentParser()
parser.add_argument('--spliceu_ref_human_dir', type=str,
                    help='human spliceu reference directory path')
parser.add_argument('--simpleaf_outdir', type=str,
                    help='simpleaf quant output dir')
parser.add_argument('--frag_len_model_path', type=str,
                    help='fragment length model path')
parser.add_argument('--err_model_path', type=str, help='error model path')
parser.add_argument('--out_dir', type=str, help='output dir')
# parser.add_argument('--spliceu_txome_path', type=str,
#                     help='spliceu txome fasta file')
# parser.add_argument('--t2g_path', type=str, help='t2g path')
# parser.add_argument('--gtf_path', type=str, help='gtf path')
# parser.add_argument('--human_polyA_path', type=str, help='human polyA path')
# parser.add_argument('--spliceu_fa_path', type=str,
#                     help='spliceu reference sequence path')

args = parser.parse_args()

spliceu_ref_human_dir = args.spliceu_ref_human_dir
simpleaf_outdir = args.simpleaf_outdir
frag_len_model_path = args.frag_len_model_path
err_model_path = args.err_model_path
out_dir = args.out_dir
# get the path of spliceu reference files from the spliceu_ref_human_dir
spliceu_txome_fa_path = os.path.join(spliceu_ref_human_dir, "spliceu.fa")
t2g_path = os.path.join(spliceu_ref_human_dir, "t2g_3col.tsv")
gtf_path = os.path.join(spliceu_ref_human_dir, "spliceu.gtf")
human_polyA_path = os.path.join(spliceu_ref_human_dir, "tx_polya.bed")
# spliceu_txome_path = args.spliceu_txome
# t2g_path = args.t2g_path
# gtf_path = args.gtf_path
# human_polyA_path = args.human_polyA_path
# spliceu_fa_path = args.spliceu_fa_path

spliceu_txome_fa_path = "/fs/nexus-projects/sc_frag_len/nextflow/data/spliceu_ref/refdata-gex-GRCh38-2020-A/spliceu.fa"
err_model_path = '/fs/nexus-projects/sc_frag_len/nextflow/simulation/iss_project/err_model.npz'
out_dir = "/fs/nexus-projects/sc_frag_len/nextflow/simulation/91bp_simu/test"
simpleaf_outdir = '/fs/nexus-projects/sc_frag_len/nextflow/simulation/process_real_data/simpleaf_workdir/workflow_output/simpleaf_quant/af_quant'
t2g_path = "/fs/nexus-projects/sc_frag_len/nextflow/data/spliceu_ref/refdata-gex-GRCh38-2020-A/t2g_3col.tsv"
gtf_path = '/fs/nexus-projects/sc_frag_len/nextflow/data/spliceu_ref/refdata-gex-GRCh38-2020-A/spliceu.gtf'
human_polyA_path = '/fs/nexus-projects/sc_frag_len/nextflow/data/spliceu_ref/refdata-gex-GRCh38-2020-A/tx_polya.bed'
frag_len_model_path = "/fs/nexus-projects/sc_frag_len/nextflow/fit_model/Jan8_newout/spline_model.pkl"


os.makedirs(out_dir, exist_ok=True)

# {S_U_count}------ We need to calc the number of S/U reads for each gene
fry_quant = pyroe.load_fry(frydir=simpleaf_outdir, output_format='raw')

col_sum_S = fry_quant.layers['spliced'].sum(axis=0).tolist()[0]
col_sum_U = fry_quant.layers['unspliced'].sum(axis=0).tolist()[0]
gene_list = fry_quant.var.index.values
del fry_quant
S_U_count = {gene: (int(s), int(u))
             for gene, s, u in zip(gene_list, col_sum_S, col_sum_U) if int(s) + int(u) != 0}

t2g_df = pd.read_csv(t2g_path, sep="\t", header=None,
                     names=["txid", "gid", "splicing_status"])
t2g = t2g_df.set_index("txid").to_dict()["gid"]
t2s = t2g_df.set_index("txid").to_dict()["splicing_status"]
single_exon_tx = set([tx for tx in t2g.keys() if not tx.endswith(
    '-U') and tx + '-U' not in t2g.keys()])

# {Stx2ee,Utx2epos}------- We collect exon info, in order to check if the sampled r2 is valid
# for each tx, we get the ee junctions site by cumsum the sorted exon length
# we get the exon sites for each U-tx by substracting the start site of 1st exon
start_time = time.time()
gtf_df = pr.read_gtf(gtf_path)
gtf_df['exon_rank'].df.astype(int)
exon_info = gtf_df[(gtf_df.Feature == 'exon')][[
    'transcript_id', 'gene_id', 'exon_rank']]
exon_info = exon_info[~exon_info.transcript_id.str.endswith("-U")].df
exon_info = exon_info.astype({"exon_rank": int})
Stx2ee = {}
Utx2epos = {}
for tx, tx_exon_info in exon_info.groupby('transcript_id', observed=True):
    # we sort the exon by exon_rank
    # if on negative strand, the exon with the largest genomic pos is the 1st exon
    tx_exon_info = tx_exon_info.sort_values(by='exon_rank')
    ee_list = (tx_exon_info['End'] - tx_exon_info['Start']).cumsum()
    Stx2ee[tx] = list(ee_list)
    if tx not in single_exon_tx:
        # save the exon start site(dereived from S-tx) for each U-tx
        # from absolute geno pos -> relative tx pos
        if tx_exon_info['Strand'].astype(str).tolist()[0] == '+':
            start_idx = tx_exon_info['Start'].iloc[0]
            Utx2epos[tx+'-U'] = list(zip(list(tx_exon_info['Start'] -
                                              start_idx), list(tx_exon_info['End']-start_idx)))
        else:
            start_idx = tx_exon_info['End'].iloc[0]
            Utx2epos[tx+'-U'] = list(zip(list(start_idx-tx_exon_info['End']),
                                     list(start_idx - tx_exon_info['Start'])))

# 152 sec
end_time = time.time()
print(
    f"--- calc Stx2ee,Utx2epos, Used {end_time - start_time} seconds ---")

start_time = time.time()
# {t2pa} ------- Get the polyA site for each tx

# file---> tx_id  polyA_start  polyA_end
t2pa = {}
with open(human_polyA_path, 'r') as file:
    for line in file:
        row = line.strip().split('\t')
        t2pa.setdefault(row[0], []).append(int(row[1]))
# add the terminal polyA to the polyA list
for tx in set(t2pa.keys()) & set(Stx2ee.keys()):
    t2pa[tx].append(Stx2ee[tx][-1])
for tx in set(t2pa.keys()) & set(Utx2epos.keys()):
    t2pa[tx].append(Utx2epos[tx][-1][1])
# Used: 20 sec
end_time = time.time()
print(
    f"--- calc t2pa, Used {end_time - start_time} seconds ---")
# ------- We only sample from genes that count > 0
non_zero_gene_list = set(S_U_count.keys())
zero_gene_list = set(gene_list)-set(non_zero_gene_list)

# {gene_tx_list}------- For each gene, we sampling 1 tx from tx-S list
start_time = time.time()
valid_tx_list = set()
# notes: only valid gene nad valid tx in gene)tx_list
gene_tx_list = {}
i = 0
with open(t2g_path, 'r') as file:
    for line in file:
        i += 1
        row = line.strip().split('\t')
        txid, gid, splicing_status = row[0], row[1], row[2]
        # discard tx that have no counts
        if gid in zero_gene_list:
            continue
        # discard tx that have no polyA site
        if txid not in t2pa.keys():
            continue
        # all tx = len(tx_list)= len(tx2ee) = 392673
        # vs partial: len(t2pa) = 261685
        valid_tx_list.add(txid)
        if splicing_status == 'U':
            gene_tx_list.setdefault(gid, [[], []])[0].append(txid)
        else:
            gene_tx_list.setdefault(gid, [[], []])[1].append(txid)

# 60 sec now 1.2 seconds
end_time = time.time()
print(f"--- calc gene_tx_list, Used {end_time - start_time} seconds ---")

# ------- Prepare fragment length sampling
with open(frag_len_model_path, 'rb') as file:
    spline = pickle.load(file)
x_test = np.arange(0, 1000)
frag_len_weight = interpolate.splev(x_test, spline)

start_time = time.time()
count = 0
# for each gene, we sampling: 1 tx -> 1 polyA site -> 1 fragment length -> r2
# #         <---r2--->......(tlen).......<s-AAAAA-e>
# #  |-----------------------------|--------------------|--------|
# #  0                               len1               len2     len3


def sample_r_pos_from_tx(tx_list, t2pa, frag_len_weight):
    """
    Given a gene, this function takes the genes' transcript name list, transcripts' polyA site list, and fragment length weight list
    Then, the function randomly choose a transcript, a polyA site, and a fragment length
    Finally, the function returns the transcript name, r2 start position, r1 end position, and a bool indicating if the r2 is the last fragment of the transcript
    """
    tx = random.choices(tx_list, k=1)[0]
    pa_index = random.randrange(0, len(t2pa[tx]))
    pa_start = t2pa[tx][pa_index]
    sampled_frag_len = random.choices(
        x_test, weights=frag_len_weight, k=1)[0]

    r2_tx_start = pa_start - sampled_frag_len + 1
    r1_tx_end = pa_start
    return tx, r2_tx_start, r1_tx_end, pa_index == len(t2pa[tx])-1


def write_rec(read_id, r1_tx_start, r1_tx_end, r2_tx_start, r2_tx_end, tx_seq, error_model, r1_out_file, r2_out_file, read_len):
    # get the tx ref seq
    # create read1 and read2
    r1_seq = tx_seq[r1_tx_start:r1_tx_end]
    r2_seq = tx_seq[r2_tx_start: r2_tx_end]

    # check if the read length is what we want
    r1_len = read_len - 60
    if len(r1_seq) == r1_len and len(r2_seq) == read_len:

        r1_rec = SeqRecord(
            'A'*60 + r1_seq.reverse_complement(),
            id=read_id + '/1',
            description='')
        # simulate the read
        r1_rec = error_model.introduce_error_scores(
            r1_rec, 'forward')
        r1_rec.seq = error_model.mut_sequence(
            r1_rec, 'forward')

        r1_rec = r1_rec[60:]

        SeqIO.write(sequences=r1_rec,
                    handle=r1_out_file, format="fastq")

        r2_rec = SeqRecord(
            r2_seq,
            id=read_id+'/2',
            description='')
        # simulate the read
        r2_rec = error_model.introduce_error_scores(
            r2_rec, 'forward')
        r2_rec.seq = error_model.mut_sequence(
            r2_rec, 'forward')

        SeqIO.write(sequences=r2_rec, handle=r2_out_file, format="fastq")


def sample_reads_from_gene_count(batch_id, out_dir, gname_list, gene_count_dict, gene_tx_list, read_len, t2pa, spliceu_txome, frag_len_weight, Stx2ee, Utx2epos, error_model, max_num_fail_per_read=500, allow_fundamental_ambiguity=False):
    """
    Given a gene, this function samples the required number of spliced and unspliced reads
    """
    MAX_ITER = 100000

    r1_out_path = os.path.join(out_dir, f"r1_{batch_id}.fastq.gz")
    r2_out_path = os.path.join(out_dir, f"r2_{batch_id}.fastq.gz")
    # if out file exists, delete it

    # if out file exists, delete it
    if os.path.exists(r1_out_path):
        os.remove(r1_out_path)
    if os.path.exists(r2_out_path):
        os.remove(r2_out_path)

    r1_out_file = bgzf.BgzfWriter(r1_out_path, "wb")
    r2_out_file = bgzf.BgzfWriter(r2_out_path, "wb")

    batch_results = []
    for gname in gname_list:
        if gene_tx_list.get(gname) is None:
            continue

        num_s, num_u = gene_count_dict[gname]
        tx_list_s_and_u = gene_tx_list[gname]
        # Spliced
        is_spliced = True
        num_fail_s = 0
        num_success_s = 0
        num_iter = 0
        max_num_fail = num_s * max_num_fail_per_read
        reads2ref_map = []
        # given a gene, we call the sample_r_pos_from_tx to get a read position
        # if the read is ambiguous (exonic), we succeed, otherwise we fail
        # we repeat this process until we get the required number of reads
        if len(tx_list_s_and_u[is_spliced]) != 0:
            while len(reads2ref_map) < num_s and num_fail_s < max_num_fail and num_iter < MAX_ITER:
                num_iter += 1
                # we sample a read
                tx, r2_tx_start, r1_tx_end, is_from_tail = sample_r_pos_from_tx(
                    tx_list_s_and_u[is_spliced], t2pa, frag_len_weight)

                r2_tx_end = r2_tx_start + read_len

                # we check if the read is within the transcript
                if (r2_tx_start < 0) or (r1_tx_end > Stx2ee[tx][-1]):
                    num_fail_s += 1
                    continue
                # we make sure the read do not span ee junctions
                #               ee#1                     ee#2
                # tx ============|=========================|==============>>
                #       r2start----->r2end
                elif any((r2_tx_start < i) and (i < r2_tx_end) for i in Stx2ee[tx]):
                    num_fail_s += 1
                    continue
                else:
                    if allow_fundamental_ambiguity:
                        rname = f"{tx}-S-{len(reads2ref_map)}"
                        write_rec(rname, r1_tx_end - (read_len - 60), r1_tx_end, r2_tx_start, r2_tx_start +
                                  read_len, spliceu_txome[tx], error_model, r1_out_file, r2_out_file, read_len)
                        reads2ref_map.append((rname, tx))
                    # we do not want to include terminal polyA reads
                    elif not is_from_tail:
                        # check if r1 end and r2 end is in the same exon
                        #               ee#1                     ee#2
                        # tx ============|=========================|==============>>
                        #                  r2----->    <-----r1
                        if any((r2_tx_start < i) and (i < r1_tx_end) for i in Stx2ee[tx]):
                            rname = f"{tx}-S-{len(reads2ref_map)}"
                            write_rec(rname, r1_tx_end - (read_len - 60), r1_tx_end, r2_tx_start, r2_tx_start +
                                      read_len, spliceu_txome[tx], error_model, r1_out_file, r2_out_file, read_len)
                            reads2ref_map.append((rname, tx))

                    num_success_s += 1

        # Unspliced
        is_spliced = False
        num_fail_u = 0
        num_success_u = 0
        num_iter = 0
        max_num_fail = num_u * max_num_fail_per_read
        if len(tx_list_s_and_u[is_spliced]) != 0:
            while len(reads2ref_map) < num_u and num_fail_u < max_num_fail and num_iter < MAX_ITER:
                num_iter += 1
                tx, r2_tx_start, r1_tx_end, is_from_tail = sample_r_pos_from_tx(
                    tx_list_s_and_u[is_spliced], t2pa, frag_len_weight)
                r2_tx_end = r2_tx_start + read_len

                # we check if the read is within the transcript
                if (r2_tx_start < 0) or (r2_tx_start + read_len > Utx2epos[tx][-1][1]):
                    num_fail_u += 1
                    continue
                # we make sure the read stay in an exon
                if any((pos[0] <= r2_tx_start) and (r2_tx_end <= pos[1]) for pos in Utx2epos[tx]):
                    if allow_fundamental_ambiguity:
                        rname = f"{tx}-{len(reads2ref_map)}"
                        write_rec(rname, r1_tx_end - (read_len - 60), r1_tx_end,
                                  r2_tx_start, r2_tx_start + read_len, spliceu_txome[tx], error_model, r1_out_file, r2_out_file, read_len)

                        reads2ref_map.append((rname, tx))
                    # we do not want to include terminal polyA reads
                    elif not is_from_tail:
                        # we do not want to include reads whose polA and r2 are in the same exon
                        #    e1start     e2start                   e2end
                        # tx |============|=========================|==============>>
                        #                      r2----->    <-----r1
                        if not any((e_start <= r2_tx_start) and (r1_tx_end <= e_end) for e_start, e_end in Utx2epos[tx]):
                            rname = f"{tx}-{len(reads2ref_map)}"
                            write_rec(rname, r1_tx_end - (read_len - 60), r1_tx_end,
                                      r2_tx_start, r2_tx_start + read_len, spliceu_txome[tx], error_model, r1_out_file, r2_out_file, read_len)

                            reads2ref_map.append((rname, tx))
                    num_success_u += 1
                else:
                    num_fail_u += 1
                    continue

        print("gene: ", gname, "num_reads", len(reads2ref_map), "num_s: ", num_s, "final_s: ", num_success_s, "num_u: ", num_u,
              "final_u: ", num_success_u, "num_fail_s: ", num_fail_s, "num_fail_u: ", num_fail_u)
        batch_results.append(
            [reads2ref_map, num_success_s+num_success_u, num_s + num_u])

    r1_out_file.close()
    r2_out_file.close()
    return batch_results

# In this function, we take all genes and required items to call the following parallel


def sample_reads_from_gene_count_parallel(gene_count_dict, gene_tx_list, t2pa, spliceu_txome, frag_len_weight, Stx2ee, Utx2epos, error_model, out_dir, num_threads=31, read_len=91, max_num_fail_per_read=500, allow_fundamental_ambiguity=False, random_seed=1):

    start_time = time.time()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    random.seed(random_seed)
    gene_list = list(gene_count_dict.keys())
    random.shuffle(gene_list)
    # Calculate the size of each sublist (except possibly the last one)
    sublist_size = math.ceil(len(gene_list) / num_threads)

    # Create the sublists
    sublists = [gene_list[sublist_size*i:min(sublist_size*(i+1), len(gene_list))]
                for i in range(num_threads)]

    reads_by_gene_list = Parallel(n_jobs=num_threads)(delayed(sample_reads_from_gene_count)(batch_id, out_dir, gname_list, gene_count_dict, gene_tx_list,
                                                                                            read_len, t2pa, spliceu_txome, frag_len_weight, Stx2ee, Utx2epos, error_model, max_num_fail_per_read, allow_fundamental_ambiguity) for batch_id, gname_list in enumerate(sublists))
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time} seconds ---")
    print(f" simulated for all genes")

    # we put all batched reads into the same file
    # first open a final fastq.gz file
    r1_out_path = os.path.join(out_dir, f"sim_rlen{read_len}_r1.fastq.gz")
    r2_out_path = os.path.join(out_dir, f"sim_rlen{read_len}_r2.fastq.gz")
    # if out file exists, delete it
    if os.path.exists(r1_out_path):
        os.remove(r1_out_path)
    if os.path.exists(r2_out_path):
        os.remove(r2_out_path)

    r1_out_file = bgzf.BgzfWriter(r1_out_path, "wb")
    r2_out_file = bgzf.BgzfWriter(r2_out_path, "wb")

    for batch_id in range(num_threads):
        r1_batch_path = os.path.join(out_dir, f"r1_{batch_id}.fastq.gz")
        r2_batch_path = os.path.join(out_dir, f"r2_{batch_id}.fastq.gz")
        with bgzf.BgzfReader(r1_batch_path) as r1_batch_file:
            for record in SeqIO.parse(r1_batch_file, "fastq"):
                SeqIO.write(sequences=record,
                            handle=r1_out_file, format="fastq")
        with bgzf.BgzfReader(r2_batch_path) as r2_batch_file:
            for record in SeqIO.parse(r2_batch_file, "fastq"):
                SeqIO.write(sequences=record,
                            handle=r2_out_file, format="fastq")
        os.remove(r1_batch_path)
        os.remove(r2_batch_path)

    r1_out_file.close()
    r2_out_file.close()

    return reads_by_gene_list


spliceu_txome = {record.id: record.seq for record in SeqIO.parse(
    spliceu_txome_fa_path, "fasta")}

error_model = kde.KDErrorModel(err_model_path)

simulated_data_not_allow_fundamental_ambiguity_out_dir = "/fs/nexus-projects/sc_frag_len/nextflow/simulation/simulated_data_not_allow_fundamental_ambiguity"
reads2ref_map_list = sample_reads_from_gene_count_parallel(
    S_U_count, gene_tx_list, t2pa, spliceu_txome, frag_len_weight, Stx2ee, Utx2epos, error_model, simulated_data_not_allow_fundamental_ambiguity_out_dir, num_threads=31, read_len=91, max_num_fail_per_read=10000, allow_fundamental_ambiguity=False, random_seed=1)

simulated_data_allow_fundamental_ambiguity_out_dir = "/fs/nexus-projects/sc_frag_len/nextflow/simulation/simulated_data_allow_fundamental_ambiguity"
reads2ref_map_list = sample_reads_from_gene_count_parallel(
    S_U_count, gene_tx_list, t2pa, spliceu_txome, frag_len_weight, Stx2ee, Utx2epos, error_model, simulated_data_allow_fundamental_ambiguity_out_dir, num_threads=31, read_len=91, max_num_fail_per_read=10000, allow_fundamental_ambiguity=True, random_seed=1)
