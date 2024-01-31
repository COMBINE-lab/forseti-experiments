#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict
import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import random
import json
from Bio import SeqIO
import regex

parser = argparse.ArgumentParser(
    description="Count and summarize non-zero length reads.")
parser.add_argument("bam", type=str, help="analyze this bam file")
parser.add_argument("gse_id", type=str, help="sample gse id")
parser.add_argument("srr_id", type=str, help="sample srr id")
parser.add_argument("fp_file", type=str, help="fastp json file")
parser.add_argument("chem_version", type=str, help="chemistry_version")
parser.add_argument("genome_file", type=str, help="genome FASTA file")
parser.add_argument("savedir", type=str, help="where to save the results")

args = parser.parse_args()
# args = parser.parse_args([
#     "/fs/nexus-projects/sc_frag_len/nextflow/workflow_output/preprocessing/run_star/GSE135922/SRR9990661/pe_Aligned.sortedByCoord.out.bam",
#     "GSE135922",
#     "SRR9990661",
#     "/fs/nexus-projects/sc_frag_len/nextflow/workflow_output/preprocessing/fastp/GSE135922/SRR9990661/SRR9990661_fastp.json",
#     "v3",
#     "/fs/nexus-projects/sc_frag_len/nextflow/data/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
#     "/fs/nexus-projects/sc_frag_len/nextflow/Jan4_seSTAR_check/"
# ])
# args = parser.parse_args([
#     "/fs/nexus-projects/sc_frag_len/nextflow/workflow_output/preprocessing/run_star/GSE125970/SRR8513797/pe_Aligned.sortedByCoord.out.bam",
#     "GSE125970",
#     "SRR8513797",
#     "/fs/nexus-projects/sc_frag_len/nextflow/workflow_output/preprocessing/fastp/GSE125970/SRR8513797/SRR8513797_fastp.json",
#     "v2",
#     "/fs/nexus-projects/sc_frag_len/nextflow/data/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
#     "/fs/nexus-projects/sc_frag_len/nextflow/Jan4_seSTAR_check/GSE125970"
# ])

gse_id = args.gse_id
srr_id = args.srr_id
fp_file = args.fp_file
genome_file = args.genome_file
savedir = args.savedir
print(args.bam)


save = pysam.set_verbosity(0)
bamfile = pysam.AlignmentFile(args.bam, "rb")
pysam.set_verbosity(save)

# Reading the FASTA file
genome = {record.id: record.seq for record in SeqIO.parse(
    genome_file, "fasta")}

tlens = []
# rdr = bamfile.fetch(until_eof=True)

"""
Intended logic:
add back the soft clipped seq, ignore skip N, get the new tlen
"""


def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        # if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
        # continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def correct_tlen(read1, read2, r1_len, r2_len):
    '''
    1. Add back soft clipped seq by STAR
    2. Handling shared skipped region(represent by N in CIGAR) 
    '''
    polyt = ''
    bg = ''
    ri = random.randint(-500, 500)
    # R1 maps on the forward strand, R2 maps on the reverse strand
    if read1.flag == 99 and read2.flag == 147:
        left_read = read1
        right_read = read2
        left_read_fix_len = r1_len
        right_read_fix_len = r2_len

        # find the 30 bases
        # we are working on genomic sequence. 
        # the orientation of the 30 bases follows R1
        # so they are always the reverse of R2 orientation
        # therefore, we want to do RC to make then follow the R2 orientation
        polyt_end = read1.reference_start

        if read1.cigartuples[0][0] == 4:
            polyt_end -= read1.cigartuples[0][1]
        polyt = genome[read1.reference_name][(
            polyt_end-30):polyt_end].reverse_complement()
        bg = genome[read1.reference_name][(
            polyt_end-ri-30):(polyt_end-ri)].reverse_complement()

    # R1 maps on the reverse strand, R2 maps on the forward strand
    elif read2.flag == 163 and read1.flag == 83:
        left_read = read2
        right_read = read1
        left_read_fix_len = r2_len
        right_read_fix_len = r1_len

        # find the 30 bases
        polyt_start = read1.reference_end

        if read1.cigartuples[-1][0] == 4:
            polyt_start += read1.cigartuples[-1][1]
        polyt = genome[read1.reference_name][(
            polyt_start):(polyt_start + 30)]
        bg = genome[read1.reference_name][(
            polyt_start+ri):(polyt_start+ri+30)]

    '''logic
    find the distance(dis) between: reference end of LEFT reads and reference start of RIGHT reads
    new tlen = dis+ len(R1)+len(R2).
        + special case: for distance = negative value:
          if r1 and r2 share same skipped region(1 or 2)(r1,r2 have same length of N region). add len(N) for dis.
    '''
    # -------Find LEFT_read real end (by adding back soft clip)
    left_read_end = left_read.reference_end
    if left_read.cigartuples[-1][0] == 4:  # 4 for soft clip
        left_read_end += left_read.cigartuples[-1][1]

    # -------Find Right_read real start
    right_read_start = right_read.reference_start
    if right_read.cigartuples[0][0] == 4:  # 4 for soft clip
        right_read_start -= right_read.cigartuples[0][1]

    dis = right_read_start - left_read_end

    def find_N_info_cigar(cigar_tuple):
        find_N_list = []
        for i in range(len(cigar_tuple)):
            if cigar_tuple[i][0] == 3:
                find_N_list.append(cigar_tuple[i][1])
        return find_N_list

    # -------new tlen
    new_tlen = dis + left_read_fix_len + right_read_fix_len

    if 'N' in left_read.cigarstring and 'N' in right_read.cigarstring:
        find_left_N = find_N_info_cigar(left_read.cigartuples)
        find_right_N = find_N_info_cigar(right_read.cigartuples)
        i_max = min(len(find_left_N), len(find_right_N))
        if find_left_N[-1] == find_right_N[0]:
            N_len = find_left_N[-1]
            new_tlen += N_len
        elif i_max > 1 and find_left_N[-1] == find_right_N[1] and find_left_N[-2] == find_right_N[0]:
            new_tlen += find_left_N[-1]
            new_tlen += find_left_N[-2]
        elif i_max > 2 and find_left_N[-1] == find_right_N[2] and find_left_N[-2] == find_right_N[1] and find_left_N[-3] == find_right_N[0]:
            new_tlen += find_left_N[-1]
            new_tlen += find_left_N[-2]
            new_tlen += find_left_N[-3]

    return new_tlen, polyt, bg


f_tlens = open(os.path.join(savedir, "{}_new_tlens.txt".format(srr_id)), "w")

fp_file = open(fp_file)
fp_report = json.load(fp_file)
read1_rest_length = fp_report["read1_before_filtering"]["total_cycles"]
read2_length = fp_report["read2_before_filtering"]["total_cycles"]


polyt_seq = set()
polya_seq = set()
polyt_bg_seq = set()
polya_bg_seq = set()

num_reads = 0

for read1, read2 in read_pair_generator(bamfile):
    num_reads += 1
    new_tlen = 0

    new_tlen, polyt, bg = correct_tlen(
        read1, read2, read1_rest_length, read2_length)

    # we use a very loose filter here
    if "TTT" in polyt:
        if len(polyt) == 30:
            polyt_seq.add(polyt)
        if len(bg) == 30:
            polyt_bg_seq.add(bg)
    if "AAA" in polyt:
        if len(polyt) == 30:
            polya_seq.add(polyt)
        if len(bg) == 30:
            polya_bg_seq.add(bg)

    tlens.append(new_tlen)
    f_tlens.write(str(new_tlen)+"\n")

f_tlens.close()

poly = {
    "polyt_seq": polyt_seq,
    "polya_seq": polya_seq,
    "polyt_bg_seq": polyt_bg_seq,
    "polya_bg_seq": polya_bg_seq,
}

# define the directory to save the priming site sequences
psdir = os.path.join(savedir, "priming_site_seqs")
os.makedirs(psdir, exist_ok=True)

# TODO: Save the polyT to file or just pass it?
for lst, cont in poly.items():
    file_name = os.path.join(psdir, f"{lst}.txt")
    with open(file_name, 'w') as file:
        for item in cont:
            file.write(f"{item}\n")

tlens = np.array(tlens)

tlens_80_100 = tlens[(80 < tlens) & (tlens < 100)]
tlens_100_500 = tlens[(100 < tlens) & (tlens < 500)]
tlens_small = tlens[tlens < 500]
tlens_0_1000 = tlens[tlens < 1000]

tlens_80_100_ratio = tlens_80_100.size / tlens.size
tlens_100_500_ratio = tlens_100_500.size / tlens.size
tlens_small_ratio = tlens_small.size / tlens.size
tlens_0_1000_ratio = tlens_0_1000.size/tlens.size

tlens_medium = tlens[(500 < tlens) & (tlens < 3000)]
tlens_medium_ratio = tlens_medium.size / tlens.size
tlens_large = tlens[3000 < tlens]
tlens_large_ratio = tlens_large.size / tlens.size

with open(os.path.join(savedir, "ratios-{}.txt".format(gse_id)), "w") as file:
    file.writelines(["tlens_small_ratio: {} \n".format(tlens_small_ratio),
                     "tlens_medium_ratio: {}\n".format(tlens_medium_ratio),
                     "tlens_large_ratio: {}\n".format(tlens_large_ratio),
                     ])


def plot_fig_classic(nums, bins, xlab, ylab, suptitle, title, save):
    plt.figure()
    plt.hist(nums, bins)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.suptitle(suptitle)
    plt.savefig(os.path.join(savedir, save))
    plt.close()


plot_fig_classic(tlens_80_100, 19, "Length", "Count", "Template lengths (80-100) - " +
                 gse_id, tlens_80_100_ratio, "classic-tlen-{}-80-100.png".format(gse_id))
plot_fig_classic(tlens_100_500, 398,  "Length", "Count", "Template lengths (100-500) - " +
                 gse_id, tlens_100_500_ratio, "classic-tlen-{}-100-500.png".format(gse_id))
plot_fig_classic(tlens_small, 499, "Length", "Count", "Template lengths (small) - " +
                 gse_id, tlens_small_ratio, "classic-tlen-{}-small.png".format(gse_id))


def plot_fig(nums, range_text):
    sns.set_style('darkgrid')
    plt.figure()
    ax = sns.distplot(nums, kde=True)
    ax.set(title=f"Template lengths ({range_text}) - " +
           gse_id, xlabel='Fragment length', ylabel='Density')
    plt.savefig(os.path.join(savedir, f"tlen-{gse_id}-{range_text}.png"))
    plt.show()


plot_fig(tlens, 'all')
plot_fig(tlens_80_100, '80-100')
plot_fig(tlens_100_500, '100-500')
plot_fig(tlens_0_1000, '0-1000')
plot_fig(tlens_small,  'small')
plot_fig(tlens_medium, 'medium')
plot_fig(tlens_large, 'large')
