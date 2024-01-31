#!/usr/bin/env python3
import argparse
import pysam
import pandas as pd
import os
from collections import defaultdict
import pyranges as pr
import time
import pickle
start_time = time.time()

parser = argparse.ArgumentParser(
    description="classify the trancripts source of read")
parser.add_argument("GSE", type=str, help="sample gse id")
parser.add_argument("SRR", type=str, help="sample srr id")
parser.add_argument("se_txome_bam", type=str,
                    help="paired-end transcriptome based bam file")
parser.add_argument("pe_txome_bam", type=str,
                    help="paired-end transcriptome based bam file")
parser.add_argument("r1_r2_exonic", type=str,
                    help="A bed file contianing exonic read1 & read2")
parser.add_argument("t2g_file", type=str,
                    help="A tsv file contianing the transcript to gene mapping")
parser.add_argument("exon_dist_thresh", type=int, default=1000,
                    help="A tsv file contianing the transcript to gene mapping")
parser.add_argument("outdir", type=str, help="where to save the results")

args = parser.parse_args()

# args = parser.parse_args(["GSE130636",
#                           "SRR9004344",
#                           f"se_Aligned.toTranscriptome.out.bam",
#                           f"pe_Aligned.toTranscriptome.out.bam",
#                           f"r1_r2_exonic.bed",
#                           "spliceu_dir/t2g_3col.tsv",
#                           f"SRR9004344"])
# work_out = "/fs/nexus-projects/sc_frag_len/nextflow/find_r1/GSE125970"
# args = parser.parse_args(["GSE125970",
#                           "SRR8513797",
#                           f"{work_out}/se_Aligned.toTranscriptome.out.bam",
#                           f"{work_out}/pe_Aligned.toTranscriptome.out.bam",
#                           f"{work_out}/r1_r2_exonic.bed",
#                           "/fs/nexus-projects/sc_frag_len/nextflow/find_r1/SRR9990661/refdata-gex-GRCh38-2020-A/t2g_3col.tsv",
#                           f"{work_out}/dh_jan17"])

GSE = args.GSE
SRR = args.SRR
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

# read bam header to create the tid2name dict
pe_bam = pysam.AlignmentFile(args.pe_txome_bam, "rb")
tid2name = {tid: tname for tid, tname in enumerate(pe_bam.references)}
tname2id = {tname: tid for tid, tname in enumerate(pe_bam.references)}
pe_bam.close()

# read in the t2g file
t2g_df = pd.read_csv(args.t2g_file, sep="\t", header=None,
                     names=["tid", "gid", "splicing_status"])
t2g = t2g_df.set_index("tid").to_dict()["gid"]
t2s = t2g_df.set_index("tid").to_dict()["splicing_status"]
del t2g_df
# 1---------Get single-exon transcripts from gtf file
# we read in the gtf file

# if a spliced tx does not have its unspliced version, it is a single exon tx
# we might want to include them in the test case because they will be more likely to 
# take the reads bc they are short and we give them a polya tail
single_exon_tx = set([tx for tx in t2g.keys() if not tx.endswith('-U') and tx + '-U' not in t2g.keys()])

# get the fraction of single exon tx in all tx-S: 
# number of single exon tx / number of single-exon tx + number of spliced tx
single_exon_tx_freq = len(single_exon_tx) / ((len(tname2id) - len(single_exon_tx))/2 + len(single_exon_tx))

"""
We want to find the reads that maps to a single gene and can get a definitive splicing status because of its read1
that is, 
1. they should be map to the single gene
2. their exonic_r2 should be exonic
3. for spliced reads, their exonic_r1 should be either exon-exon junctional reads or contained within another exon of the same transcript
4. for unspliced reads, their exonic_r1 should be either exon-intron junctional reads, or contained within the intron of the same transcript. Moreover, the exonic_r1 should not be in the exon of any other transcript.

So, we first get the set of references of each exonic_r1 and exonic_r2 separately, and then compare them
"""

# we read process each line of the tsv file to fill the exonic_r1 and exonic_r2 dict

exonic_r1 = {}
exonic_r2 = {}

# we read in the r1_r2_exonic file
# THE 3,13,14,15 columns are: read_id, exon_start, exon_end, mapped_tx_id
# exonic_r1 and exonic_r2 are identified according to the end of the read name
with open(args.r1_r2_exonic, 'r') as file:
    for line in file:
        line = line.strip().split('\t')
        # split the read name to know if it is read1 or read2
        read_id_full = line[3].split('/')
        rname = read_id_full[0]
        rid = read_id_full[1]
        # if read1, add the tx to exonic_r1 dict
        if rid == '1':
            if exonic_r1.get(rname) is None:
                exonic_r1[rname] = set([(tname2id[line[15]], int(line[13]), int(line[14]))])
            else:
                exonic_r1[rname].add((tname2id[line[15]], int(line[13]), int(line[14])))
        # if read2, add the tx to exonic_r2 dict
        else:
            if exonic_r2.get(rname) is None:
                exonic_r2[rname] = set([(tname2id[line[15]], int(line[13]), int(line[14]))])
            else:
                exonic_r2[rname].add((tname2id[line[15]], int(line[13]), int(line[14])))

# from bam file, I want to get the reference list of each read
pe_read_to_ref = {rname: set() for rname in exonic_r2.keys()} # keys: read_id, values: set of references
pe_bam = pysam.AlignmentFile(args.pe_txome_bam, "rb")
for read in pe_bam.fetch(until_eof=True):
    if read.is_read2:
        read_id = read.query_name
        if read_id in pe_read_to_ref:
            pe_read_to_ref[read_id].add(read.reference_id)
pe_bam.close()
# len(pe_read_to_ref) = 4,706,773
# we get unique gene mapped reads by checking if all references of the reads mapped to the same gene
unigene_reads = {rname: t2g[tid2name[list(refs)[0]]] for rname, refs in pe_read_to_ref.items() if len(set([t2g[tid2name[tid]] for tid in refs])) == 1}
# len(unigene_reads)
# 3,596,464
##################
# classify reads
##################
# Note that all processed reads mapped to some transcripts by STAR
# so, if we see that the exonic_r1 of a read is not mapped to any exons (no in the exonic_r2 dict), its either unspliced or ee junctional reads
unspliced_or_eejunction = set(exonic_r2.keys()) - set(exonic_r1.keys()) & set(unigene_reads.keys())
# len(unspliced_or_eejunction) 186,647
double_e_or_ambiguous = set(exonic_r2.keys()) & set(exonic_r1.keys()) & set(unigene_reads.keys())
# len(double_e_or_ambiguous) 3,409,817

# we first want to get the unspliced reads
# if the reads' reference are all unspliced, then it is unspliced
# otherwise, it across an ee junction
# we get the reference list of each read from pe_read_to_ref
# then we check if all references of the read are unspliced
unspliced_reads = set()
for rname in unspliced_or_eejunction:
    paired_ref = pe_read_to_ref[rname]
    # we want to make sure that the transcript exonic_r2 mapped is the same as the transcripts exonic_r1 and exonic_r2 both mapped
    # TODO: test if we need to remove terminal exonic reads, i.e., use == or &
    if all([tid2name[tid].endswith('-U') for tid in paired_ref]) and set([tname2id.get(tid2name[tx]+"-U") for tx, start, end in exonic_r2[rname]]) & paired_ref:
        unspliced_reads.add(rname)
# 6749

ee_junction_spliced_reads = set()
# we first want to process the exon-exon junctional reads
# get candicate pool
ee_junctions = unspliced_or_eejunction - unspliced_reads

for rname in ee_junctions:
    # we want to make sure that the transcript exonic_r2 mapped is the same as the transcripts exonic_r1 and exonic_r2 both mapped
    # TODO: test if we need to remove terminal exonic reads
    paired_ref = pe_read_to_ref[rname]
    if all([not tid2name[tid].endswith('-U') for tid in paired_ref]) and set([tx for tx, start, end in exonic_r2[rname]]) == paired_ref:
        ee_junction_spliced_reads.add(rname)
# 51821

# then, we get the reads whose exonic_r1 and exonic_r2 are in the different exons of the same transcripts
# we also want to make sure the two exons is 1000 based away from each other 
double_e_spliced_reads = set()
for rname in double_e_or_ambiguous:
    # we make sure the reference list of exonic_r1 and exonic_r2 are the same
    # however, the position of the reference are different
    paired_refs = pe_read_to_ref[rname]
    r1_refs = {tid: (start, end) for tid, start, end in exonic_r1[rname] if tid in paired_refs}
    r2_refs = {tid: (start, end) for tid, start, end in exonic_r2[rname] if tid in paired_refs}
    
    if r1_refs.keys() == r2_refs.keys():
        if all([r1_exon_pos != r2_refs[tid] and abs(r1_exon_pos[1]-r2_refs[tid][0]) < args.exon_dist_thresh for tid, r1_exon_pos in r1_refs.items()]):
            double_e_spliced_reads.add(rname)
# 99660

# we write the reads to two files
se_bam = pysam.AlignmentFile(args.se_txome_bam, "rb")
double_e_spliced_reads_file = pysam.AlignmentFile(
    os.path.join(outdir, "r2_exonic_r1_another_exon_spliced_txome.bam"), "wb", template=se_bam)
ee_junction_spliced_reads_file = pysam.AlignmentFile(
    os.path.join(outdir, "r2_exonic_r1_ee_junction_spliced_txome.bam"), "wb", template=se_bam)
unspliced_reads_file = pysam.AlignmentFile(
    os.path.join(outdir, "r2_exonic_r1_unspliced_txome.bam"), "wb", template=se_bam)
for read in se_bam.fetch(until_eof=True):
    read_id = read.query_name
    if read_id in double_e_spliced_reads:
        double_e_spliced_reads_file.write(read)
    elif read_id in ee_junction_spliced_reads:
        ee_junction_spliced_reads_file.write(read)
    elif read_id in unspliced_reads:
        unspliced_reads_file.write(read)

double_e_spliced_reads_file.close()
ee_junction_spliced_reads_file.close()
unspliced_reads_file.close()

endtime = time.time()
print("total time: ", endtime - start_time)

# we write the read to gene mapping in the unigene_reads to a tsv file
with open(os.path.join(outdir, "unigene_reads_to_gene.tsv"), "w") as file:
    for rname, gene in unigene_reads.items():
        file.write(f"{rname}\t{gene}\n")

