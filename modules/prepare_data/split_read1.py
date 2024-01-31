#!/usr/bin/env python3

import tqdm
import dnaio
import argparse
import pandas as pd
import json
import os

parser = argparse.ArgumentParser(description="Split reads.")
parser.add_argument("SRR", type=str,
                     help="sample SRR ID")
parser.add_argument("read1_path", type=str,
                     help="path to read1 fastq file")
parser.add_argument("chemistry",type=str,help="chemistry version")
parser.add_argument("out_dir",type=str,help="output directory. Usually named by GSE ID")

args = parser.parse_args()
os.makedirs(args.out_dir, exist_ok = True)

# barcode is either 10 (v2) or 12 (v3) bases
barcode_i = 16

if args.chemistry == 'v2':
    umi_i = 10
else:
    umi_i = 12
umi_i += barcode_i # UMI is 16 bases

polyT_i = umi_i + 32 # polyT is 30 bases + VN is 2

with dnaio.open(args.read1_path) as reader, \
    dnaio.open(os.path.join(args.out_dir, "_".join([args.SRR,"1.fastq.gz"])),
            os.path.join(args.out_dir, "_".join([args.SRR,"1_pe.fastq.gz"])),
            mode="w") as writer1:
    
    for i, (record1) in enumerate(tqdm.tqdm(reader)):
        writer1.write(record1[:umi_i], record1[polyT_i:])

