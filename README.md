# forseti-experiments

This repo includes the workflow pipeline used for generating the results in the manuscript "Forseti: A mechanistic and predictive model of the splicing status of scRNA-seq reads" (link coming soon). 

**To ease the future testing and validation, we have included the trained binding affinity model, and the output of the evaluation scripts in the `workflow_output` directory.** 

The pipeline is implemented in Nextflow and you will need to have Nextflow installed in your system to run the pipeline. If not, please check https://www.nextflow.io/docs/latest/getstarted.html for details. The pipeline is designed to run on a SLURM cluster, but can be easily adapted to other cluster systems. 

## Introduction
This workflow can be devided into six main steps:
1. **Prepeare the data:** In this step we fetch the sequencing reads from 10 selected paired-end scRNA-seq experiments and their corresponding reference sets from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads), and build the _spliceu_ augmented reference, in which each transcript has a spliced and a unspliced version.
2. **Process the data:** In this step, we process the sequencing reads to compute the fragment length of the reads, and obtain the potentail priming windows of oligo(dT) primers indicated by the paired-end alignments.
3. **Fit and evaluate the key components:** In this step, we use the emperical priming windows and fragment lenths obtained from the previous step to (1) train a multi-layer perceptron (MLP) using the priming windows to predict the binding affinity of oligo(dT) primer for given sequences the context of scRNA-seq, and (2) fit an emperical fragment length distribution using a cubic spline to represent the emperical distribution of the cDNA fragment length in scRNA-seq.
4. **Simulate data:** Using the the _spliceu_ augmented reference generated in step 1, we simulate two sets of sequencing reads with and without fundamental ambiguity, respectively. Here, fundamentally ambiguious reads are the reads that come from either polyA tail priming and the read2 is contained within the terminal exon, or the reads where the original cDNA fragment arising from an exon. As a mechanistic model, Forseti is unable to distinguish between these two sources of fundamentally ambiguous reads.
5. **Evaluate Forseti:** We evaluate the performance of Forseti on the two sets of simulated reads, as well as two experimental test datasets randomly selected from the 10 processed datasets. 

## Running the pipeline
To run the pipeline, we need to make sure we have Nextflow installed in your system. Please make sure that you have 3 TB of disk space to run the workflow. Then, we need to fetch the repo from GitHub using the following **shell** command:

```bash
git clone https://github.com/COMBINE-lab/forseti-experiments.git 
cd forseti-experiments
```

Now, we are in the directory of the workflow. Based on the slurm setting of your environment, you will need to change the specific configuration about the slurm setting in the `nextflow.config` file. For example, you will need to change the `process` section to fit your slurm configuration, especially the `clusterOptions` flag. You might also want to set the Memory and number of threads to an appropirate value. We recommend you to use as many threads as possible to speed up the pipeline.

With the configuration set up, we now just need to call the 6 steps introduced above individually using the following commands. **NOTICE THAT** the following commands need to be execute one by one. Please wait the previous step to finish before running the next step.

```bash

nextflow main.nf -resume --mode prepare_data
nextflow main.nf -resume --mode process_data
nextflow main.nf -resume --mode fit_models
nextflow main.nf -resume --mode simulation
nextflow main.nf -resume --mode evaluation

```

### Tips:
1. Nextflow cannot make sure the fetched FASTQ files are not corrupted (see [this issue](https://github.com/nextflow-io/nextflow/issues/2491)). If you find the pipeline failed at the `prepare_data` step (because we check the fetched file internally), you will want to check `data/datasets` and remove the lines corresponding to the successfully processed files (SRR numbers in the directory) from `input_files/sample_url_sheet.csv` and re-run the pipeline.
    
    For example, if you see that `data/datasets/GSE122357/SRR8039659`, then you want to remove the line #11 in `input_files/sample_url_sheet.csv` from the file. The line is 
    ```
    GSE122357,SRR8039659,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR803/009/SRR8039659/SRR8039659_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR803/009/SRR8039659/SRR8039659_2.fastq.gz,Mouse,Brain,v2,0,151
    ``` 

2. We specified `-resume` in the command line to avoid recomputing the results. If you want Nextflow to forget the previous results and rerun a step from the beginning, you can remove the `-resume` flag.







