# forseti-experiments

This repo includes the workflow pipeline used for generating the results in the manuscript "Forseti: A mechanistic and predictive model of the splicing status of scRNA-seq reads" (link coming soon). 

**To ease the future testing and validation, we have included the trained binding affinity model, and the output of the evaluation scripts in the `workflow_output` directory.** 

The pipeline is implemented in Nextflow and you will need to have Nextflow installed in your system to run the pipeline. If not, please check https://www.nextflow.io/docs/latest/getstarted.html for details. The pipeline is designed to run on a SLURM cluster, but can be easily adapted to other cluster systems. 

## Introduction
This workflow can be devided into six main steps:
1. Prepeare the data, including the sequencing reads from scRNA-seq experiments, and the spliceu augmented reference, in which each transcript has a spliced and a unspliced version.
2. Process the data, including computing the fragment length of the reads, and obtaining the potential priming windows for alignment.
3. fit the model, including training the binding affinity model (a multi-layer perceptron, MLP) for the binding affinity of oligo(dT) primer in the context of scRNA-seq, and fitting the fragment length distribution model using a cubic spline representing the distribution of the cDNA fragment length in scRNA-seq.
4. Simulate data: we simulate reads from the _spliceu_ augmented reference generated in step 1, including sequencing reads with and without fundamental ambiguity. Here, fundamentally ambiguious reads are the reads that come from either polyA tail priming or the whole cDNA fragment arising from an exon. As a mechanistic model, Forseti is unable to distinguish between these two sources of fundamentally ambiguous reads.
5. Evaluate Forseti on the simulated data, and the two experimental test datasets randomly selected from the processed datasets.

## Running the pipeline
To run the pipeline, we need to make sure we have Nextflow installed in your system. Then, we need to fetch the repo from GitHub using the following **shell** command:

```bash
git clone https://github.com/COMBINE-lab/forseti-experiments.git 
cd forseti-experiments
```

Now, we are in the directory of the workflow. We need to call the 6 steps introduced above individually. For example, to run the first step, we can use the following shell command. **NOTICE THAT** the following commands need to be execute one by one. Please wait the previous step to finish before running the next step.

```bash

nextflow main.nf -resume --mode prepare_data
nextflow main.nf -resume --mode process_data
nextflow main.nf -resume --mode fit_model
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







