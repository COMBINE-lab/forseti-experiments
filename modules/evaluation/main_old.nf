
// Run starsolo using the spliceu reference and output the transcriptome BAM file.
process run_starsolo {
    publishDir "$params.output_dir/$GSE/run_starsolo", mode: 'symlink', pattern: "$SRR"
    input:
        tuple val(species),
            val(tissue),
            val(chemistry), 
            val(GSE),
            val(SRR),
            path(fastp_json),     
            path(genome_path),
            path("fastp_dir"),
            path("index_dir")
    output:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path(fastp_json),
            path("$SRR/Aligned.sortedByCoord.out.bam"),
            path("$SRR/Aligned.toTranscriptome.out.bam"), emit: run_starsolo
        path "$SRR"


    script:
        outdir="$SRR"
        umilen = if (chemistry == "v2") 10 else 12

        """
        mkdir -p $outdir

        reads1="\$(find -L fastp_dir -name "*_1.fastq.gz" -type f | sort| paste -sd,)"

        reads2="\$(find -L fastp_dir -name "*_2.fastq.gz" -type f | sort| paste -sd,)"

        ulimit -n 2048

        STAR --runThreadN $params.num_threads --readFilesIn \$reads2 \$reads1 \
        --readFilesCommand zcat \
        --genomeDir index_dir --outFileNamePrefix ${SRR} \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_UMI_Simple \
        --soloUMIlen $umilen \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist $params.whitelist_dir/10x_${chemistry}_permit.txt \
        --soloFeatures GeneFull \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --limitIObufferSize 50000000 50000000 \
        --quantMode TranscriptomeSAM 
        """        
}

// Yuan. Jan 23 
process simulation {
    // TODO: publishDir "$params.output_dir/simulation/err_model", mode: 'symlink', pattern: 
    input:
    ouput:
    script:

        """
        # set up env and simpleaf
        # conda create env -f simpleaf.yml
        # conda activate simpleaf

        bash run_simpleaf.sh
        ##### TODO:  How to modify jsonette?

        1k_pbmc_star_dir=${params.output_dir}/simulation/star_out
        mkdir -p $1k_pbmc_star_dir

        data_dir=${params.output_dir}/simulation/simpleaf_workdir/data/human_pbmc_1k_v3_fastqs

        # build err model with iss require aligned bam file with MD tag and paired-end in same read length
        STAR --runThreadN 8 \
        --readFilesIn ${data_dir}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz ${data_dir}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
        --readFilesCommand zcat \
        --genomeDir ${params.data_dir}/star_spliceu_index/refdata-gex-GRCh38-2020-A \
        --outFileNamePrefix ${1k_pbmc_star_dir}/real_data_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMtlen 2 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --alignSplicedMateMapLminOverLmate 0 \

        samtools sort ${1k_pbmc_star_dir}/real_data_Aligned.sortedByCoord.out.bam -o ${1k_pbmc_star_dir}/real_data_sorted.bam
        samtools index ${1k_pbmc_star_dir}/real_data_sorted.bam

        # # build err model with iss
        wd="/fs/nexus-projects/sc_frag_len/nextflow/simulation"
        iss model -b ${1k_pbmc_star_dir}/real_data_sorted.bam -o err_model

        # # Simulate reads
        spliceu_ref_human_dir=${params.data_dir}/spliceu_ref/refdata-gex-GRCh38-2020-A
        out_dir=${params.output_dir}/simulation/simulated_data 
        mkdir -p $out_dir

        python simulator.py \
            --spliceu_ref_human_dir $spliceu_ref_human_dir \
            --simpleaf_outdir ${params.output_dir}/workflow_output \
            --frag_len_model_path \
            --err_model_path ${1k_pbmc_star_dir}/err_model.npz \
            --out_dir $out_dir 

    """
}

