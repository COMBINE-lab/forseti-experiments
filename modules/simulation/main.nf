process sim_simpleaf {
    input:
        path spliceu_dir
    output:
        path "simpleaf_quant/af_quant", emit: af_quant
        path "${fastq_dir}", emit: fastq_dir
    
    script:
    fastq_dir = "fastqs"
    // TODO: CHANGE THE INDEX DIR
    """
    export ALEVIN_FRY_HOME="af_home"
    mkdir -p ${fastq_dir}

    simpleaf set-paths
    ulimit -n 2048

    # download FASTQ files
    wget -qO- https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar | tar xf - --strip-components=1 -C ${fastq_dir}

    wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz
    gunzip 3M-february-2018.txt.gz

    # index
    #simpleaf index \
    #    --output simpleaf_index \
    #    --ref-seq ${spliceu_dir}/spliceu.fa \
    #    --threads ${params.num_threads}


    # quant
    simpleaf quant --chemistry 10xv3 \
        --output simpleaf_quant \
        --expected-ori fw \
        --index /fs/nexus-projects/sc_frag_len/nextflow/work/87/ae659e88fda54234c016e546c4ae71/simpleaf_index/index \
        --reads1 ${fastq_dir}/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,${fastq_dir}/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
        --reads2 ${fastq_dir}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,${fastq_dir}/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
        --resolution cr-like \
        --t2g-map ${spliceu_dir}/t2g_3col.tsv \
        --unfiltered-pl 3M-february-2018.txt \
        --threads ${params.num_threads}
    """

}

process sim_err_model {
    input:
        path fastq_dir
    output:
        path "err_model.npz"

    script:
        star_index = "${params.data_dir}/star_spliceu_index/refdata-gex-GRCh38-2020-A"
        star_out = "star_out"
        // TODO: add star in the conda env
        """
        # build err model with iss require aligned bam file with MD tag and paired-end in same read length
        STAR --runThreadN ${params.num_threads} \
        --readFilesIn ${fastq_dir}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz ${fastq_dir}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
        --readFilesCommand zcat \
        --genomeDir ${star_index} \
        --outFileNamePrefix ${star_out}/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMtlen 2 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --alignSplicedMateMapLminOverLmate 0

        samtools index ${star_out}/Aligned.sortedByCoord.out.bam

        # build err model with iss
        iss model -b ${star_out}/Aligned.sortedByCoord.out.bam -o err_model

        """
}

process sim_simulate_reads_allow_fundamental_ambiguity {
    publishDir "${params.output_dir}/simulation/allow_fundamental_ambiguity", mode: 'symlink'
    input: 
        path af_quant
        path err_model
        path "spliceu_dir"
        path star_index_dir
        path spline_model_path  
    output:
        path "${out_dir}", emit: star_dir
        path "${out_dir}/Aligned.toTranscriptome.out.bam", emit: txome_bam
        path "${sim_fastq_dir}", emit: fastqs

    script:
        sim_fastq_dir = "fastqs"
        out_dir = "star_out"

        """
        python $moduleDir/simulator.py \
        spliceu_dir/spliceu.fa spliceu_dir/spliceu.gtf \
        spliceu_dir/tx_polya.bed \
        spliceu_dir/t2g_3col.tsv \
        ${err_model} \
        ${af_quant} \
        ${spline_model_path} \
        ${sim_fastq_dir} \
        ${params.num_threads_large} \
        ${params.simulation.read_len} \
        ${params.simulation.max_num_fail_per_read}\
        ${params.simulation.random_seed}\
        --allow-fundamental-ambiguity

        STAR --runThreadN ${params.num_threads} \
            --readFilesIn ${sim_fastq_dir}/sim_rlen${params.simulation.read_len}_r2.fastq.gz \
            --readFilesCommand zcat \
            --genomeDir ${star_index_dir} \
            --outFileNamePrefix ${out_dir}/ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --quantTranscriptomeBan Singleend
        """
}


process sim_simulate_reads_not_allow_fundamental_ambiguity {
    publishDir "${params.output_dir}/simulation/not_allow_fundamental_ambiguity", mode: 'symlink'
    input: 
        path af_quant
        path err_model
        path "spliceu_dir"
        path star_index_dir
        path spline_model_path  
    output:
        path "${out_dir}", emit: star_dir
        path "${out_dir}/Aligned.toTranscriptome.out.bam", emit: txome_bam
        path "${sim_fastq_dir}", emit: fastqs

    script:
        sim_fastq_dir = "fastqs"
        out_dir = "star_out"
        
        """
        python $moduleDir/simulator.py \
        spliceu_dir/spliceu.fa spliceu_dir/spliceu.gtf \
        spliceu_dir/tx_polya.bed \
        spliceu_dir/t2g_3col.tsv \
        ${err_model} \
        ${af_quant} \
        ${spline_model_path} \
        ${sim_fastq_dir} \
        ${params.num_threads_large} \
        ${params.simulation.read_len} \
        ${params.simulation.max_num_fail_per_read}\
        ${params.simulation.random_seed}

        STAR --runThreadN ${params.num_threads} \
            --readFilesIn ${sim_fastq_dir}/sim_rlen${params.simulation.read_len}_r2.fastq.gz \
            --readFilesCommand zcat \
            --genomeDir ${star_index_dir} \
            --outFileNamePrefix ${out_dir}/ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --quantTranscriptomeBan Singleend
        """
}


workflow simulation {
    main:
        // define paths
        spliceu_dir = file("${params.data_dir}/spliceu_ref/refdata-gex-GRCh38-2020-A", checkIfExists: true)
        star_index_dir = file("${params.data_dir}/star_spliceu_index/refdata-gex-GRCh38-2020-A", checkIfExists: true)
        spline_model_path = file("${params.output_dir}/models/spline_model/spline_model.pkl", checkIfExists: true)
        
        // processes
        sim_simpleaf(spliceu_dir)
        sim_err_model(sim_simpleaf.out.fastq_dir)
        sim_simulate_reads_allow_fundamental_ambiguity(
            sim_simpleaf.out.af_quant, 
            sim_err_model.out, 
            spliceu_dir,
            star_index_dir,
            spline_model_path
        )
        sim_simulate_reads_not_allow_fundamental_ambiguity(
            sim_simpleaf.out.af_quant, 
            sim_err_model.out, 
            spliceu_dir,
            star_index_dir,
            spline_model_path
        )

    emit:
        afa = sim_simulate_reads_allow_fundamental_ambiguity.out.txome_bam
        not_afa = sim_simulate_reads_not_allow_fundamental_ambiguity.out.txome_bam
}