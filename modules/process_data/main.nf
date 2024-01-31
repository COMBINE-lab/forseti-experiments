import groovy.json.JsonOutput

// Treat the cDNA insert of the full length scRNA-seq reads as paired end reads and run star aligner to get their alignments. 
// We ask STAR to report the BAM file using both genome and transcriptome coordinates.
// To ease our downstream analysis, we will only keep genomicly unimapped reads that mapped in proper pair.
process run_star {
    publishDir "$params.output_dir/process_data/run_star/$GSE/${SRR}", mode: 'symlink'
    // afterScript "rm -rf ${SRR}_pe ${SRR}_se"
    // species, tissue, chemistry, GSE, SRR, read1, read2, spliceu_dir, index_dir
    input:
        tuple val(species),
            val(tissue),
            val(chemistry), 
            val(GSE),
            val(SRR),
            path("data_dir"),
            path("spliceu_dir"),
            path("index_dir")
    output:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path("data_dir"),
            path("spliceu_dir"),
            path("pe_Aligned.sortedByCoord.out.bam"),
            path("pe_Aligned.toTranscriptome.out.bam"),
            path("se_Aligned.sortedByCoord.out.bam"),
            path("se_Aligned.toTranscriptome.out.bam"), emit: run_star
        path("${SRR}_pe")
        path("${SRR}_se")

    script:
        """
        # we first map paired-end reads
        STAR --runThreadN $params.num_threads \
        --readFilesIn data_dir/paired_end_filtered/${SRR}_1.fastq.gz data_dir/paired_end_filtered/${SRR}_2.fastq.gz \
        --readFilesCommand zcat \
        --genomeDir index_dir \
        --outFileNamePrefix ${SRR}_pe/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMtlen 2 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --alignSplicedMateMapLminOverLmate 0 \
        --quantMode TranscriptomeSAM \
        --quantTranscriptomeBan Singleend


        mv ${SRR}_pe/Aligned.toTranscriptome.out.bam pe_Aligned.toTranscriptome.out.bam
        samtools view --require-flags 2 -d NH:1 ${SRR}_pe/Aligned.sortedByCoord.out.bam -o pe_Aligned.sortedByCoord.out.bam
        
        # we then map read2 as single-end reads
        STAR --runThreadN $params.num_threads \
        --readFilesIn data_dir/paired_end_filtered/${SRR}_2.fastq.gz \
        --readFilesCommand zcat \
        --genomeDir index_dir \
        --outFileNamePrefix ${SRR}_se/ \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --quantMode TranscriptomeSAM \
        --quantTranscriptomeBan Singleend

        mv ${SRR}_se/Aligned.toTranscriptome.out.bam se_Aligned.toTranscriptome.out.bam
        samtools view -d NH:1 ${SRR}_se/Aligned.sortedByCoord.out.bam -o se_Aligned.sortedByCoord.out.bam

        """
    
    stub:
        """
        touch se_Aligned.sortedByCoord.out.bam
        touch se_Aligned.toTranscriptome.out.bam
        touch pe_Aligned.sortedByCoord.out.bam
        touch pe_Aligned.toTranscriptome.out.bam
        
        """
}

process frag_len_dist {
    publishDir "$params.output_dir/process_data/frag_len_dist/${GSE}", mode: 'symlink'
    // afterScript "find $workDir -name '${reads1}' -type f -exec rm -rf {} + && find $workDir -name '${reads2}' -type f -exec rm -rf {} +"

    input:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path("data_dir"),
            path("spliceu_dir"),
            path(pe_genome_bam),
            path(pe_txome_bam), // not used
            path(se_genome_bam), // not used
            path(se_txome_bam) // not used
    output:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path("$SRR")

    script:
        outdir="$SRR"

        """
        # just for rerun and rerun
        mkdir -p $outdir

        cp $moduleDir/correct_tlen_plot.py .
        python correct_tlen_plot.py ${pe_genome_bam} ${GSE} ${SRR} data_dir/paired_end_filtered/${SRR}_fastp.json ${chemistry} spliceu_dir/genome.fa ${outdir}
        """
    
    stub:
        """
        mkdir -p $SRR
        touch $SRR/SRR6494015_new_tlens.txt
        touch $SRR/polyt_seq.txt
        touch $SRR/bg_seq.txt
        touch $SRR/polya_rc_seq.txt
        touch $SRR/polyt_bg_seq.txt
        touch $SRR/polya_bg_seq.txt
        """ 
}

process write_json {
    label 'local'
    publishDir "$params.output_dir/process_data/frag_len_dist/${GSE}/${SRR}", mode: 'symlink'

    input:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path(frag_len_dist)

    output:
        path "${SRR}.json", emit: dataset_description
        
    exec:
        // from here
        // https://groups.google.com/g/nextflow/c/tp_b1p0DBE4?pli=1
        def json = JsonOutput.toJson(  [species: "${species}", 
                                        tissue: "${tissue}", 
                                        chemistry: "${chemistry}", 
                                        GSE: "${GSE}", 
                                        SRR: "${SRR}"
                                        ]
                                    )
        task.workDir.resolve("${SRR}.json").text = json
}

process definitive_reads_by_r1 {
    publishDir "$params.output_dir/process_data/definitive_reads_by_r1/${GSE}", mode: 'symlink'
    afterScript "rm -rf r1_r2_exonic.bed"
    input:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path(data_dir),
            path(spliceu_dir),
            path(pe_genome_bam),
            path(pe_txome_bam), // not used
            path(se_genome_bam), // not used
            path(se_txome_bam) // not used
    output:
        tuple val(species),
            val(tissue),
            val(chemistry),
            val(GSE),
            val(SRR),
            path("$SRR")

    script:
        """
        bedtools intersect \
        -a $pe_genome_bam \
        -b $spliceu_dir/exon_by_tx.bed \
        -f 1.0 \
        -wo \
        -bed \
        > r1_r2_exonic.bed

        cp $moduleDir/classify_r1.py .

        python classify_r1.py ${GSE} ${SRR} $se_txome_bam $pe_txome_bam r1_r2_exonic.bed $spliceu_dir/t2g_3col.tsv ${params.exon_dist_thresh} ${SRR}

        """
    
    stub:
        """
        mkdir -p $SRR
        touch $SRR/r2_exonix_r1_spliced.bam
        touch $SRR/r2_exonix_r1_unspliced.bam
        """
}

workflow process_data {
    sample_sheet = Channel
            .fromPath(params.input_files.sample_sheet)
            .splitCsv(header:true, sep:",", strip: true)
            .map{ row-> tuple(row.Species,
                                row.Tissue,
                                row.Chromium_version,
                                row.GSE,
                                row.SRR,
                                "${params.data_dir}/datasets/${row.GSE}/${row.SRR}")
            }
            
    ref_sheet = Channel
            .fromPath(params.input_files.ref_sheet)
            .splitCsv(header:true, sep:",", strip: true)
            .map{ row-> tuple(row.species,
                                "${params.data_dir}/spliceu_ref/${row.ref_name}",
                                "${params.data_dir}/star_spliceu_index/${row.ref_name}"
                                )
            }

    run_star_input = sample_sheet.combine(ref_sheet, by: 0)
    // species, tissue, chemistry, GSE, SRR, data_dir, spliceu_dir, index_dir
    run_star(run_star_input)
    frag_len_dist(run_star.out.run_star)
    definitive_reads_by_r1(run_star.out.run_star)
    write_json(frag_len_dist.out)
}
