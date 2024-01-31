import groovy.json.JsonOutput


// build spliced + unspliced reference for each transcript. Here each isoform of each gene corresponds to two reference sequences, one is the spliced version, the concatenated exons, and the other is the unspliced version, the transcript body. The unspliced version of transcript has a single exon containing the transcript body.
// Meanwhile, we will also find the polyA/T sites for each transcript, and output a bed file containing the polyA/T sites of transcripts. 
// A polyA/T site is a polyA/T k-mer with at most m mismatches to the reference genome. We also require that the mismatches are not at the starting and ending base of the k-mer.
// Notice that the polyA tail of each transcript is not included in the bed file. The position of polyA tail corresponds to the transcript length in the SAM/BAM file header.

process spliceu_ref {
    label 'r'
    publishDir "$params.data_dir/spliceu_ref", mode: 'move'

    input:
        tuple val(species),
            val(ref_name),
            path(genome_path),
            path(gtf_path)
    
    output:
        tuple val(species),
            val(ref_name), emit: spliceu_ref
        path "${ref_name}"

    script:
        """
        cp $moduleDir/build_spliceu.R .
        Rscript build_spliceu.R ${genome_path} ${gtf_path} ${ref_name} ${params.num_threads} ${params.polya_min_length} ${params.polya_max_mismatch}
        mv $genome_path ${ref_name}/genome.fa
        mv $gtf_path ${ref_name}/genes.gtf
        """

    stub:
        """
        mkdir ${ref_name}
        touch ${ref_name}/spliceu.gtf
        touch ${ref_name}/spliceu.fa
        touch ${ref_name}/exon_by_tx.bed
        touch ${ref_name}/tx_polya.bed
        touch ${ref_name}/tx_polyt.bed
        touch ${ref_name}/t2g_3col.tsv
        touch ${ref_name}/genome.fa
        touch ${ref_name}/genes.gtf
        """
}

// This process builds the start index using the spliceu gene annotations. By utilizing this reference, we will be able to find the compatible spliced and unspliced transcripts for each read. 
// I have confirmed that spliced reads will not show compatibility with unspliced transcripts :)

process star_index_spliceu {
    publishDir "${params.data_dir}/star_spliceu_index", mode: 'move'

    input:
        tuple val(species),
            val(ref_name)
    
    output:
        tuple val(species),
            val(ref_name),
            path("${ref_name}")

    script:
        genome_path = file("${params.data_dir}/spliceu_ref/${ref_name}/genome.fa")
        spliceu_gtf_path = file("${params.data_dir}/spliceu_ref/${ref_name}/spliceu.gtf") 
        """
        STAR --runThreadN $params.num_threads \
        --runMode genomeGenerate --genomeDir ${ref_name} \
        --genomeFastaFiles $genome_path \
        --sjdbGTFfile $spliceu_gtf_path \
        --sjdbOverhang 100
        """
    
    stub:
        """
        mkdir ${ref_name}
        """
}

 process split_reads() {
    input:
        tuple val(species),
            val(tissue),
            val(chemistry), 
            val(GSE),
            val(SRR),
            path(r1_path),
            path(r2_path)
    
    output:
        tuple val(species),
            val(tissue),
            val(chemistry), 
            val(GSE),
            val(SRR),
            path("${SRR}_1.fastq.gz"),
            path("${SRR}_2.fastq.gz"),
            path("${SRR}_1_pe.fastq.gz"), emit: split_reads

    script:
        """
        cp $moduleDir/split_read1.py .

        python split_read1.py ${SRR} ${SRR}_1.fastq.gz $chemistry ${SRR}
        rm -f \$(readlink -f $r1_path)
        mv \$(readlink -f $r2_path) ${SRR}_2.fastq.gz
        mv ${SRR}/${SRR}_1.fastq.gz ${SRR}_1.fastq.gz
        mv ${SRR}/${SRR}_1_pe.fastq.gz ${SRR}_1_pe.fastq.gz
        """
    
    stub:
        """
        mkdir ${SRR}
        touch ${SRR}/${SRR}_1.fastq.gz
        """
}

// In this step, we filter low quality reads. We disabled adapter trimming becuase the biological part of read 1 does not at the beginning of read1.
process fastp {
    publishDir "$projectDir/data/datasets/${GSE}", mode: 'move'
    // afterScript "find ${workDir}/stage* -name '${r2_path}' -type f -exec rm -rf {} +"

    // we do not need the pe file anymore
    input:
        tuple val(species),
            val(tissue),
            val(chemistry), 
            val(GSE),
            val(SRR),
            path(cbumi),
            path(r2_path),
            path(r1_path)
    output:
        path("${SRR}")

    script:
        """
        mkdir -p ${SRR}/fastp_filtered

        fastp -i $r1_path -I $r2_path \
        -o ${SRR}/fastp_filtered/${SRR}_1.fastq.gz -O ${SRR}/fastp_filtered/${SRR}_2.fastq.gz \
        --qualified_quality_phred $params.fastp.qualified_quality_phred \
        --unqualified_percent_limit $params.fastp.unqualified_percent_limit \
        --low_complexity_filter --complexity_threshold $params.fastp.complexity_threshold \
        -h ${SRR}/paired_end_filtered/${SRR}_fastp.html \
        -j ${SRR}/paired_end_filtered/${SRR}_fastp.json \
        --disable_adapter_trimming \
        --length_required 10

        rm -f \$(readlink -f $r1_path)

        mv \$(readlink -f $r2_path) ${SRR}/${SRR}_2.fastq.gz
        mv \$(readlink -f $cbumi) ${SRR}/${SRR}_1.fastq.gz
        """

    stub:
        """
        mkdir -p $SRR
        touch ${SRR}/${SRR}_fastp.json
        touch ${SRR}/${SRR}_1.fastq.gz
        touch ${SRR}/${SRR}_2.fastq.gz
        """
}

workflow prepare_data {
    sample_sheet = Channel
            .fromPath(params.input_files.sample_sheet)
            .splitCsv(header:true, sep:",", strip: true)
            .map{ row-> tuple(row.Species,
                                row.Tissue,
                                row.Chromium_version,
                                row.GSE,
                                row.SRR,
                                row.URL_R1,
                                row.URL_R2)
            }

    ref_sheet = Channel
            .fromPath(params.input_files.ref_sheet)
            .splitCsv(header:true, sep:",", strip: true)
            .map{ row-> tuple(row.species,
                                row.ref_name,
                                "${projectDir}/${row.genome_path}",
                                "${projectDir}/${row.gtf_path}")
            }
    spliceu_ref(ref_sheet)
    star_index_spliceu(spliceu_ref.out.spliceu_ref)
    // split_reads(sample_sheet)
    // fastp(split_reads)
}
 
