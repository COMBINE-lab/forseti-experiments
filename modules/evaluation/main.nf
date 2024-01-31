// in this script, we run the evalution pipeline for the simulated and experimental evalution datasets.

// we use GSE122357 and GSE125188 as the evaluation sets
// These two sets were selected as test sets when fitting the spline model.

// The two simulation datsets

process eval_sim_allow_fundamental_ambiguity {
    publishDir "${params.output_dir}/evaluation", mode: 'symlink'
    input:
        path spliceu_dir
        path mlp_path
        path spline_path
        path txome_bam_path
    output:
        path "${out_dir}"

    script:
        out_dir = "eval_sim_allow_fundamental_ambiguity"
        """
        cat <<EOF > params.json
        {   
            "t2g_path": "${spliceu_dir}/t2g_3col.tsv",
            "txome_path": "${spliceu_dir}/spliceu.fa",
            "mlp_path": "${mlp_path}",
            "spline_path": "${spline_path}",
            "max_frag_len": ${params.eval.max_frag_len},
            "num_threads": ${params.num_threads_large},
            "out_dir": "${out_dir}",
            "txome_bam_path": "${txome_bam_path}"
        }
        EOF

        cp $moduleDir/evaluation_sim.ipynb .
        jupyter nbconvert --to notebook --execute evaluation_sim.ipynb --output ${out_dir}.ipynb
        mv ${out_dir}.ipynb ${out_dir}
        """
}

process eval_sim_not_allow_fundamental_ambiguity {
    publishDir "${params.output_dir}/evaluation", mode: 'symlink'
    input:
        path spliceu_dir
        path mlp_path
        path spline_path
        path txome_bam_path
    output:
        path "${out_dir}"

    script:
    out_dir = "eval_sim_not_allow_fundamental_ambiguity"
    """
    cat <<EOF > params.json
    {   
        "t2g_path": "${spliceu_dir}/t2g_3col.tsv",
        "txome_path": "${spliceu_dir}/spliceu.fa",
        "mlp_path": "${mlp_path}",
        "spline_path": "${spline_path}",
        "max_frag_len": ${params.eval.max_frag_len},
        "num_threads": ${params.num_threads_large},
        "out_dir": "${out_dir}",
        "txome_bam_path": "${txome_bam_path}"
    }
    EOF

    cp $moduleDir/evaluation_sim.ipynb .
    jupyter nbconvert --to notebook --execute evaluation_sim.ipynb --output ${out_dir}.ipynb
    mv ${out_dir}.ipynb ${out_dir}
    """
}

process eval_real_GSE122357 {
    publishDir "${params.output_dir}/evaluation", mode: 'symlink'
    input:
        path spliceu_dir
        path mlp_path
        path spline_path
        path data_dir
    output:
        path "eval_real_${gse}"

    script:
        gse = "GSE122357"
        out_dir = "eval_real_${gse}"
        """
        cat <<EOF > params.json
        {
            "GSE": "${gse}",
            "data_dir": "${data_dir}",
            "spliceu_dir": "${spliceu_dir}",
            "mlp_path": "${mlp_path}",
            "spline_path": "${spline_path}",
            "max_frag_len": ${params.eval.max_frag_len},
            "num_threads": ${params.num_threads_large},
            "out_dir": "${out_dir}"
        }
        EOF

        cp $moduleDir/evaluation_real.ipynb .
        jupyter nbconvert --to notebook --execute evaluation_real.ipynb --output ${out_dir}.ipynb
        mv ${out_dir}.ipynb ${out_dir}
    """
}

process eval_real_GSE125188 {
    publishDir "${params.output_dir}/evaluation", mode: 'symlink'
    input:
        path spliceu_dir
        path mlp_path
        path spline_path
        path data_dir
    output:
        path "eval_real_${gse}"

    script:
        gse = "GSE125188"
        out_dir = "eval_real_${gse}"
        """
        cat <<EOF > params.json
        {
            "GSE": "${gse}",
            "data_dir": "${data_dir}",
            "spliceu_dir": "${spliceu_dir}",
            "mlp_path": "${mlp_path}",
            "spline_path": "${spline_path}",
            "max_frag_len": ${params.eval.max_frag_len},
            "num_threads": ${params.num_threads_large},
            "out_dir": "${out_dir}"
        }
        EOF

        cp $moduleDir/evaluation_real.ipynb .
        jupyter nbconvert --to notebook --execute evaluation_real.ipynb --output ${out_dir}.ipynb
        mv ${out_dir}.ipynb ${out_dir}
        """
}


workflow evaluation {
    main:
        // define the paths
        real_data_dir = file("${params.output_dir}/process_data/definitive_reads_by_r1", checkIfExists: true)
        sim_afa_txome_bam = file("${params.output_dir}/simulation/allow_fundamental_ambiguity/star_out/Aligned.toTranscriptome.out.bam", checkIfExists: true)
        sim_not_afa_txome_bam = file("${params.output_dir}/simulation/not_allow_fundamental_ambiguity/star_out/Aligned.toTranscriptome.out.bam", checkIfExists: true)
        human_spliceu_dir = file("${projectDir}/data/spliceu_ref/refdata-gex-GRCh38-2020-A", checkIfExists: true)
        mouse_spliceu_dir = file("${projectDir}/data/spliceu_ref/refdata-gex-mm10-2020-A", checkIfExists: true)
        mlp = file("${params.output_dir}/models/mlp_model/mlp.pkl", checkIfExists: true)
        spline = file("${params.output_dir}/models/spline_model/spline_model.pkl", checkIfExists: true)

        // run evalution scripts
        // human
        eval_sim_allow_fundamental_ambiguity(
            human_spliceu_dir, 
            mlp, 
            spline, 
            sim_afa_txome_bam
        )
        eval_sim_not_allow_fundamental_ambiguity(
            human_spliceu_dir,
            mlp, 
            spline, 
            sim_not_afa_txome_bam
        )
        eval_real_GSE125188(human_spliceu_dir, mlp, spline, real_data_dir)

        // mouse
        eval_real_GSE122357(mouse_spliceu_dir, mlp, spline, real_data_dir)
}