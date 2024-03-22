import groovy.json.JsonOutput

 process fit_spline {
    publishDir "$params.output_dir/models", mode: 'symlink', pattern: 'spline_model'
    input:
        path parent_dir
        path sample_sheet
    
    output:
        path "spline_model"
        path "spline_model/spline_model.pkl", emit: spline
    """
    $moduleDir/fit_spline.py ${params.output_dir} $sample_sheet spline_model
    """

 }

 process fit_spline_ipynb {
    publishDir "$params.output_dir/models", mode: 'symlink', pattern: 'spline_model'
    input:
        path parent_dir
        path sample_sheet
    
    output:
        path "spline_model"
        path "spline_model/spline_model.pkl", emit: spline
    """
    cat <<EOF > params.json
    {   
        "parent_dir": "${parent_dir}",
        "PE_sheet": "${sample_sheet}",
        "outdir": "spline_model",
        "random_seed": ${params.spline.random_seed}
    }
    EOF
    cp $moduleDir/fit_spline.ipynb .
    jupyter nbconvert  --execute --to html fit_spline.ipynb --output fit_spline.ipynb
    mv fit_spline.ipynb spline_model
    """
 }


 process fit_mlp {
    publishDir "$params.output_dir/models", mode: 'symlink'
    input:
        path(sample_sheet)
    
    output:
        path("mlp_model")
    """
    $moduleDir/fit_mlp.py $sample_sheet mlp_model
    """

 }
 process fit_mlp_ipynb {
    publishDir "$params.output_dir/models", mode: 'symlink', pattern: 'mlp_model'
    input:
        path parent_dir
        path sample_sheet
    
    output:
        path("mlp_model")
        path("mlp_model/mlp.pkl"), emit: mlp

    """
    cat <<EOF > params.json
    {   
        "parent_dir": "${parent_dir}",
        "PE_sheet": "${sample_sheet}",
        "outdir": "mlp_model",
        "random_seed": ${params.mlp.random_seed},
        "snr_len": ${params.mlp.snr_len},
        "snr_mismatch": ${params.mlp.snr_mismatch}
    }
    EOF

    cp $moduleDir/fit_mlp.ipynb .
    jupyter nbconvert  --execute --to notebook fit_mlp.ipynb --output fit_mlp.ipynb
    mv fit_mlp.ipynb mlp_model
    """
 }

workflow fit_models {
    main:
    // fit_spline(params.input_files.sample_sheet)
    // fit_mlp(params.input_files.sample_sheet)
    fit_spline_ipynb(params.output_dir, params.input_files.sample_sheet)
    fit_mlp_ipynb(params.output_dir, params.input_files.sample_sheet)

    emit:
    spline = fit_spline_ipynb.out.spline
    mlp = fit_mlp_ipynb.out.mlp
}
