import groovy.json.JsonOutput

include {prepare_data} from "./modules/prepare_data"
include {process_data} from "./modules/process_data"
include {fit_models} from "./modules/fit_models"
include {simulation} from "./modules/simulation"
include {evaluation} from "./modules/evaluation"

workflow {
    if (params.mode == null) {
        throw new RuntimeException("Please provide a run mode. Either 'full' or 'process'")
    } else if (params.mode == 'full') {
        println "Running full pipeline"
        prepare_data().collect()
        process_data().collect()
        fit_model().collect()
        simulation().collect()
        evaluation.collect()
    } else if (params.mode == 'prepare_data') {
        println "Preparing data"
        prepare_data()
    } else if (params.mode == 'process_data') {
        println "Processing data"
        process_data()
    } else if (params.mode == 'fit_models') {
        println "Fitting models"
        fit_models()
    } else if (params.mode == 'simulation') {
        println "Running simulation"
        simulation()
    } else if (params.mode == 'evaluation') {
        println "Evaluating models"
        evaluation()
    } else {
        throw new RuntimeException("The provided run mode, ${params.mode} is invalid; Cannot proceed.")
    }

}
 
//  workflow fit_model {
//     fit_splines(params.input_files.sample_sheet)
//     fit_mlp()
//  }

//  process fit_spline {
//     publishDir "$params.output_dir/models", mode: 'symlink'
//     input:
//         path(sample_sheet)
    
//     output:
//         path("spline_model")
//     """
//     $projectDir/scripts/fit_spline.py ${params.output_dir} $sample_sheet spline_model
//     """

//  }

//  process fit_mlp {
//     publishDir "$params.output_dir/models", mode: 'symlink'
//     input:
//         path(sample_sheet)
    
//     output:
//         path("mlp_model")
//     """
//     $projectDir/scripts/fit_mlp.py $sample_sheet mlp_model
//     """

//  }