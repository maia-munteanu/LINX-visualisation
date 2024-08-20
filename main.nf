nextflow.enable.dsl=2

// Define parameters
params.samples_dir = "/g/strcombio/fsupek_home/mmunteanu/LINX/Samples" // Main directory with sample folders
params.samples_list = "/g/strcombio/fsupek_home/mmunteanu/LINX/samples.txt" // File with the list of samples to process
params.ensembl_cache_dir = "/g/strcombio/fsupek_home/mmunteanu/LINX/HMFtools-Resources_dna_pipeline_v5_34_37_hmf_dna_pipeline_resources.37_v5.34/v5_34/ref/37/common/ensembl_data" // Directory with Ensembl cache
params.output_dir = "/g/strcombio/fsupek_home/mmunteanu/LINX/output" // Output directory for plots
params.singularity_image = "/g/strcombio/fsupek_home/mmunteanu/linx_circos_r_container.sif" // Path to the Singularity image

// Workflow definition
workflow {
    main {
        // Channel to read sample names from file
        Channel.fromPath(params.samples_list)
            .splitText()
            .map { it.trim() }
            .filter { it }
            .set { sample_names }

        // Process each sample
        sample_names
            .map { sample_name ->
                file("${params.samples_dir}/${sample_name}")
            }
            .filter { file(it).exists() }
            .set { sample_dirs }

        // Run LINX plotting for each sample
        process_linx_plot(sample_dirs)
    }
}

// Process definition
process process_linx_plot {

    // Specify Singularity container
    container params.singularity_image

    // Input channels
    input:
    path sample_dir from sample_dirs

    // Output channels
    output:
    path "${params.output_dir}/${sample_dir.name}" into plot_output

    // Script to execute
    script:
    """
    mkdir -p ${params.output_dir}/${sample_dir.name}

    java -cp /opt/linx/linx_v1.25.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample ${sample_dir.name} \
        -ensembl_data_dir ${params.ensembl_cache_dir} \
        -ref_genome_version 37 \
        -plot_out ${params.output_dir}/${sample_dir.name}/plot \
        -data_out ${params.output_dir}/${sample_dir.name}/data \
        -vis_file_dir ${params.output_dir}/${sample_dir.name}/ \
        -circos /usr/local/circos-0.69-9/bin/circos \
        -threads 8
    """
}

