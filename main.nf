nextflow.enable.dsl=2

// Parameters
params.samples_file = "/g/strcombio/fsupek_home/mmunteanu/LINX/samples.txt" // File with the list of samples to process
params.output_directory = "/g/strcombio/fsupek_home/mmunteanu/LINX/output"
params.ensembl_cache_dir = "/g/strcombio/fsupek_home/mmunteanu/LINX/HMFtools-Resources_dna_pipeline_v5_34_37_hmf_dna_pipeline_resources.37_v5.34/v5_34/ref/37/common/ensembl_data" // Directory with Ensembl cache
params.singularity_image = "/g/strcombio/fsupek_home/mmunteanu/linx_circos_r_container.sif" // Path to the Singularity image
params.ref = "37"
params.threads = "8"

// Workflow
workflow {
    samples = Channel.fromPath(params.samples_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true).map { row -> [ row.sample, path(row.input) ] } 
    linx_circos_plot(samples)
}

// Process
process linx_circos_plot {
    tag { sample }
    container params.singularity_image

    input:
    tuple val(sample), path(input)

    output:
    path "plot" into plot_output
    path "data" into data_output

    publishDir "${params.output_directory}/${sample}", mode: 'move', pattern: 'plot'
    publishDir "${params.output_directory}/${sample}", mode: 'move', pattern: 'data'


    script:
    """
    java -cp /opt/linx/linx_v1.25.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample ${sample} \
        -vis_file_dir ${input} \
        -circos /usr/local/circos-0.69-9/bin/circos \
        -ensembl_data_dir ${params.ensembl_cache_dir} \
        -ref_genome_version ${params.ref} \
        -plot_out ./plot \
        -data_out ./data \
        -threads ${params.threads}
    """
}

