nextflow.enable.dsl = 2

workflow {
    samples = Channel.fromPath(params.samples_file, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row -> [ row.sample, file(row.input) ] }
    
    linx_circos_plot(samples)
}

process linx_circos_plot {
    tag { sample }

    input:
    tuple val(sample), path(input)

    output:
    path "plot"
    path "data"

    publishDir "${params.output_directory}/${sample}", mode: 'move', pattern: 'plot'
    publishDir "${params.output_directory}/${sample}", mode: 'move', pattern: 'data'

    script:
    """
    mkdir -p tmp
    export TMPDIR=\$(pwd)/tmp
    java -cp /opt/linx/linx_v1.25.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample ${sample} \
        -vis_file_dir ${input} \
        -circos /usr/local/circos-0.69-9/bin/circos \
        -ensembl_data_dir ${params.ensembl_cache_dir} \
        -ref_genome_version ${params.ref} \
        -plot_out ./plot \
        -data_out ./data \
        -threads ${params.threads}
    rm -rf tmp  
    """
}

