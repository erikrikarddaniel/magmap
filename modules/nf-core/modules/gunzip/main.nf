// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUNZIP {
    tag "$archive"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path archive

    output:
    path "$gunzip",       emit: gunzip
    path "versions.yml" , emit: versions

    script:
    gunzip       = archive.toString() - '.gz'
    """
    gunzip \\
        -f \\
        -c \\
        $options.args \\
        $archive > $gunzip

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
