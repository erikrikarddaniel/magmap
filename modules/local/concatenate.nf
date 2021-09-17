// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CONCATENATE {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pigz=2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0"
    }

    input:
    path files

    output:
    path "file.gz", emit: file

    script:
    cpus = Math.floor(task.cpus/2).toInteger()
    
    // Should the use of options.args for filtering be made explicit by calling it options.filter?
    """
    for f in $files; do unpigz -c -p $cpus \$f; done ${options.args} | pigz -c -p $cpus > file.gz
    """
}
