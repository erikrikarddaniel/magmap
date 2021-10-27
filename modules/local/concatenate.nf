// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CONCATENATE {
    label 'process_long'
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
    val  outfile
    path files

    output:
    path "${outfilename}", emit: file

    script:
    cpus    = Math.floor(task.cpus/2).toInteger()

    // If input is gzipped and options.args is set to something, i.e. filtering will be done, unzip input
    catcmd  = ( files[0] =~ /.gz$/ && options.args ) ? "unpigz -c -p $cpus" : "cat"

    // If input is gzipped, but not unzipped by catcmd, don't zip again
    outcmd  = ( files[0] =~ /.gz$/ && catcmd == 'cat' ) ? '' : '| pigz -c -p $cpus'

    // Create a temporary file name to avoid collisions in a second call
    outfilename = ( outfile != '' ) ? outfile : File.createTempFile('outfile', '.gz').getName()
    
    """
    ${catcmd} $files ${options.args} ${outcmd} > $outfilename

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
