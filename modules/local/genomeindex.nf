// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENOMEINDEX {
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
    path gffs

    output:
    path "${outfilename}", emit: genomes2id
    path "versions.yml"  , emit: versions

    script:
    outfilename = File.createTempFile('outfile', '.gz').getName()
    cpus        = Math.floor(task.cpus/2).toInteger()
    
    """
    echo "genome\tID" | gzip -c > ${outfilename}
    for f in ${gffs}; do
        unpigz -c -p ${cpus} \$f | grep -o 'ID=[A-Z0-9_]\\+' | sed "s/^/\$f\\t/; s/ID=//; s/.gff.gz//" | pigz -c -p ${cpus} >> ${outfilename}
    done

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
