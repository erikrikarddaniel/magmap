// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FETCH_NCBI {
    tag "$ncbi_accs"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "anaconda::wget=1.20.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/wget:1.20.1"
    } else {
        container "quay.io/biocontainers/wget:1.20.1"
    }

    input:
    path ncbi_accs

    output:
    path "*.fna.gz"    , emit: fnas         // Contigs
    path "*.gff.gz"    , emit: gffs
    path "*.faa.gz"    , emit: faas
    path "failures.txt"        , emit: failures
    path "versions.yml"        , emit: versions

    script:
    
    """
    mkdir genomes/

    # Fetch file indexes
    wget -O 00refseq.index  ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
    wget -O 10genbank.index ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

    echo "The following files were not found at NCBI" > failures.txt
    for a in \$(cat $ncbi_accs); do
        echo "--> \$a" >> log
        d=\$(grep \$a *.index | cut -f 20 | head -n 1)
        for s in gff.gz _genomic.fna.gz faa.gz; do
            echo "      \$s" >> log
            \$( wget \${d}/*\${s} )
            if [ \$(ls \${a}*\${s} | wc -l) -eq 0 ]; then
                echo "Found no \${d}/*\${s}" >> failures.txt
            fi
        done
        rm -f *from_genomic*.fna.gz 
    done

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( wget --version | grep '^GNU' | sed 's/GNU Wget //' | sed 's/ .*//' )
    END_VERSIONS
    """
}
