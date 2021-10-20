// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FETCH_NCBI {
    tag "$ncbi_accs"
    label 'process_long'
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
    path "*.fna.gz"    , emit: fnas          // Contigs
    path "failures.txt", emit: failures
    path "versions.yml", emit: versions

    script:
    
    """
    mkdir genomes/

    # Fetch file indexes
    wget -O 00refseq.index  ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
    wget -O 10genbank.index ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

    echo "Contig files for the following accessions were not found at NCBI" > failures.txt
    for a in \$(cat $ncbi_accs); do
        d=\$(grep \$a *.index | cut -f 20 | head -n 1 | sed 's/https/ftp/')
        echo "--> \$a: \$d <--"
        \$(wget \${d}/*_genomic.fna.gz)
        rm -f *from_genomic*.fna.gz 
        if [ \$(ls *_genomic.fna.gz | wc -l) ]; then
            echo "\$a"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( wget --version | grep '^GNU' | sed 's/GNU Wget //' | sed 's/ .*//' )
    END_VERSIONS
    """
}
