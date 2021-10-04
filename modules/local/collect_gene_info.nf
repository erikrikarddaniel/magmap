// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COLLECT_GENE_INFO {
    tag 'genes.tsv.gz'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0"
    }

    input:
    path files

    output:
    path "genes.tsv.gz"          , emit: genefile
    path "versions.yml"          , emit: versions

    script:
    def software = getSoftwareName(task.process)
    
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    setDTthreads($task.cpus)


    # Read all gff files and separate out the last, ';'-separated field into separate columns
    tibble(f = Sys.glob('*.gff.gz')) %>%
        mutate(
            data = purrr::map(
                f, 
                function(file) {
                    fread(
                        cmd = sprintf("gunzip -c %s | grep -E '\t'", file), sep = '\t',
                        col.names = c('chromosome', 'gcaller', 'type', 'from', 'to', 'b', 'strand', 'c', 'rest') 
                    ) %>%
                    select(-b, -c) %>%
                    as_tibble() %>%
                    separate_rows(rest, sep = ';') %>%
                    separate(rest, c('t', 'v'), sep = '=') %>%
                    pivot_wider(names_from = t, values_from = v)
                }
            )
        ) %>%
        mutate(genome = str_remove(f, ".gff.gz")) %>%
        select(-f) %>%
        relocate(genome) %>%
        unnest(data) %>%
        write_tsv("genes.tsv.gz")

    write(
        sprintf(
            "${getProcessName(task.process)}:\n    R: %s.%s\n    dplyr: %s\n    dtplyr: %s\n    tidyr: %s\n    data.table: %s\n",
            R.Version()\$major, R.Version()\$minor,
            packageVersion('dplyr'),
            packageVersion('dtplyr'),
            packageVersion('tidyr'),
            packageVersion('data.table')
        ),
        'versions.yml'
    )
    """
}
