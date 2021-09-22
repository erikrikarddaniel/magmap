// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COLLECTDATA {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1"
    } else {
        container "quay.io/biocontainers/r-tidyverse:1.2.1"
    }

    input:
    //tuple val(samples), path(fc_tables)
    path('*.featureCounts.txt')

    output:
    path "counts.tsv.gz", emit: counts
    //path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    //println "DEBUG: fc_tables: ${fc_tables}"
    
    """
    #!/usr/bin/env Rscript

    #library(data.table)
    #library(dtplyr)
    library(readr)
    library(dplyr)
    library(stringr)


    tibble(f = Sys.glob('*.featureCounts.txt')) %>%
        mutate(
            d = purrr::map(
                f, 
                function(file) {
                    #read_tsv(file, col_types = 'cciicii', skip = 1) %>%
                    read.delim(file, skip = 1, sep = '\t') %>%
                        tidyr::pivot_longer(ncol(.), names_to = 'sample', values_to = 'count') %>%
                        filter(count > 0) %>%
                        mutate(sample = str_remove(sample, '.sorted.bam'))
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        mutate(r = count/Length) %>%
        group_by(sample) %>% mutate(tpm = r/sum(r) * 1e6) %>% ungroup() %>%
        select(-f, -r) %>%
        write_tsv("counts.tsv.gz")
    """
}
