// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TAXONOMY_GENOME_QUALITY {
    label 'process_m'
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
    path checkm_tsv
    path gtdbtk_arc
    path gtdbtk_bac
    path genes

    output:
    path "genomes.tsv" , emit: genomes
    path "versions.yml", emit: versions

    script:
    
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    setDTthreads($task.cpus)

    # Read the featureCounts files and subset to a set of unique ids
    genomes <- fread(cmd = "gunzip -c ${genes}") %>% lazy_dt() %>%
        distinct(genome) %>%
        as_tibble()

    # CheckM results
    checkm <- fread("${checkm_tsv}", sep = "\\t") %>%
        as_tibble()

    # GTDB-Tk results
    gtdbtk <- fread("${gtdbtk_arc}", sep = "\\t", colClasses = list(character = 1:20)) %>% lazy_dt() %>%
        union(
            fread("${gtdbtk_bac}", sep = "\\t", colClasses = list(character = 1:20)) %>% lazy_dt()
        ) %>%
        as_tibble()

    checkm %>%
        transmute(
            genome = `Bin Id`, completeness = Completeness, contamination = Contamination, strain_heterogeneity = `Strain heterogeneity`
        ) %>%
        full_join(
            gtdbtk %>%
            transmute(genome = user_genome, classification) %>%
            mutate(classification = str_remove_all(classification, '[a-z]__')) %>%
            separate(classification, c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';'),
            by = 'genome'
        ) %>%
        semi_join(genomes, by = 'genome') %>%
        write_tsv("genomes.tsv")

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
