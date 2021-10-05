// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COLLECT_STATS {
    label 'process_low'
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
    val samples
    path trimlogs
    path idxstats

    output:
    path "overall_stats.tsv"     , emit: overall_stats
    //path "versions.yml"          , emit: versions

    script:
    
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    tibble(sample = c("${samples.join('", "')}")) %>%
        mutate(
            t = map(
                sample, 
                function(s) {
                    fread(
                        cmd = sprintf("grep 'Reads written (passing filters)' %s*trimming_report.txt | sed 's/.*: *//' | sed 's/ .*//'", s), 
                        sep = ',',
                        col.names = c('n_trimmed')
                    )
                }
            ),
            i = map(
                sample,
                function(s) {
                    fread(cmd = sprintf("grep -v '^*' %s*idxstats", s), sep = '\\t', col.names = c('chr', 'length', 'n_mapped', 'n_unmapped')) %>%
                        lazy_dt() %>%
                        summarise(n_mapped = sum(n_mapped), n_unmapped = sum(n_unmapped)) %>%
                        as_tibble()
                }
            )
        ) %>%
        unnest(t, i) %>%
        write_tsv('overall_stats.tsv')
    """
}
/**
    #!/usr/bin/env Rscript

    library(data.table)
    #library(dtplyr)
    library(readr)
    library(dplyr)
    library(stringr)

    setDTthreads($task.cpus)

    tibble(f = Sys.glob('*_trimming_report.txt')) %>%
        mutate(sample = str_remove(f, '.fastq.*')) %>%
        write_tsv('overall_stats.tsv')

    write(
        sprintf(
            "${getProcessName(task.process)}:\n    R: %s.%s\n    dplyr: %s\n    dtplyr: %s\n    data.table: %s\n",
            R.Version()\$major, R.Version()\$minor,
            packageVersion('dplyr'),
            packageVersion('dtplyr'),
            packageVersion('data.table')
        ),
        'versions.yml'
    )
    **/
