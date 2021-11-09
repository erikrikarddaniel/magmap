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
    val  samples
    path trimlogs
    path idxstats
    path fcs
    path bbduks

    output:
    path "overall_stats.tsv"     , emit: overall_stats
    path "versions.yml"          , emit: versions

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

    TYPE_ORDER = c('n_pairs_trimmed', 'n_pairs_non_contaminated', 'n_pairs_mapped', 'n_pairs_unmapped', 'n_pairs_features', 'n_pairs_cds')

    # Collect stats for each sample, create a table in long format that can be appended to
    t <- tibble(sample = c("${samples.join('", "')}")) %>%
        mutate(
            # N. after trimming
            t = map(
                sample, 
                function(s) {
                    fread(
                        cmd = sprintf("grep 'Reads written (passing filters)' %s*trimming_report.txt | sed 's/.*: *//' | sed 's/ .*//' | sed 's/,//g'", s), 
                        sep = ',',
                        col.names = c('n_pairs_trimmed')
                    )
                }
            ),
            i = map(
                sample,
                function(s) {
                    fread(cmd = sprintf("grep -v '^*' %s*idxstats", s), sep = '\\t', col.names = c('chr', 'length', 'n_pairs_mapped', 'n_pairs_unmapped')) %>%
                        lazy_dt() %>%
                        summarise(n_pairs_mapped = sum(n_pairs_mapped)/2, n_pairs_unmapped = sum(n_pairs_unmapped)/2) %>%
                        as_tibble()
                }
            )
        ) %>%
        unnest(c(t, i)) %>%
        pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v') %>%
        union(
            # Total observation after featureCounts
            tibble(file = Sys.glob('counts*.tsv.gz')) %>%
                mutate(d = map(file, function(f) fread(cmd = sprintf("gunzip -c %s", f), sep = '\\t'))) %>%
                as_tibble() %>%
                unnest(d) %>%
                group_by(sample) %>% summarise(n_pairs_features = sum(count), .groups = 'drop') %>%
                pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v')
        ) %>%
        union(
            # Total observation after featureCounts
            fread(cmd = "gunzip -c counts.CDS.tsv.gz", sep = '\\t') %>%
                as_tibble() %>%
                group_by(sample) %>% summarise(n_pairs_cds = sum(count), .groups = 'drop') %>%
                pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v')
        )

    # Add in stats from BBDuk, if present
    for ( f in Sys.glob('*.bbduk.log') ) {
        s = str_remove(f, '.bbduk.log')
        t <- t %>% union(
            fread(cmd = sprintf("grep 'Result:' %s | sed 's/Result:[ \\t]*//; s/ reads.*//'", f), col.names = c('v')) %>%
                as_tibble() %>%
                mutate(sample = s, m = 'n_pairs_non_contaminated', v = v/2)
        )
    }

    # Write the table in wide format
    t %>% 
        mutate(m = parse_factor(m, levels = TYPE_ORDER, ordered = TRUE)) %>%
        arrange(sample, m) %>%
        pivot_wider(names_from = m, values_from = v) %>%
        write_tsv('overall_stats.tsv')

    write(
        sprintf(
            "${getProcessName(task.process)}:\n    R: %s.%s\n    dplyr: %s\n    dtplyr: %s\n    data.table: %s\n    readr: %s\n    purrr: %s\n    tidyr: %s\n    stringr: %s\n",
            R.Version()\$major, R.Version()\$minor,
            packageVersion('dplyr'),
            packageVersion('dtplyr'),
            packageVersion('data.table'),
            packageVersion('readr'),
            packageVersion('purrr'),
            packageVersion('tidyr'),
            packageVersion('stringr')
        ),
        'versions.yml'
    )
    """
}
