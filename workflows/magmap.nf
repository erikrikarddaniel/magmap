/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMagmap.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Populate the ch_genome_fnas channel either from the params.genome_fna_dir and params.genome_fna_ext combination or the genome_fnas list
if ( params.genome_fna_dir && params.genome_fna_ext ) {
    Channel
        .fromPath(params.genome_fna_dir + params.genome_fna_ext)
        .ifEmpty { exit 1, "Cannot find any files matching \"${params.genome_fna_dir}${params.genome_fna_ext}\"\nPlease revise the genome directory (\"--genome_fna_dir\") and extension (\"--genome_fna_ext\")." }
        .set { ch_genome_fnas }
} else if ( params.genome_fnas ) {
    Channel
        .fromPath(Arrays.asList(params.genome_fnas.split(',')))
        .set { ch_genome_fnas }
} else if ( params.ncbi_accessions ) {
    Channel
        .fromPath( params.ncbi_accessions )
        .set {  ch_ncbi_accessions }
} else {
    exit 1, "You need to specifiy either a combination of genome fna directory and genome fna extension (\"--genome_fna_dir\" and \"--genome_fna_ext\") or a comma separated list of genome fna files (\"--genome_fnas\")."
}

// Populate the ch_genome_gffs channel either from the params.genome_gff_dir and params.genome_gff_ext combination or the genome_gffs list
if ( params.genome_gff_dir && params.genome_gff_ext ) {
    Channel
        .fromPath(params.genome_gff_dir + params.genome_gff_ext)
        .ifEmpty { exit 1, "Cannot find any files matching \"${params.genome_gff_dir}${params.genome_gff_ext}\"\nPlease revise the genome directory (\"--genome_gff_dir\") and extension (\"--genome_gff_ext\")." }
        .set { ch_genome_gffs }
} else if ( params.genome_gffs ) {
    Channel
        .fromPath(Arrays.asList(params.genome_gffs.split(',')))
        .set { ch_genome_gffs }
} else if ( ! params.ncbi_accessions ) {
    exit 1, "You need to specifiy either a combination of genome gff directory and genome gff extension (\"--genome_gff_dir\" and \"--genome_gff_ext\") or a comma separated list of genome gff files (\"--genome_gffs\")."
}

// If params.ncbi_accessions is set, turn off counting of non-CDS features
count_cds = params.count_cds
if ( ! params.ncbi_accessions ) {
    count_rrna  = params.count_rrna
    count_trna  = params.count_trna
    count_tmrna = params.count_tmrna
} else {
    count_rrna  = false
    count_trna  = false
    count_tmrna = false
}

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

collect_featurecounts_options_cds   = modules['collect_featurecounts_cds']
collect_featurecounts_options_rrna  = modules['collect_featurecounts_rrna']
collect_featurecounts_options_trna  = modules['collect_featurecounts_trna']
collect_featurecounts_options_tmrna = modules['collect_featurecounts_tmrna']

collect_stats_options               = modules['collect_stats']
collect_gene_info_options           = modules['collect_gene_info']

//
// MODULE: Local to the pipeline
//
include { COLLECT_STATS                                        } from '../modules/local/collect_stats.nf'      addParams( options: collect_stats_options )
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_CDS   } from '../modules/local/collect_featurecounts' addParams( options: collect_featurecounts_options_cds )
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_RRNA  } from '../modules/local/collect_featurecounts' addParams( options: collect_featurecounts_options_rrna )
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TRNA  } from '../modules/local/collect_featurecounts' addParams( options: collect_featurecounts_options_trna )
include { COLLECT_FEATURECOUNTS as COLLECT_FEATURECOUNTS_TMRNA } from '../modules/local/collect_featurecounts' addParams( options: collect_featurecounts_options_tmrna )

include { COLLECT_GENE_INFO     } from '../modules/local/collect_gene_info' addParams( options: collect_gene_info_options )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
fetch_ncbi_options                  = modules['fetch_ncbi']
prodigal_options                    = modules['prodigal']

bbmap_index_options                 = modules['bbmap_index']
bbmap_index_options.args           += params.usemodulo ? Utils.joinModuleArgs(["usemodulo"]) : ''

include { INPUT_CHECK         } from '../subworkflows/local/input_check'         addParams(options: [:])
include { FETCH_NCBI_PRODIGAL } from '../subworkflows/local/fetch_ncbi_prodigal' addParams(fetch_ncbi_options: fetch_ncbi_options, prodigal_options: prodigal_options)
include { CREATE_BBMAP_INDEX  } from '../subworkflows/local/create_bbmap_index'  addParams(bbmap_index_options: bbmap_index_options)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options       = modules['multiqc']
multiqc_options.args     += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

// Move most of these (all except bbmap_align) to addParams() call directly?
bbduk_options                            = modules['bbduk']
bbmap_align_options                      = modules['bbmap_align']
bbmap_align_options.args                += params.usemodulo ? Utils.joinModuleArgs(["usemodulo"]) : ''
concatenate_gff_options                  = modules['concatenate_gff']
subread_featurecounts_options_cds        = modules['subread_featurecounts_cds']
subread_featurecounts_options_rrna       = modules['subread_featurecounts_rrna']
subread_featurecounts_options_trna       = modules['subread_featurecounts_trna']
subread_featurecounts_options_tmrna      = modules['subread_featurecounts_tmrna']

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC            } from '../modules/nf-core/modules/fastqc/main'               addParams( options: modules['fastqc'] )
include { MULTIQC           } from '../modules/nf-core/modules/multiqc/main'              addParams( options: multiqc_options   )
include { BBMAP_BBDUK       } from '../modules/nf-core/modules/bbmap/bbduk/main'          addParams( options: bbduk_options )
include { BBMAP_ALIGN       } from '../modules/nf-core/modules/bbmap/align/main'          addParams( options: bbmap_align_options )
include { BAM_SORT_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_samtools'            addParams( sort_options: modules['samtools_sort_genomes'], index_options: modules['samtools_index_genomes'], stats_options: modules['samtools_index_genomes']      )
include { CONCATENATE as CONCATENATE_GFF } from '../modules/local/concatenate'            addParams( options: concatenate_gff_options )
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS   } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options_cds )
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_RRNA  } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options_rrna )
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TRNA  } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options_trna )
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TMRNA } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options_tmrna )
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_files : ['_versions.yml':'']] )

//
// SUBWORKFLOW: Adapted from rnaseq!
//
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''

include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MAGMAP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: if asked to, fetch genomes from NCBI, then call ORFs with Prodigal
    //
    if ( params.ncbi_accessions ) {
        FETCH_NCBI_PRODIGAL(ch_ncbi_accessions)
        ch_genome_fnas = FETCH_NCBI_PRODIGAL.out.fnas
        ch_genome_gffs = FETCH_NCBI_PRODIGAL.out.gffs.map { it[1] }
    }

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    CREATE_BBMAP_INDEX ( ch_genome_fnas )
    ch_versions = ch_versions.mix(CREATE_BBMAP_INDEX.out.versions)

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_map_reads = BBMAP_BBDUK.out.reads
        ch_versions  = ch_versions.mix(BBMAP_BBDUK.out.versions)
    } else {
        ch_map_reads = FASTQC_TRIMGALORE.out.reads
    }

    //
    // MODULE: Run BBMap
    //
    BBMAP_ALIGN ( ch_map_reads, CREATE_BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: sort, index and run stats on genomes bam
    //
    BAM_SORT_SAMTOOLS ( BBMAP_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    //
    // MODULE: Concatenate gff files
    //
    CONCATENATE_GFF ( 'genomes.gff.gz', ch_genome_gffs.collect() )

    //ch_featurecounts = SAMTOOLS_SORT.out.bam
    BAM_SORT_SAMTOOLS.out.bam
        .combine(CONCATENATE_GFF.out.file)
        .set { ch_featurecounts }

    //
    // MODULE: Run featureCounts
    // 
    // Should be a subworkflow...
    ch_cds_counts = Channel.empty()
    if ( count_cds ) {
        FEATURECOUNTS_CDS ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

        COLLECT_FEATURECOUNTS_CDS   ( FEATURECOUNTS_CDS.out.counts.collect   { it[1] } )
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS_CDS.out.versions)
        ch_cds_counts = COLLECT_FEATURECOUNTS_CDS.out.counts
    }

    ch_rrna_counts = Channel.empty()
    if ( count_rrna ) {
        FEATURECOUNTS_RRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_RRNA.out.versions)

        COLLECT_FEATURECOUNTS_RRNA  ( FEATURECOUNTS_RRNA.out.counts.collect  { it[1] } )
        ch_rrna_counts = COLLECT_FEATURECOUNTS_RRNA.out.counts
    }

    ch_trna_counts = Channel.empty()
    if ( count_trna ) {
        FEATURECOUNTS_TRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TRNA.out.versions)

        COLLECT_FEATURECOUNTS_TRNA  ( FEATURECOUNTS_TRNA.out.counts.collect  { it[1] } )
        ch_trna_counts = COLLECT_FEATURECOUNTS_TRNA.out.counts
    }

    ch_tmrna_counts = Channel.empty()
    if ( count_tmrna ) {
        FEATURECOUNTS_TMRNA ( ch_featurecounts )
        ch_versions = ch_versions.mix(FEATURECOUNTS_TMRNA.out.versions)

        COLLECT_FEATURECOUNTS_TMRNA ( FEATURECOUNTS_TMRNA.out.counts.collect { it[1] } )
        ch_tmrna_counts = COLLECT_FEATURECOUNTS_TMRNA.out.counts
    }

    ch_fcs = Channel.empty()
    ch_fcs = ch_fcs.mix(
        ch_cds_counts,
        ch_rrna_counts,
        ch_trna_counts,
        ch_tmrna_counts
    ).collect()

    //
    // MODULE: Run collect_gene_info
    //
    COLLECT_GENE_INFO ( ch_genome_gffs.collect() )
    ch_versions = ch_versions.mix(COLLECT_GENE_INFO.out.versions)

    //
    // MODULE: Overall statistics
    //
    COLLECT_STATS(
        FASTQC_TRIMGALORE.out.trim_log.map { meta, fastq -> meta.id }.collect(),
        FASTQC_TRIMGALORE.out.trim_log.map { meta, fastq -> fastq[0] }.collect(),
        BAM_SORT_SAMTOOLS.out.idxstats.collect()  { it[1] },
        ch_fcs
    )
    ch_versions = ch_versions.mix(COLLECT_STATS.out.versions)

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMagmap.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
