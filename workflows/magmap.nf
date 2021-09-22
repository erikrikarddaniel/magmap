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
} else {
    exit 1, "You need to specifiy either a combination of genome gff directory and genome gff extension (\"--genome_gff_dir\" and \"--genome_gff_ext\") or a comma separated list of genome gff files (\"--genome_gffs\")."
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

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { COLLECTDATA } from '../modules/local/collectdata' addParams( options: [:] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { CREATE_BBMAP_INDEX } from '../subworkflows/local/create_bbmap_index' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options       = modules['multiqc']
multiqc_options.args     += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

bbmap_align_options       = modules['bbmap_align']
samtools_sort_options     = modules['samtools_sort']
concatenate_gff_options   = modules['concatenate_gff']
subread_featurecounts_options = modules['subread_featurecounts_cds']

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC        } from '../modules/nf-core/modules/fastqc/main'               addParams( options: modules['fastqc'] )
include { MULTIQC       } from '../modules/nf-core/modules/multiqc/main'              addParams( options: multiqc_options   )
include { BBMAP_ALIGN   } from '../modules/nf-core/modules/bbmap/align/main' addParams( options: bbmap_align_options )
include { SAMTOOLS_SORT } from '../modules/nf-core/modules/samtools/sort/main'        addParams( options: samtools_sort_options )
include { CONCATENATE as CONCATENATE_GFF } from '../modules/local/concatenate'        addParams( options: concatenate_gff_options )
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options )

//
// SUBWORKFLOW: Adapted from rnaseq!
//
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MAGMAP {

    ch_software_versions = Channel.empty()

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
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Concatenate the genome fasta files and create a BBMap index
    //
    CREATE_BBMAP_INDEX ( ch_genome_fnas )
    ch_software_versions = ch_software_versions.mix(CREATE_BBMAP_INDEX.out.bbmap_version.ifEmpty(null))

    //
    // MODULE: Run BBMap
    //
    BBMAP_ALIGN ( FASTQC_TRIMGALORE.out.reads, CREATE_BBMAP_INDEX.out.index )
    ch_software_versions = ch_software_versions.mix(BBMAP_ALIGN.out.version.ifEmpty(null))

    //
    // MODULE: Sort BBMap output
    //
    SAMTOOLS_SORT ( BBMAP_ALIGN.out.bam )
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_SORT.out.version.ifEmpty(null))

    //
    // MODULE: Concatenate gff files
    //
    CONCATENATE_GFF ( ch_genome_gffs.collect() )

    //ch_featurecounts = SAMTOOLS_SORT.out.bam
    SAMTOOLS_SORT.out.bam
        .combine(CONCATENATE_GFF.out.file)
        .set { ch_featurecounts }

    //
    // MODULE: Run featureCounts
    //
    FEATURECOUNTS_CDS ( ch_featurecounts )
    ch_software_versions = ch_software_versions.mix(FEATURECOUNTS_CDS.out.version.ifEmpty(null))

    //
    // MODULE: Run collectdata
    //
    COLLECTDATA ( FEATURECOUNTS_CDS.out.counts.collect { it[1] } )
    ch_software_versions = ch_software_versions.mix(COLLECTDATA.out.r_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(COLLECTDATA.out.dplyr_version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
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
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
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
