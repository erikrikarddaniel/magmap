//
// Create a BBMap index out of a set of fasta nucleotide files.
//

params.bbmap_index_options = [:]

include { CONCATENATE } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: false ] ] )
include { BBMAP_INDEX } from '../../modules/nf-core/modules/bbmap/index/main' addParams( options: params.bbmap_index_options )

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        CONCATENATE('genomes.fna.gz', ch_genome_fnas.collect()) 
        BBMAP_INDEX(CONCATENATE.out.file)

    emit: 
    index         = BBMAP_INDEX.out.index
    versions      = BBMAP_INDEX.out.versions
}
