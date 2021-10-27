//
// Create a BBMap index out of a set of fasta nucleotide files.
//

params.bbmap_index_options = [:]

include { CONCATENATE as FIRST_CONCATENATE  } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: false ] ] )
include { CONCATENATE as SECOND_CONCATENATE } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: false ] ] )
include { BBMAP_INDEX                       } from '../../modules/nf-core/modules/bbmap/index/main' addParams( options: params.bbmap_index_options )

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        FIRST_CONCATENATE  ('', ch_genome_fnas.collate(1000)) 
        SECOND_CONCATENATE ('genomes.fna.gz', FIRST_CONCATENATE.out.file.collect())
        BBMAP_INDEX        (SECOND_CONCATENATE.out.file)

    emit: 
    index         = BBMAP_INDEX.out.index
    versions      = BBMAP_INDEX.out.versions
}
