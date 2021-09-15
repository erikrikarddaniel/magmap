include { CONCATENATE } from '../../modules/local/concatenate' addParams( )
include { BBMAP_INDEX } from '../../modules/nf-core/modules/bbmap/index/main' addParams( )

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        CONCATENATE(ch_genome_fnas.collect()) 
        BBMAP_INDEX(CONCATENATE.out.file)

    emit: 
    index = BBMAP_INDEX.out.index
}
