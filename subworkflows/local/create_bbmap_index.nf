include { CONCATENATE } from '../../modules/local/concatenate' addParams( )
include { BBMAP_INDEX } from '../../modules/nf-core/modules/bbmap/index/main' addParams( )

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        CONCATENATE(ch_genome_fnas.collect()) //.file.map { ch_genomes }
        BBMAP_INDEX(CONCATENATE.out.file)

    emit: 
    BBMAP_INDEX.out.index
}
