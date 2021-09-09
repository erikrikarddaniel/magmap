include { CONCATENATE } from '../../modules/local/concatenate' addParams( )

workflow CREATE_BBMAP_INDEX {
    take: ch_genome_fnas

    main:
        CONCATENATE(ch_genome_fnas.collect()) //.file.map { ch_genomes }

    emit: 
    CONCATENATE.out.file
}
