//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CONCATENATE as FIRST_CONCATENATE  } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: false ], args: "| grep -E '\\t'" ] )
include { CONCATENATE as SECOND_CONCATENATE } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: false ] ] )
include { GENOMEINDEX                       } from '../../modules/local/genomeindex' addParams( [ options: [ publish_files: false ] ] )
include { CONCATENATE as GINDEX_CONCATENATE } from '../../modules/local/concatenate' addParams( [ options: [ publish_files: true  ] ] )

workflow CONCATENATE_GFFS {
    take: ch_genome_gffs

    main:
        FIRST_CONCATENATE  ('', ch_genome_gffs.collate(1000)) 
        SECOND_CONCATENATE ('genomes.gff.gz', FIRST_CONCATENATE.out.file.collect())
        GENOMEINDEX        (ch_genome_gffs.collate(1000)) 
        GINDEX_CONCATENATE ('genomeindex.tsv.gz', GENOMEINDEX.out.genomes2id.collect())

    emit: 
    gff    = SECOND_CONCATENATE.out.file
    gindex = GINDEX_CONCATENATE.out.file
}
