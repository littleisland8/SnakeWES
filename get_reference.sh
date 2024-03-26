wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O resources/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd resources

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

cd -
