#  metagenomic shotgun sequencing 
./bin/metahit-modules/read_qc.sh -1 test_data/sim_sg_R1.fastq -2 test_data/sim_sg_R2.fastq -o output/readqc/sg -t 4

#  Hi-C sequencing 
./bin/metahit-modules/read_qc.sh -1 test_data/sim_hic_r1.fq -2 test_data/sim_hic_r2.fq -o output/readqc/hic -t 4