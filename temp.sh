python bin/metahit-scripts/bin_refinement.py --method metacc \
--contig_file output/normalization/raw/contig_info.csv \
--hic_matrix output/normalization/raw/contact_matrix_metacc.npz \
--output output/bins \
--fasta output/assembly/final_assembly.fasta \
--num_gene 100 \
--seed 42

./bin/metahit-modules/normalization.sh normcc \
--contig_file /home/hh/Documents/metahit1018/metahit/output/normalization/raw/contig_info.csv \
--contact_matrix_file /home/hh/Documents/metahit1018/metahit/output/normalization/raw/contact_matrix_metacc.npz \
--output /home/hh/Documents/metahit1018/metahit/output/normalization/normcc \
--thres 1


bin/metahit-modules/raw_contact1.sh --bam "output/alignment/sorted_map.bam" \
    --fasta "output/assembly/final_assembly.fasta" \
    --out "./output1" \
    --enzyme "HindIII"

