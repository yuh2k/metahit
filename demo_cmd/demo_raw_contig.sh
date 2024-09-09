bash ./bin/metahit-modules/raw_contig.sh --bam ./output/alignment/MAP_SORTED.bam --fasta ./output/assembly/final_assembly.fasta --out ./output/downstream

python ./tools/npz_view_tool.py --file ./output/downstream/hic_contact_matrix.npz
