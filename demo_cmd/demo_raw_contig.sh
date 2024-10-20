bash ./bin/metahit-modules/raw_contig.sh --bam ./output/alignment/MAP_sorted.bam --fasta ./output/assembly/final_assembly.fasta --out ./output/downstream

python ./npz_view_tool.py --file ./output/normalization/contact_matrix.npz
python ./npz_view_tool.py --file ./output/normalization/normcc_contact_matrix.npz

