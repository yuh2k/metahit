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



bin/metahit-modules/bin_refinement.sh output/assembly/final_assembly.fasta output/alignment/sorted_map.bam output/bins \                                    
-t 10 --enzyme DpnII --metacc-min-len 1000 --metacc-min-signal 2 --bin3c-min-len 1000 \
--bin3c-min-signal 5 --thres 0.05 --cover 
bin/metahit-modules/bin_refinement.sh output/assembly/final_assembly.fasta output/alignment/sorted_map.bam output/bins \
-t 10 --enzyme DpnII --metacc-min-len 1000 --metacc-min-signal 2 --bin3c-min-len 1000 \
--bin3c-min-signal 5 --thres 0.05 --cover


echo "[INFO] Running Bin Refinement Process..."
./metahit.py bin_refinement --fasta "output/assembly/final_assembly.fasta" \
  --bam "output/alignment/sorted_map.bam" \
  --output "output/bins" \
  -t 10 \
  --enzyme DpnII \
  --metacc-min-len 1000 \
  --metacc-min-signal 2 \
  --bin3c-min-len 1000 \
  --bin3c-min-signal 1 \
  --thres 0.01 \
  --cover