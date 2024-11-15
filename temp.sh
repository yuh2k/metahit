./bin/metahit-modules/bin_refinement.sh --fasta output/assembly/final_assembly.fasta \
                    --contig_file output/normalization/raw/contig_info.csv \
                    --contact_matrix_file output/normalization/raw/contact_matrix_user.npz \
                    --output output_dir \
                    --marker_file output_dir/contigs.hmmout \
                    --threads 20