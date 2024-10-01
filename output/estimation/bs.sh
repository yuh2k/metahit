#!/bin/bash
echo "Note: This script is designed to run with the amount of memory detected by BBMap."
echo "      If Samtools crashes, please ensure you are running on the same platform as BBMap,"
echo "      or reduce Samtools' memory setting (the -m flag)."
echo "Note: Please ignore any warnings about 'EOF marker is absent'; this is a bug in samtools that occurs when using piped input."
samtools view -bShu output/estimation/SG_map.sam | samtools sort -m 1263M -@ 3 - -o output/estimation/SG_map_sorted.bam
samtools index output/estimation/SG_map_sorted.bam
