#!/usr/bin/env python
from __future__ import print_function
import sys
# This script summarizes the statistics of each bin by parsing
# the checkm_folder/storage/bin_stats_ext.tsv file of the CheckM output


if len(sys.argv) == 3:
    binner = sys.argv[2]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
elif len(sys.argv) == 4:
    source = {}
    for line in open(sys.argv[3]):
        cut = line.strip().split("\t")
        source[cut[0]] = cut[7]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
else:
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize")


for line in open(sys.argv[1]):
    if "Name" in line:
        continue
    content = line.strip().split("\t")
    # print(content)

    # if dic["Completeness"]<20 or dic["Contamination"]>10: continue

    if len(sys.argv) == 3:
        print("\t".join([content[0], str(content[1])[:5],
                         str(content[2])[:5], str(content[9])[:5],
                         content[3], str(content[6]),
                         str(content[8]), binner]))

    elif len(sys.argv) == 4:
        print("\t".join([content[0], str(content[1])[:5],
                         str(content[2])[:5], str(content[9])[:5],
                         content[3], str(content[6]),
                         str(content[8]), source[content[0]]]))

    else:
        print("\t".join([content[0], str(content[1])[:5],
                         str(content[2])[:5], str(content[9])[:5],
                         content[3], str(content[6]),
                         str(content[8])]))
