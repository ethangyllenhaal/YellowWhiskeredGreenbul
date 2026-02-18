#!/bin/bash

# takes output as command line argument
# run like ./combined_stats.sh path/to/out

echo -e "Name\tMean Depth\tSites genotyped (Millions)\tCleaned reads (Millions)" > $1

# for each sample, output name, mean depth, millions of sites genotypes, and millions of cleaned R1s in that order
for i in $(cat sample_list_full.txt);
do
    paste <(echo "${i}") \
      <(awk 'NR==5' stats/${i}.stats | cut -f4 -d ' ') \
      <(echo $(($(awk 'NR==9' stats/${i}.stats) / 1000000))) \
      <(echo $(($(awk 'NR==11' stats/${i}.stats) / 1000000))) >> $1
done
