This section by Swapnil Boyane, usage notes below:

1. Use 01_filter_phylo.sh to fileter the SNPs for phylogeny and allow variant and invariants sites with missing data.
2. Open an R session and run 02_phylo_100kbp.r (configured for 100kbp windows). This creates the tree_100kbp/ directory containing the phylo100kbp_array.sh submission script.
3. Keep  "create_fasta.r",  "create_fasta_from_vcf.r" in "tree_100kbp" directory.
4. Submit "phylo100kbp_array.sh" to array job for running.
5. Use 03_combine_trees_greeenbul.r to combine all individual RAxML_bipartitions.tre files into a single file. 
6. Run 04_maximum_clade_credibility_tree (via DendroPy) to generate a Maximum Clade Credibility tree 
7. 05_astral_tree : Generate a species tree using ASTRAL. 
