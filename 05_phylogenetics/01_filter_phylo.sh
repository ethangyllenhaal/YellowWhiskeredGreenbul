#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter_0.7miss
#SBATCH --nodes=1 --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-33

source activate bcftools

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} scaffolds.txt | tail -n1 )


# define main working directory
workdir=/lustre/scratch/sboyane/01_greenbuls


#SNP and invariant site output, 30% max missing data, no indels
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --max-missing 0.7 --keep keeplist.txt  --max-alleles 2 --remove-indels --recode --recode-INFO-all --out ${workdir}/09_phylogeny/${input_array}_RagTag

# bgzip and tabix index files 

~/anaconda3/bin/bgzip ${workdir}/09_phylogeny/${input_array}_RagTag.recode.vcf

#tabix
~/anaconda3/bin/tabix -p vcf ${workdir}/09_phylogeny/${input_array}_RagTag.recode.vcf.gz

