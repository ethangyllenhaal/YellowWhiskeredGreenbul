#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter_ABBA_BABA
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-34

source ~/anaconda3/etc/profile.d/conda.sh
conda activate vcftools

# This step should be done after vcfs have been made for each chromosome
# Requires a scaffolds.txt file

# define main working directory
workdir=/whatever/you/put/the/scaffolds_file/in

# define input files from helper file
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filtering for ABBA_BABA
vcftools --vcf ${workdir}/${region_array}.vcf NOT DONE YET! 




--remove excluded_samples.txt --max-missing 0.6 --minGQ 20 --minDP 6 --max-meanDP 40 --min-alleles 2 --max-alleles 2 --mac 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/11_ABBA_BABA/${region_array}__ABBA_BABA

Filter vcfs with an array job that uses vcftools (34 arrays total, 1 per chromosome) and the scaffolds.txt file to delimit chromosomes. 
Filter to only include biallelic sites. Max mean depth should be 20. The W chromosome (NC_089173.1) should already be removed at this point. 
Also, make sure these new vcfs are bgzipped and indexed via tabix the end (Script not included). 
