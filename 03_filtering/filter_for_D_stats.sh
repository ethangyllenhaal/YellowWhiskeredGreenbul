#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter_ABBA_BABA
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-35
#SBATCH --mail-user=arrice@ttu.edu
#SBATCH --mail-type=ALL

source ~/anaconda3/etc/profile.d/conda.sh
conda activate vcftools

input_array=$( head -n${SLURM_ARRAY_TASK_ID} scaffolds.txt | tail -n1 )

workdir=/lustre/scratch/arrice/greenbul3

# filtering for ABBA_BABA
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --remove bad_sample.txt \
--max-missing 0.7 --min-alleles 2 --max-meanDP 20 --max-alleles 2 --mac 2 --max-maf 0.49 \
--remove-indels --recode --recode-INFO-all --out ${workdir}/ABBA_BABA/${input_array}__ABBA_BABA

bgzip ${workdir}/ABBA_BABA/${input_array}__ABBA_BABA.recode.vcf

tabix -p vcf ${workdir}/ABBA_BABA/${input_array}__ABBA_BABA.recode.vcf.gz
