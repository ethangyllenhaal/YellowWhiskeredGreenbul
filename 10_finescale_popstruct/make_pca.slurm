#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=filter_greenbul
#SBATCH --output=_filter_out_%j
#SBATCH --error=_filter_error_%j
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu

# submitted from /lustre/scratch/egyllenh/greenbul/
dir=$SLURM_SUBMIT_DIR
cd $dir

# conda env with vcftools
source activate bcftools-env

vcftools --vcf $dir/analysis_vcfs/greenbul_combined_noZ.vcf \
         --remove-indels --out $dir/analysis_vcfs/east22_noZ_90_10k \
         --min-alleles 2 --max-alleles 2 --mac 3 --keep east22.txt \
         --max-missing .90 --max-meanDP 20 --thin 10000 --recode

vcftools --vcf $dir/analysis_vcfs/greenbul_combined_noZ.vcf \
         --remove-indels --out $dir/analysis_vcfs/central22_noZ_90_10k \
         --min-alleles 2 --max-alleles 2 --mac 3 --keep central22.txt \
         --max-missing .90 --max-meanDP 20 --thin 10000 --recode

vcftools --vcf $dir/analysis_vcfs/greenbul_combined_noZ.vcf \
         --remove-indels --out $dir/analysis_vcfs/west5_noZ_90_10k \
         --min-alleles 2 --max-alleles 2 --mac 3 --keep west5.txt \
         --max-missing .90 --max-meanDP 20 --thin 10000 --recode
