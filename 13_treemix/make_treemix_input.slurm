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

dir=$SLURM_SUBMIT_DIR
cd $dir

# env with plink1 and vcftools, also plink2treemix.py from Treemix in directory
source activate plink1-env

# at least 5 alleles to ensure not fixed in virens, 100% complete
vcftools --vcf $dir/analysis_vcfs/greenbul_combined_noZ.vcf \
         --remove-indels --out $dir/analysis_vcfs/greenbul51_noZ_100_mac5_10k \
         --min-alleles 2 --max-alleles 2 --mac 5 --keep sample_list_filter.txt \
         --max-missing 0.9 --max-meanDP 20 --thin 10000 --recode

conda deactivate

# script from here, by Joana Meier: https://github.com/speciationgenomics/scripts/blob/master/vcf2treemix.sh
sh vcf2treemix.sh analysis_vcfs/greenbul51_noZ_100_mac5_10k.recode.vcf greenbul51.clust
