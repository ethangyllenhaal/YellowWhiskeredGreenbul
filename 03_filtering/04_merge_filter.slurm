#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=filter_greenbul
#SBATCH --output=_array_filter_out
#SBATCH --error=_array_filter_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --array=1-34

source activate bcftools-env
dir=$SLURM_SUBMIT_DIR

# picks a chromosome based on job array ID
region_array=$( head -n${SLURM_ARRAY_TASK_ID} $dir/scaffolds.txt | tail -n1 )
# merge across individuals for chromosome for filtered VCF
bcftools merge --threads 4 -m id --regions ${region_array} \
	 $dir/03_vcf/*.vcf.gz > \
	 $dir/combined_vcfs/intervals/greenbul_combined_${region_array}.vcf

# lightly filter for summary stats, here for males onle
vcftools --vcf $dir/combined_vcfs/intervals/greenbul_combined_${region_array} --recode --recode-INFO-all \
        --out $dir/analysis_vcfs/intervals/greenbulM24_70_${region_array} --keep keep_male \
        --max-missing 0.7 --remove-indels --max-alleles 2 --max-maf 0.49 --max-meanDP 20

# bgzip and index
bgzip $dir/analysis_vcfs/intervals/greenbulM24_70_${region_array}.recode.vcf
tabix $dir/analysis_vcfs/intervals/greenbulM24_70_${region_array}.recode.vcf.gz

# all individuals
vcftools --vcf $dir/combined_vcfs/intervals/greenbul_combined_${region_array}.vcf --recode --recode-INFO-all \
        --out $dir/analysis_vcfs/intervals/greenbul49_70_${region_array} --keep greenbul49.txt \
        --max-missing 0.7 --remove-indels --max-alleles 2 --max-maf 0.49 --max-meanDP 20

# bgzip and index
bgzip $dir/analysis_vcfs/intervals/greenbul49_70_${region_array}.recode.vcf
tabix $dir/analysis_vcfs/intervals/greenbul49_70_${region_array}.recode.vcf.gz

# filter for general purpose, 90% complete with 10kbp thinning
vcftools --vcf $dir/combined_vcfs/intervals/greenbul_combined_${region_array}.vcf --recode --recode-INFO-all \
	 --out $dir/analysis_vcfs/intervals/greenbul49_90_10k_${region_array} --keep greenbul49.txt \
	 --max-missing 0.9 --mac 2 --max-alleles 2 --max-maf 0.49 --max-meanDP 20 --thin 10000

