#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=filter_greenbul
#SBATCH --output=_filter_out_%j
#SBATCH --error=_filter_error_%j
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu

dir=$SLURM_SUBMIT_DIR
cd $dir

# environment with GATK4
source activate gatk-env

# Make a file with a list of paths for GatherVcfs to use
> $dir/analysis_vcfs/gather_list
while read interval; do
      echo $dir/combined_vcfs/intervals/greenbul_combined_${interval}.vcf >> \
      $dir/analysis_vcfs/gather_list
done < $dir/scaffolds_noZ.txt

# Run GatherVcfs to combine intervals
gatk GatherVcfs \
     -I $dir/analysis_vcfs/gather_list \
     -O $dir/analysis_vcfs/greenbul_combined_noZ.vcf

conda deactivate

# environment with vcftools
source activate bcftools-env

# make 1k thinned vcf for "high quality" good site (mean dp <= 1.5x overall mean)
vcftools --vcf $dir/analysis_vcfs/greenbul_combined_noZ.vcf \
         --remove-indels --out $dir/analysis_vcfs/greenbul47_full1k_noZ_90 \
	 --thin 1000 \
         --max-alleles 2 --keep keep_gadma.txt \
         --max-missing .90 --max-meanDP 15 --recode

# Filter to biallelic sites without further thinning
vcftools --vcf $dir/analysis_vcfs/greenbul47_full1k_noZ_90.recode.vcf \
         --remove-indels --out $dir/analysis_vcfs/greenbul47_biallelic1k_noZ_90 \
         --min-alleles 2 --max-alleles 2 --keep keep_gadma.txt \
         --max-missing .90 --max-meanDP 15 --recode

vcftools --vcf $dir/analysis_vcfs/greenbul47_full1k_noZ_90.recode.vcf \
         --remove-indels --out $dir/analysis_vcfs/greenbul47_invariant1k_noZ_90 \
         --max-alleles 1 --keep keep_gadma.txt \
         --max-missing .90 --max-meanDP 15 --recode

conda deactivate

source activate easySFS

# run easySFS, first in preview mode then for real based on the preview
# note that preview mode has to be run interactivately if popmap and vcf don't match
# projection based on 10k thinned dataset with Z was --proj=34,38,8

# sequence length: 112778 variant, 632399 total

# not thinned
/home/egyllenh/easySFS/easySFS.py \
    -i $dir/analysis_vcfs/greenbul47_biallelic1k_noZ_90.recode.vcf \
    --total-length 632399 \
    -p $dir/gadma/popmap_gadma -o gadma/easysfs_noZ_thin1k -a -f --proj=34,38,8
