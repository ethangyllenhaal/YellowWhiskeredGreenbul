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
echo"
source activate gatk-env

# Make a file with a list of paths for GatherVcfs to use
> $dir/analysis_vcfs/gather_list
while read interval; do
      echo $dir/analysis_vcfs/intervals/greenbul49_90_10k_${interval}.recode.vcf >> \
      $dir/analysis_vcfs/gather_list
done < $dir/scaffolds_noZ.txt

# Run GatherVcfs to combine intervals
gatk GatherVcfs \
     -I $dir/analysis_vcfs/gather_list \
     -O $dir/analysis_vcfs/greenbul49_combined_noZ_90_10k.vcf

conda deactivate
"
source activate plink-env

# Set up a re-named VCF to make it compatible with PLINK/bed2diffs
# Just replaced chromosome names with 1-34
cp $dir/analysis_vcfs/greenbul49_combined_noZ_90_10k.vcf $dir/analysis_vcfs/greenbul49_combined_rename_noZ_90_10k.vcf
sed -i "s/_RagTag//g" $dir/analysis_vcfs/greenbul49_combined_rename_noZ_90_10k.vcf
while IFS=',' read -ra rename; do
    old="${rename[0]}"
    new="${rename[1]}"
    sed -i "s/$old/$new/g" $dir/analysis_vcfs/greenbul49_combined_rename_noZ_90_10k.vcf
done < $dir/chrom_names_noZ

# Make input with PLINK2
plink2 --chr-set 33 --vcf $dir/analysis_vcfs/greenbul49_combined_rename_noZ_90_10k.vcf \
       --make-bed --out $dir/eems/greenbul49_noZ_90_10k


# switch conda envs and run bed2diffs
conda deactivate
source activate eems-env

/home/egyllenh/eems/bed2diffs/src/bed2diffs_v1 --bfile $dir/eems/greenbul49_noZ_90_10k --nthreads 4
