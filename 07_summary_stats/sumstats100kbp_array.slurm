#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=sumstats
#SBATCH --nodes=1 --ntasks=1
#SBATCH --output=_stat_out
#SBATCH --error=_stat_error
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-498

source activate bcftools-env

# Set the number of runs that each SLURM task should do
PER_TASK=19

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=$END_NUM; run++ )); do
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run

	chrom_array=$( head -n${run} stats_helper_chrom.txt | tail -n1 )

	start_array=$( head -n${run} stats_helper_start.txt | tail -n1 )

	end_array=$( head -n${run} stats_helper_end.txt | tail -n1 )

	gunzip -cd /home/egyllenh/scratch/redo_greenbul/analysis_vcfs/intervals/greenbul49_70_NC_089172.1_RagTag.recode.vcf.gz | grep "#" > /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	tabix /home/egyllenh/scratch/redo_greenbul/analysis_vcfs/intervals/greenbul49_70_${chrom_array}.recode.vcf.gz ${chrom_array}:${start_array}-${end_array} >> /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf > /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf

	Rscript calculate_windows.r /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf greenbul_popmap.txt
    Rscript calculate_windows_fine.r /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf greenbul_popmap_fine.txt
	Rscript calculate_windows_veryfine.r /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${s\
tart_array}__${end_array}.simple.vcf greenbul_popmap_veryfine.txt

	rm /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf
	rm /home/egyllenh/scratch/redo_greenbul/winstats_noZ/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf
done
