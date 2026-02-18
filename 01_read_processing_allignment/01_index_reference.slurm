#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=index_greenbul
#SBATCH --output=_index_out
#SBATCH --error=_index_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu

# working directory path
dir=/lustre/scratch/egyllenh/redo_greenbul/eula_ragtag_correct_chr0/
# reference name in that path
ref=ragtag.scaffold

cd $dir

# environment with bwa, samtools, and picard
source activate bcftools-env

# BWA index
bwa index ${ref}.fasta

# samtools fasta index (most human readable)
samtools faidx ${ref}.fasta

# picard sequence dictionary
picard CreateSequenceDictionary \
R=${ref}.fasta \
O=${ref}.dict


