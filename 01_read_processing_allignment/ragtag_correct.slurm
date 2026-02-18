#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=ragtag_greenbul
#SBATCH --output=_ragtag_out_%j
#SBATCH --error=_ragtag_error_%j
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu

source activate ragtag-env

cd $SLURM_SUBMIT_DIR

# correct the reference, breaking up scaffolds that disagree with the original ref
ragtag.py correct sylvia_ref/GCF_009819655.1_bSylAtr1.pri_genomic.fna \
	  eula_fragmented/GCA_041108415.1_ASM4110841v1_genomic.fna \
          -t 8 -u -o ./eula_ragtag_correct_chr0

# scaffold the corrected reference to the same ref
ragtag.py scaffold sylvia_ref/GCF_009819655.1_bSylAtr1.pri_genomic.fna \
	  eula_ragtag_correct_chr0/ragtag.correct.fasta \
          -t 8 -u -C -o ./eula_ragtag_correct_chr0
