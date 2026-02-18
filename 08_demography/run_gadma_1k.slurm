#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=gadma_greenbul
#SBATCH --output=_gadma_out_%j
#SBATCH --error=_gadma_error_%j
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu

# submitted from /lustre/scratch/egyllenh/greenbul_tinker/gadma/
dir=$SLURM_SUBMIT_DIR
cd $dir

source activate gadma-env

# initial run
# gadma -p params_file_1k -e extra_params_file

# resumed run
gadma --resume noZ_10pop_thin1k
