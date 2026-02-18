#!/bin/bash

#SBATCH --chdir=./
#SBATCH --job-name=sumstats
#SBATCH --nodes=1 --ntasks=1
#SBATCH --output=_stats_out
#SBATCH --error=_stats_error
#SBATCH --partition nocona
#SBATCH --time=4:00:00
#SBATCH --mail-user=egyllenh@ttu.edu

# now for veryfine

 # add headers for stats files
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_heterozygosity.txt
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_fst.txt
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_Dxy.txt
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_Tajima.txt
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_theta.txt
grep 'pop1' sumstats_100k_veryfine/NC_089160.1_RagTag__100001__200000__stats.txt > ../combined_sumstats/window_noZ100k_veryfine_pi.txt

# pull out only the relevant stat for each file
for i in $( ls sumstats_100k_veryfine/*stats.txt ); do grep 'Fst' $i >> ../combined_sumstats/window_noZ100k_veryfine_fst.txt; done
for i in $( ls sumstats_100k_veryfine/*stats.txt ); do grep 'Tajima_D' $i >> ../combined_sumstats/window_noZ100k_veryfine_Tajima.txt; done
for i in $( ls sumstats_100k_veryfine/*stats.txt ); do grep 'Dxy' $i >> ../combined_sumstats/window_noZ100k_veryfine_Dxy.txt; done
for i in $( ls sumstats_100k_veryfine/*stats.txt ); do grep 'theta' $i >> ../combined_sumstats/window_noZ100k_veryfine_theta.txt; done
for i in $( ls sumstats_100k_veryfine/*stats.txt ); do grep -w 'pi' $i >> ../combined_sumstats/window_noZ100k_veryfine_pi.txt; done
