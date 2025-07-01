#!/bin/bash
#SBATCH --account eDNA
#SBATCH --cpus-per-task 6
#SBATCH --mem 20g
#SBATCH --array=1-xxx%xxx
#SBATCH --time=02:00:00
#SBATCH --error=fastsimcoal_Ne_20_80_200S_different_scenarios_20Chr_15Mb_5rep.%A_%a.e
#SBATCH --output=fastsimcoal_Ne_20_80_200S_different_scenarios_20Chr_15Mb_5rep.%A_%a.o
#SBATCH --job-name=fastsimcoal_Ne_20_80_200S_different_scenarios_20Chr_15Mb_5rep
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com


cd /home/yzliu/eDNA/faststorage/yzliu/DK_proj/steps/systematic_fsc_test3/100000Ne_5rep

par_folder=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/fastsimcoal/para_test3/100000Ne

par=$(ls $par_folder/ft_sim_100000Ne*1E_*05i_20Chr_15Mb.par | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)

## anders advice to use scratch local space
cd $TEMPDIR
cp $par .
mkdir out
cd out

source /home/yzliu/miniforge3/etc/profile.d/conda.sh
## https://stackoverflow.com/questions/61915607/commandnotfounderror-your-shell-has-not-been-properly-configured-to-use-conda
conda activate fastsimcoal2

fsc27093 -i ../*.par -n 5 -q -x -s 0 -m --foldedSFS -c 6 -g -G -k 300000000

## copy files to faststorage
eDAN_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/steps/systematic_fsc_test3/100000Ne_5rep
cp -r * $eDAN_DIR
