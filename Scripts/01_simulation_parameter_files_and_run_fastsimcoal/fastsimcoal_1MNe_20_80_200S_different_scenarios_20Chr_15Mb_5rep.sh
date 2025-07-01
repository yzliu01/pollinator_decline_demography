#!/bin/bash
#SBATCH --account eDNA
#SBATCH --cpus-per-task 10
#SBATCH --mem 100g
#SBATCH --array=1-xx%xx
#SBATCH --time=05-10:14:00
#SBATCH --error=fastsimcoal_1MNe_20_80_200S_different_scenarios_20Chr_15Mb_5rep.%A_%a.e
#SBATCH --output=fastsimcoal_1MNe_20_80_200S_different_scenarios_20Chr_15Mb_5rep.%A_%a.o
#SBATCH --job-name=fastsimcoal_1MNe_20_80_200S_different_scenarios_20Chr_15Mb_5rep
#SBATCH --mail-type=all #begin,end,fail,all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com #send email notification

cd /home/yzliu/eDNA/faststorage/yzliu/DK_proj/steps/systematic_fsc_test3/100_replicates

par_folder=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/fastsimcoal/para_test3/100_replicates

par=$(ls $par_folder/ft_sim_1000000Ne*20Chr_15Mb.par | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)

## anders advice to use scratch local space
cd $TEMPDIR
cp $par .
mkdir out
cd out

source /home/yzliu/miniforge3/etc/profile.d/conda.sh
## https://stackoverflow.com/questions/61915607/commandnotfounderror-your-shell-has-not-been-properly-configured-to-use-conda
conda activate fastsimcoal2

fsc27093 -i ../*.par -n 100 -q -x -s 0 -m --foldedSFS -c 6 -k 300000000


## copy files to faststorage
eDAN_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/steps/systematic_fsc_test3/100_replicates
cp -r * $eDAN_DIR

exit 0
