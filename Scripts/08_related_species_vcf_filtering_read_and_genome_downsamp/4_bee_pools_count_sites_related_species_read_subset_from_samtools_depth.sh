#!/bin/bash
##SBATCH --account EcoGenetics
##SBATCH --partition normal
##SBATCH --cpus-per-task 8
##SBATCH --mem-per-cpu 128G
#SBATCH --account eDNA
#SBATCH --mem 50G
#SBATCH --time=06:00:00
#SBATCH --array=1%1
#SBATCH --error=4_bee_pools_fun_count_sites_related_species_read_subset_samtools_depth.%A_%a.e
#SBATCH --output=4_bee_pools_fun_count_sites_related_species_read_subset_samtools_depth.%A_%a.o
#SBATCH --job-name=bee_fun_count_sites_related_species_read_subset_samtools_depth
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com

##SLURM_ARRAY_TASK_ID=2

DP_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/bam/bam_stats/samtools_depth

cd $DP_dir

## originally here
#/home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/bee_proj_script/data_process/data_2025/steps/base_overlap_clip_cov_prop_assessment/scripts/
#1_count_sites_real_data_samtools_depth.10_16_species.sh

## files are in sorted order (sort -V)
Hae_related_DP_file_list=(
Andhae.New_REF_AndFul.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHat.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_BomPas.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
)

Mar_related_DP_file_list=(
Andmar.New_REF_AndBic.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndHae.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndHat.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_BomPas.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
)

Pas_related_DP_file_list=(
Bompas.New_REF_ApisMel.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bompas.New_REF_BomCon.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bompas.New_REF_BomHor-.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bompas.New_REF_BomPas.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bompas.New_alt_REF_BomMus.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
)

Vet_related_DP_file_list=(
Bomvet.New_REF_ApisMel.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bomvet.New_REF_BomCon.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bomvet.New_REF_BomHor-.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bomvet.New_REF_BomPas.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
Bomvet.New_REF_BomVet.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
)

## subset reads
Hae_Read_Subset_DP_file_list=(
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_1.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_2.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_3.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_4.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_5.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_6.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_7.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_8.clip_interest_sites.bam.all_sites.depth
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_9.clip_interest_sites.bam.all_sites.depth
## 100%
Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth
)

Mar_Read_Subset_DP_file_list=(
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_1.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_2.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_3.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_4.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_5.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_6.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_7.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_8.clip_interest_sites.bam.all_sites.depth
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.0_9.clip_interest_sites.bam.all_sites.depth
## 100%
Andmar.New_REF_AndMar.filtered.sort.rm_marked_dups.clip_interest_sites.bam.all_sites.depth

)

#1       179     1 1.37725
#1       180     1 1.37725
#1       181     1 1.37725
#1       182     1 1.37725
#1       183     1 1.37725
#1       184     1 1.37725

#awk '{if($4 > 100) i++} END {print i}' $file



#rm *samtool_interest_depth.count

## for related species (conditional depth)
## 78 chr Hae (234 390 546)
get_Hae_related_count_DP_1_3x_5x_7x_10 () {
if [[ $Hae_related_DP_file == *"Andhae.New_REF_AndHae"* ]];then
for mean_DP_filter in 77 233 389 545 779
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Hae_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Hae_related_DP_file >> $Hae_related_species.samtool_interest_depth.all_sites.1_3_5_7_10x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

else

for mean_DP_filter in 233 389 545
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Hae_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Hae_related_DP_file >> $Hae_related_species.samtool_interest_depth.all_sites.3_5_7x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

fi

}

## 80 chr Mar (240 400 560)
get_Mar_related_count_DP_1_3x_5x_7x_10 () {
if [[ $Mar_related_DP_file == *"Andmar.New_REF_AndMar"* ]];then
for mean_DP_filter in 79 239 399 559 799
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Mar_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" '{if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Mar_related_DP_file >> $Mar_related_species.samtool_interest_depth.all_sites.1_3_5_7_10x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

else

for mean_DP_filter in 239 399 559
#79 239 399 559 799
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Mar_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" '{if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Mar_related_DP_file >> $Mar_related_species.samtool_interest_depth.all_sites.3_5_7x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

fi

}

## 68 chr Pas (204 340 476)
get_Pas_related_count_DP_1_3x_5x_7x_10 () {
if [[ $Pas_related_DP_file == *"Bompas.New_REF_BomPas"* ]];then
for mean_DP_filter in 67 203 339 475 679
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Pas_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Pas_related_DP_file >> $Pas_related_species.samtool_interest_depth.all_sites.1_3_5_7_10x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

else

for mean_DP_filter in 203 339 475
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Pas_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Pas_related_DP_file >> $Pas_related_species.samtool_interest_depth.all_sites.3_5_7x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

fi

}

## 58 chr Vet (174 290 416)
get_Vet_related_count_DP_1_3x_5x_7x_10 () {
if [[ $Vet_related_DP_file == *"Bomvet.New_REF_BomVet"* ]];then
for mean_DP_filter in 57 173 289 415 579
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Vet_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Vet_related_DP_file >> $Vet_related_species.samtool_interest_depth.all_sites.1_3_5_7_10x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

else

for mean_DP_filter in 173 289 415
    do
    ## count callable site > xx_depth
    ## with up limit
    echo $Vet_related_DP_file $mean_DP_filter
    awk -v DP="$mean_DP_filter" 'BEGIN { i = 0 } {if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Vet_related_DP_file >> $Vet_related_species.samtool_interest_depth.all_sites.3_5_7x.count
    ## no up limit
    #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count

done

fi

}

#============================================ read subset site count =========================================================

## read subset for mean_DP_filter
## should use below (1x, 3x, 5x, 7x, 10x)
## 78 chr Hae (78 234 390 546 780)
get_read_subset_Hae_count_DP_1_3_5_7_10x () {

#Hae_Read_Subset_DP_file=Andhae.New_REF_AndHae.filtered.sort.rm_marked_dups.0_1.clip_interest_sites.bam.all_sites.depth

for mean_DP_filter in 77 233 389 545 779
    do
    for file_id in $(seq 1 10)
    #for file_id in $(seq 1 9)
        do
        echo $file_id
        Hae_Read_Subset_DP_file=$(echo ${Hae_Read_Subset_DP_file_list[@]} | tr ' ' '\n' | sed -n ${file_id}p)
        ## count callable site > xx_depth
        ## with up limit
        echo $Hae_Read_Subset_DP_file $mean_DP_filter
        awk -v DP="$mean_DP_filter" '{if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Hae_Read_Subset_DP_file >> Andhae.samtool_interest_depth.read_P_01_02_03_04_05_06_07_08_09_10.DP_1_3_5_7_10x.10prop.count
        ## no up limit
        #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count
    done
done
}

## 80 chr Mar (80 240 400 560 800)
get_read_subset_Mar_count_DP_1_3_5_7_10x () {
#for mean_DP_filter in 79 239 399 559 799
for mean_DP_filter in 799
    do
    for file_id in $(seq 1 10)
    #for file_id in $(seq 1 9)
        do
        echo $file_id
        Mar_Read_Subset_DP_file=$(echo ${Mar_Read_Subset_DP_file_list[@]} | tr ' ' '\n' | sed -n ${file_id}p)
        ## count callable site > xx_depth
        ## with up limit
        echo $Mar_Read_Subset_DP_file $mean_DP_filter
        awk -v DP="$mean_DP_filter" '{if($3 > DP && $3 < 1500) i++} END {print DP, i}' $Mar_Read_Subset_DP_file >> Andmar.samtool_interest_depth.read_P_01_02_03_04_05_06_07_08_09_10.DP_1_3_5_7_10x.10prop.count
        ## no up limit
        #awk -v DP="$mean_DP_filter" '{if($3 > DP) i++} END {print DP, i}' $DP_file >> $species.samtool_interest_depth.no_up.count
    done
done
}

get_read_subset_Hae_count_DP_1_3_5_7_10x
get_read_subset_Mar_count_DP_1_3_5_7_10x

exit 0

## in this order
AndMar_REF_AndMar.01P_reads.DP_1x.blueprint
AndMar_REF_AndMar.02P_reads.DP_1x.blueprint
AndMar_REF_AndMar.03P_reads.DP_1x.blueprint
AndMar_REF_AndMar.04P_reads.DP_1x.blueprint
AndMar_REF_AndMar.05P_reads.DP_1x.blueprint
AndMar_REF_AndMar.06P_reads.DP_1x.blueprint
AndMar_REF_AndMar.07P_reads.DP_1x.blueprint
AndMar_REF_AndMar.08P_reads.DP_1x.blueprint
AndMar_REF_AndMar.09P_reads.DP_1x.blueprint
AndMar_REF_AndMar.01P_reads.DP_3x.blueprint


exit 0

cat Andmar.samtool_interest_depth.read_P_01_02_03_04_05_06_07_08_09_10.DP_1_3_5_7_10x.latest.count
79 162396
79 103071465
79 137214937
79 139730053
79 140621304
79 141057748
79 141313821
79 141479641
79 141593719

## add comma
ls *1_3_5_7_10x.latest.count
Andhae=Andhae.samtool_interest_depth.read_P_01_02_03_04_05_06_07_08_09_10.DP_1_3_5_7_10x.latest.count
Andmar=Andmar.samtool_interest_depth.read_P_01_02_03_04_05_06_07_08_09_10.DP_1_3_5_7_10x.latest.count

out_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/sofwtare/stairway_plot_v2/stairway_plot_v2.1.2/stairway_plot_blueprint/pool_downsample_reads/templates

awk '{print $2}' $Andhae > $out_dir/Hae_39_number_called_sites_downsample_reads_across_genome.txt
awk '{printf "%'\''d\n", $2}' $Andhae > $out_dir/Hae_39_number_called_sites_downsample_reads_across_genome.new.txt
awk '{print $2}' $Andmar > $out_dir/Mar_40_number_called_sites_downsample_reads_across_genome.txt
awk '{printf "%'\''d\n", $2}' $Andmar > $out_dir/Mar_40_number_called_sites_downsample_reads_across_genome.new.txt

exit 0

============================================================================================================
## attendtion to count of elements
## get total count of elements

#n=$(echo ${#Hae_related_DP_file_list[@]}) # 4 elements
#n=$(echo ${#Mar_related_DP_file_list[@]}) # 5 elements
#n=$(echo ${#Pas_related_DP_file_list[@]}) # 5 elements
#n=$(echo ${#Vet_related_DP_file_list[@]}) # 5 elements

#for SLURM_ARRAY_TASK_ID in $(seq 1 $n)
#do
#echo $SLURM_ARRAY_TASK_ID
#done

    ## related species
    Hae_related_DP_file=$(echo ${ref[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
    Hae_related_species=${Hae_related_DP_file%%.filtered*}

    Mar_related_DP_file=$(echo ${ref[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
    Mar_related_species=${Mar_related_DP_file%%.filtered*}

    Pas_related_DP_file=$(echo ${ref[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
    Pas_related_species=${Pas_related_DP_file%%.filtered*}

    Vet_related_DP_file=$(echo ${ref[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
    Vet_related_species=${Vet_related_DP_file%%.filtered*}


functions=(
#get_Hae_related_count_DP_1_3x_5x_7x_10
get_Mar_related_count_DP_1_3x_5x_7x_10
#get_Pas_related_count_DP_1_3x_5x_7x_10
#get_Vet_related_count_DP_1_3x_5x_7x_10
)

arrays=(
#Hae_related_DP_file_list
Mar_related_DP_file_list
#Pas_related_DP_file_list
#Vet_related_DP_file_list
)

species_variables=(
#Hae_related
Mar_related
#Pas_related
#Vet_related
)

for index in ${!functions[@]}
do
    fun=${functions[index]}
    array=${arrays[index]}
    species_prefix=${species_variables[index]}
    declare -n ref="$array"
    ## get total item count of each array
    #eval "echo \${#${array}[@]}"
    n=$(echo ${#ref[@]})
    echo -e "Function: $fun\tAll tasks: $n"
#done

    for ((SLURM_ARRAY_TASK_ID=1;SLURM_ARRAY_TASK_ID<=n;SLURM_ARRAY_TASK_ID++))
    do

    ## related species
    DP_file=$(echo ${ref[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
    #DP_file="${ref[SLURM_ARRAY_TASK_ID-1]}"
    related_species=${DP_file%%.filtered*}
    echo -e "$DP_file\t$related_species"
    ## variables in function
    # $Vet_related_DP_file
    # $Vet_related_species
    
    ## Dynamically assign the correct global variables for the function
    ## Use \"$dp_file\" inside eval to safely assign strings with spaces or special characters.
    eval "${species_prefix}_DP_file=\"$DP_file\""
    eval "${species_prefix}_species=\"$related_species\""

    ## run function
    echo -e "Run function: $fun\tTask:$SLURM_ARRAY_TASK_ID"
    $fun

    done
    ## wait until all tasks of each function finish
    echo -e "Finish all tasks\t$fun"

done

## compile site count for each species and related
/home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/bee_proj_script/data_process/data_2025/downsample_bam_genome_1_3_5_7_10x_norm_a/
6_subset_genome_prop_1_2_3_4_5_6_8_9_10_count_snps_and_spaces_for_sfs_1_3_5_7_10x_thin_snps.new_method.sh

============================================================================================================


exit 0

