#!/bin/sh
#SBATCH --account eDNA
##SBATCH --cpus-per-task 20
#SBATCH --mem 50g
#SBATCH --array=1%1
##SBATCH --array=1-9%9
#SBATCH --time=06:00:00
#SBATCH --error=Hae_downsample_read_2025_DP_1x_3x_5x_7x_10x.%A_%a.e
#SBATCH --output=Hae_downsample_read_2025_DP_1x_3x_5x_7x_10x.%A_%a.o
#SBATCH --job-name=Hae_downsample_read_2025_DP_1x_3x_5x_7x_10x
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com

## activate (env) tools of variant_calling_mapping
source /home/yzliu/miniforge3/etc/profile.d/conda.sh
conda activate variant_calling_mapping

## pooled bees
REF_MASKED_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome/ref_masked_bed
## complement (softmasked_regions + gene_regions) bed file - preferred

New_REF_AndHae_mask_region=$REF_MASKED_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked_ref_gene.conca_sorted.bed
New_REF_AndMar_mask_region=$REF_MASKED_DIR/Andrena_marginata_GCA_963932335.1-softmasked_ref_gene.conca_sorted.bed

## ref
REF_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome
REF_AndHae=$REF_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa
REF_AndMar=$REF_DIR/Andrena_marginata_GCA_963932335.1-softmasked.fa

## vcf
concated_vcf_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/vcf/concated_vcf_each_species_REF
cd $concated_vcf_dir

## vcf from here
## /home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/bee_proj_script/data_process/data_2025/downsample_bam_genome_1_3_5_7_10x_norm_a/
## 4_5_2024_bam_and_genome_(not_to_do_with_bam)_read_downsample_concate_vcf_filtering_genome_fraction_depth_pgzip_gzip_decompress.sh
## vcf files for each percentage of downsampling

#ls Andhae*_0_*100kb_1500x_region.vcf.gz | sort -V
#Andhae.New_REF_AndHae_0_1.100kb_1500x_region.vcf.gz
#Andhae.New_REF_AndHae_0_2.100kb_1500x_region.vcf.gz

## 100% read group vcf files: from genome 100% genome sampling
## /home/yzliu/eDNA/faststorage/yzliu/DK_proj/population_genomics/bee_proj_script/data_process/data_2025_downsample/Scripts/downsample_bam/
## 6_systematic_sample_genome.1x_3x_5x_7x_10x.Hae.template.sh

MAX_BIN=(
"39"
#"40"
#"34"
#"29"
)
species=(
Hae
#Mar
#Pas
#Vet
)

vcf_list=(
Andhae.New_REF_AndHae_0_1.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_2.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_3.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_4.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_5.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_6.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_7.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_8.100kb_1500x_region.clip.vcf.gz
Andhae.New_REF_AndHae_0_9.100kb_1500x_region.clip.vcf.gz
)

## 100% data at 1,3,5,7,10x
#AndHae_P_10_1x=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.vcf.gz
#AndHae_P_10_3x=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP234_1500x_noMS.vcf.gz
#AndHae_P_10_5x=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP390_1500x_noMS.vcf.gz
#AndHae_P_10_7x=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP546_1500x_noMS.vcf.gz
#AndHae_P_10_10x=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP780_1500x_noMS.vcf.gz

## index vcf files
index_vcf(){
for vcf in ${vcf_list[@]}
do
echo $vcf
#done
ll -lh $vcf
vcf_copy=${vcf/vcf.gz/cp.vcf.gz}
#cp $vcf $vcf_copy
echo "indexing: $vcf"
bcftools index $vcf
done
}
#index_vcf

## keep biallelic snp, remove duplicates and normalize snp with long base (bcftools norm -d none/-m-snps), also remove monomorphic snps
## output sfs directory
output_SFS_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/SFS_data

# Hae
## 39 ind 78*1 1x
#depth=(78 256 390 546 780) # 256 is wrong due to wrong calculation
MAX_BIN=39
depth=(78 234 390 546 780)
depth_time=(1x 3x 5x 7x 10x)
#depth=(390 546 780)
#depth_time=(5x 7x 10x)

#vcf_file=$(ls Andhae.New_REF_AndHae_*.100kb_1500x_region.clip.vcf.gz | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)
vcf_file=$(echo ${vcf_list[@]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
## Andhae.New_REF_AndHae_0_1.100kb_1500x_region.vcf.gz
output_prefix=${vcf_file/.100kb_1500x_region.clip.vcf.gz/}
## Andhae.New_REF_AndHae_0_1.


for i in ${!depth[@]}
#for depth in 20
do
    ## tell what is doing
    echo -e "filtering\t${depth[i]}\t${depth_time[i]}"
#done
    bcftools filter --soft-filter mask --mask-file $New_REF_AndHae_mask_region $vcf_file | \
    bcftools filter --SnpGap 5:indel | \
    bcftools norm -a -f $REF_AndHae | \
    bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < "${depth[i]}" || MEAN(FMT/DP) > 1500" | \
    bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk -v max="$MAX_BIN" '
        {
            count[$2] = $1
            if ($2 > max) max = $2
        }
        END {
            for (i = 1; i <= max; i++)
                printf "%d ", (i in count ? count[i] : 0)
        }
        END { print "" }' > $output_SFS_dir/"$output_prefix"P_downsample_reads."${depth_time[i]}".equal_self.sfs
    
    #awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$output_prefix"P_downsample_reads."${depth_time[i]}".equal_self.sfs
    ## output file name example
    #Andhae.New_REF_AndHae_0_1.P_downsample_reads.1x.equal_self.sfs
    #Andhae.New_REF_AndHae_0_1.P_downsample_reads.3x.equal_self.sfs

    ## print the sfs in terminal
    cat $output_SFS_dir/"$output_prefix"P_downsample_reads."${depth_time[i]}".equal_self.sfs

## save subset vcf files
    ##bcftools view -e "MEAN(FMT/DP) < "${depth[i]}" || MEAN(FMT/DP) > 1500" \
    ##-Oz -o ./"$output_prefix"P_downsample_reads."${depth_time[i]}".vcf.gz
    ## count snp
    ##bcftools view -H $vcf_file | wc -l > $output_SFS_dir/"$output_prefix"P_downsample_reads."${depth_time[i]}".snp.count

done

exit 0

