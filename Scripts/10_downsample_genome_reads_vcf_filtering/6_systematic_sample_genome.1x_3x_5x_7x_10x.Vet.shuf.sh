#!/bin/sh
#SBATCH --account eDNA
##SBATCH --cpus-per-task 20
#SBATCH --mem 50g
#SBATCH --array=1-10%10
#SBATCH --time=1-16:00:00
#SBATCH --error=Vet_shuf_get_filtered_2025_DP_1x_3x_5x_7x_10x.%A_%a.e
#SBATCH --output=Vet_shuf_get_filtered_2025_DP_1x_3x_5x_7x_10x.%A_%a.o
#SBATCH --job-name=Vet_shuf_get_filtered_2025_DP_1x_3x_5x_7x_10x
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com

## test
##SLURM_ARRAY_TASK_ID=1

## activate (env) tools of variant_calling_mapping
source /home/yzliu/miniforge3/etc/profile.d/conda.sh
conda activate variant_calling_mapping

## run in terminal
#export OPENBLAS_NUM_THREADS=1

## pooled bees
REF_MASKED_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome/ref_masked_bed
## complement (softmasked_regions + gene_regions) bed file - preferred
## test low sequencing depth for Andrena species
New_REF_AndHae_mask_region=$REF_MASKED_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomPas_mask_region=$REF_MASKED_DIR/Bombus_pascuorum-GCA_905332965.1-softmasked_ref_gene.conca_sorted.bed
New_REF_AndMar_mask_region=$REF_MASKED_DIR/Andrena_marginata_GCA_963932335.1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomVet_mask_region=$REF_MASKED_DIR/Bombus_veteranus.hifi_asm_pl2-softmasked_ref_gene.conca_sorted.bed

## ref
REF_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome
REF_AndHae=$REF_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa
REF_BomPas=$REF_DIR/Bombus_pascuorum-GCA_905332965.1-softmasked.fa
REF_AndMar=$REF_DIR/Andrena_marginata_GCA_963932335.1-softmasked.fa
REF_BomVet=$REF_DIR/Bombus_veteranus.hifi_asm_pl2.fa

## vcf
concated_vcf_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/vcf/concated_vcf_each_species_REF
cd $concated_vcf_dir

## earlier (GQ_issue? and attention to VCF file names)
Andhae_New_REF_AndHae_VCF=concated.AndHae_New_REF_AndHae.100kb_g1500x_regions.GQ_issue_solved.vcf.gz
Andhae_New_REF_AndHae_VCF_filter=${Andhae_New_REF_AndHae_VCF/.vcf.gz/}

Andmar_New_REF_AndMar_VCF=concated.AndMar_New_REF_AndMar.100kb_g1500x_regions.all_chr.sorted.GQ_issue_solved.vcf.gz
Andmar_New_REF_AndMar_VCF_filter=${Andmar_New_REF_AndMar_VCF/.vcf.gz/}

Bompas_New_REF_BomPas_VCF=concated.BomPas_New_REF_BomPas.100kb_g1500x_regions.all_chr.sorted.GQ_issue_solved.vcf.gz
Bompas_New_REF_BomPas_VCF_filter=${Bompas_New_REF_BomPas_VCF/.vcf.gz/}

Bomvet_New_REF_BomVet_VCF=concated.BomVet_New_REF_BomVet.100kb_g1500x_regions.cp.vcf.gz
Bomvet_New_REF_BomVet_VCF_filter=${Bomvet_New_REF_BomVet_VCF/.vcf.gz/}

## index vcf files
vcf_list=(
    $Andhae_New_REF_AndHae_VCF
    $Andmar_New_REF_AndMar_VCF
    $Bompas_New_REF_BomPas_VCF
    $Bomvet_New_REF_BomVet_VCF)

index_vcf(){
for vcf in ${vcf_list[@]}
do
echo $vcf
#done
ll -lh $vcf
#vcf_copy=${vcf/vcf.gz/cp.vcf.gz}
#cp $vcf $vcf_copy
echo "indexing: $vcf"
bcftools index $vcf
done
}
#index_vcf

#bcftools index concated.BomVet_New_REF_BomVet.100kb_g1500x_regions.cp.vcf.gz

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>  NEW  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## keep biallelic snp, remove duplicates and normalize snp with long base (bcftools norm -d none/-m-snps), also remove monomorphic snps

# 2023-2024-2025
## Proportion of genome
## 4. Andmar_New_REF_AndMar_VCF
## 40 ind 40*2*1
## 1x 3x 5x 7x 10x

BED_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome/random_prop_sample_genome
BED_LIST=(
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_01.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_02.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_03.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_04.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_05.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_06.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_07.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_08.sort.bed"
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_100b.shuf_subset_09.sort.bed"
## 100% whole genome
"Bombus_veteranus.hifi_asm_pl2.fa.fai.win_whole.subset_10.bed"

)

## remove prefix text using #
#chr_n=${ref_chr_md[i]#ref_genome_md_chr_}
## remove suffix text using %
#chr_n=${chr_n%.fasta}

## bed file input
BED=$(echo ${BED_LIST[@]} | tr " " "\n"| sed -n ${SLURM_ARRAY_TASK_ID}p)
## remove prefix text using #
#OUT_PORT=${BED/#Andrena_marginata_GCA_963932335.1-softmasked.fa.fai.win_whole.subset_}
#echo $OUT_PORT
## remove suffix text using %
#OUT_PORT=${OUT_PORT/%.bed}
#echo $OUT_PORT

## portion value
PROP_LIST=(01 02 03 04 05 06 07 08 09 10)
PROP=$(echo ${PROP_LIST[@]} | tr " " "\n"| sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "P_$PROP"

depth=(58 174 290 416 580)
depth_time=(1x 3x 5x 7x 10x)

for i in ${!depth[@]}
#for depth in 20
do
echo -e "${depth[i]}\t${depth_time[i]}"
#done
bcftools view -R $BED_DIR/$BED $Bomvet_New_REF_BomVet_VCF | \
bcftools filter --soft-filter mask --mask-file $New_REF_BomVet_mask_region | \
bcftools filter --SnpGap 5:indel | \
bcftools norm -a -f $REF_BomVet | \
bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
bcftools filter -e 'AC==0 || AC == AN' | \
bcftools view -e "MEAN(FMT/DP) < "${depth[i]}" || MEAN(FMT/DP) > 1500" \
-Oz -o ./"$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP_"${depth_time[i]}"_1500x_noMS.shuf.P_"$PROP".vcf.gz
done

exit 0

