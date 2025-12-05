#!/bin/sh
#SBATCH --account eDNA
##SBATCH --cpus-per-task 20
#SBATCH --mem 10g
#SBATCH --array=3,7,12%3
##SBATCH --array=1-19%19
#SBATCH --time=06:00:00
#SBATCH --error=related_species_get_vcf_ind_DP_1_3x_5x_7x_10.%A_%a.e
#SBATCH --output=related_species_get_vcf_ind_DP_1_3x_5x_7x_10.%A_%a.o
#SBATCH --job-name=related_species_get_vcf_ind_DP_1_3x_5x_7x_10
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com

## activate (env) tools of variant_calling_mapping
source /home/yzliu/miniforge3/etc/profile.d/conda.sh
conda activate variant_calling_mapping

## pooled bees
REF_MASKED_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome/ref_masked_bed

## complement (softmasked_regions + gene_regions) bed file
New_REF_AndHae_mask_region=$REF_MASKED_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked_ref_gene.conca_sorted.bed
New_REF_AndHat_mask_region=$REF_MASKED_DIR/Andrena_hattorfiana-GCA_944738655.1-softmasked_ref_gene.conca_sorted.bed
New_REF_AndMar_mask_region=$REF_MASKED_DIR/Andrena_marginata_GCA_963932335.1-softmasked_ref_gene.conca_sorted.bed
New_REF_AndBic_mask_region=$REF_MASKED_DIR/Andrena_bicolor-GCA_960531205.1-softmasked_ref_gene.conca_sorted.bed

New_REF_AndFul_mask_region=$REF_MASKED_DIR/Andrena_fulva-GCA_946251845.1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomCon_mask_region=$REF_MASKED_DIR/Bombus_confusus-GCA_014737475.1_ASM1473747v1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomHor_mask_region=$REF_MASKED_DIR/Bombus_hortorum-GCA_905332935.1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomPas_mask_region=$REF_MASKED_DIR/Bombus_pascuorum-GCA_905332965.1-softmasked_ref_gene.conca_sorted.bed
New_REF_alt_BomMus_mask_region=$REF_MASKED_DIR/Bombus_muscorum-GCA_963971125.1-softmasked_ref_gene.conca_sorted.bed
New_REF_BomVet_mask_region=$REF_MASKED_DIR/Bombus_veteranus.hifi_asm_pl2-rm_contamination_softmasked_ref_gene.conca_sorted.bed

New_REF_ApisMel_mask_region=$REF_MASKED_DIR/Amel_HAv-GCF_003254395.2-softmasked_ref_gene.conca_sorted.bed

## ref
REF_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome

REF_AndHae=$REF_DIR/Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa
REF_AndHat=$REF_DIR/Andrena_hattorfiana-GCA_944738655.1-softmasked.fa
REF_AndFul=$REF_DIR/Andrena_fulva-GCA_946251845.1-softmasked.fa
REF_AndBic=$REF_DIR/Andrena_bicolor-GCA_960531205.1.fa
REF_AndMar=$REF_DIR/Andrena_marginata_GCA_963932335.1-softmasked.fa

REF_BomPas=$REF_DIR/Bombus_pascuorum-GCA_905332965.1-softmasked.fa
REF_BomCon=$REF_DIR/Bombus_confusus-GCA_014737475.1_ASM1473747v1-softmasked.fa
REF_BomHor=$REF_DIR/Bombus_hortorum-GCA_905332935.1-softmasked.fa
REF_BomMus=$REF_DIR/Bombus_muscorum-GCA_963971125.1.fa
REF_BomVet=$REF_DIR/Bombus_veteranus.hifi_asm_pl2.fa

REF_ApisMel=$REF_DIR/Apis_mellifera_HAv-GCF_003254395.2-softmasked.fa


## vcf_dir
concated_vcf_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/vcf/concated_vcf_each_species_REF
cd $concated_vcf_dir

Andhae_New_REF_AndHae_VCF=AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.vcf.gz
Andhae_New_REF_AndHae_VCF_filter=${Andhae_New_REF_AndHae_VCF/.vcf.gz/}
Andhae_New_REF_AndHat_VCF=AndHae_New_REF_AndHat.individual_100kb_1500x_region.clip.vcf.gz
Andhae_New_REF_AndHat_VCF_filter=${Andhae_New_REF_AndHat_VCF/.vcf.gz/}
Andhae_New_REF_AndFul_VCF=AndHae_New_REF_AndFul.individual_100kb_1500x_region.clip.vcf.gz
Andhae_New_REF_AndFul_VCF_filter=${Andhae_New_REF_AndFul_VCF/.vcf.gz/}
Andhae_New_REF_BomPas_VCF=AndHae_New_REF_BomPas.individual_100kb_1500x_region.clip.vcf.gz
Andhae_New_REF_BomPas_VCF_filter=${Andhae_New_REF_BomPas_VCF/.vcf.gz/}

Andmar_New_REF_AndMar_VCF=AndMar_New_REF_AndMar.individual_100kb_1500x_region.clip.vcf.gz
Andmar_New_REF_AndMar_VCF_filter=${Andmar_New_REF_AndMar_VCF/.vcf.gz/}
Andmar_New_REF_AndBic_VCF=AndMar_New_REF_AndBic.individual_100kb_1500x_region.clip.vcf.gz
Andmar_New_REF_AndBic_VCF_filter=${Andmar_New_REF_AndBic_VCF/.vcf.gz/}
Andmar_New_REF_AndHat_VCF=AndMar_New_REF_AndHat.individual_100kb_1500x_region.clip.vcf.gz
Andmar_New_REF_AndHat_VCF_filter=${Andmar_New_REF_AndHat_VCF/.vcf.gz/}
Andmar_New_REF_AndHae_VCF=AndMar_New_REF_AndHae.individual_100kb_1500x_region.clip.vcf.gz
Andmar_New_REF_AndHae_VCF_filter=${Andmar_New_REF_AndHae_VCF/.vcf.gz/}
Andmar_New_REF_BomPas_VCF=AndMar_New_REF_BomPas.individual_100kb_1500x_region.clip.vcf.gz
Andmar_New_REF_BomPas_VCF_filter=${Andmar_New_REF_BomPas_VCF/.vcf.gz/}

Bompas_New_REF_BomPas_VCF=BomPas_New_REF_BomPas.individual_100kb_1500x_region.clip.vcf.gz
Bompas_New_REF_BomPas_VCF_filter=${Bompas_New_REF_BomPas_VCF/.vcf.gz/}
Bompas_New_alt_REF_BomMus_VCF=BomPas_New_alt_REF_BomMus.individual_100kb_1500x_region.clip.vcf.gz
Bompas_New_alt_REF_BomMus_VCF_filter=${Bompas_New_alt_REF_BomMus_VCF/.vcf.gz/}
Bompas_New_REF_BomHor_VCF=BomPas_New_REF_BomHor.individual_100kb_1500x_region.clip.vcf.gz
Bompas_New_REF_BomHor_VCF_filter=${Bompas_New_REF_BomHor_VCF/.vcf.gz/}
Bompas_New_REF_BomCon_VCF=BomPas_New_REF_BomCon.individual_100kb_1500x_region.clip.vcf.gz
Bompas_New_REF_BomCon_VCF_filter=${Bompas_New_REF_BomCon_VCF/.vcf.gz/}
Bompas_New_REF_ApisMel_VCF=BomPas_New_REF_ApisMel.individual_100kb_1500x_region.clip.vcf.gz
Bompas_New_REF_ApisMel_VCF_filter=${Bompas_New_REF_ApisMel_VCF/.vcf.gz/}

Bomvet_New_REF_BomVet_VCF=BomVet_New_REF_BomVet.individual_100kb_1500x_region.clip.vcf.gz
Bomvet_New_REF_BomVet_VCF_filter=${Bomvet_New_REF_BomVet_VCF/.vcf.gz/}
Bomvet_New_REF_BomPas_VCF=BomVet_New_REF_BomPas.individual_100kb_1500x_region.clip.vcf.gz
Bomvet_New_REF_BomPas_VCF_filter=${Bomvet_New_REF_BomPas_VCF/.vcf.gz/}
Bomvet_New_REF_BomHor_VCF=BomVet_New_REF_BomHor.individual_100kb_1500x_region.clip.vcf.gz
Bomvet_New_REF_BomHor_VCF_filter=${Bomvet_New_REF_BomHor_VCF/.vcf.gz/}
Bomvet_New_REF_BomCon_VCF=BomVet_New_REF_BomCon.individual_100kb_1500x_region.clip.vcf.gz
Bomvet_New_REF_BomCon_VCF_filter=${Bomvet_New_REF_BomCon_VCF/.vcf.gz/}
Bomvet_New_REF_ApisMel_VCF=BomVet_New_REF_ApisMel.individual_100kb_1500x_region.clip.vcf.gz
Bomvet_New_REF_ApisMel_VCF_filter=${Bomvet_New_REF_ApisMel_VCF/.vcf.gz/}

## keep biallelic snp, remove duplicates and normalize snp with long base (bcftools norm -d none/-m-snps), also remove monomorphic snps

#--mask-file ^regions.bed
## the ^ tells bcftools to mask everything outside those regions.
#--mask-file regions.bed
## it will mask everything inside those regions.

#conda activate variant_calling_mapping
cd /home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/vcf/concated_vcf_each_species_REF
output_SFS_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/SFS_data

Andhae_New_REF_AndHae_VCF_and_sfs(){
#for depth in {234,390,546}
for depth in {78,234,390,546,780}
do

    ## output vcf at 1, 3, 5, 7, 10x
    echo -e "\n$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndHae_mask_region $Andhae_New_REF_AndHae_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndHae | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
    -Oz -o ./"$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## index vcf for next step
    ## index vcf for next step
    echo "index vcf file"
    bcftools index -f "$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## output sfs at 1, 3, 5, 7, 10x
    bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' "$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andhae_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

done
}

## Andhae_New_REF_AndFul_VCF
Andhae_New_REF_AndFul_VCF_sfs(){
for depth in {234,390,546}
do

    ## output sfs data
    echo -e "\n$Andhae_New_REF_AndFul_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndFul_mask_region $Andhae_New_REF_AndFul_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndFul | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andhae_New_REF_AndFul_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andhae_New_REF_AndFul_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andhae_New_REF_AndFul_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andhae_New_REF_AndHat_VCF
Andhae_New_REF_AndHat_VCF_sfs(){
for depth in {234,390,546}
do
    ## output sfs data
    echo -e "\n$Andhae_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndHat_mask_region $Andhae_New_REF_AndHat_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndHat | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andhae_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andhae_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andhae_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andhae_New_REF_BomPas_VCF ****************
Andhae_New_REF_BomPas_VCF_sfs(){
MAX_BIN=39
for depth in {546,}
#for depth in {234,390,546}
do
    ## output sfs data
    echo -e "\n$Andhae_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomPas_mask_region $Andhae_New_REF_BomPas_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomPas | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk -v max="$MAX_BIN" '{
        count[$2] = $1
        if ($2 > max) max = $2
    }
    END {
        for (i = 1; i <= max; i++)
        printf "%d ", (i in count ? count[i] : 0)
    }
    END { print "" }' > $output_SFS_dir/"$Andhae_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    ## not deal with sfs with 0
    #awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andhae_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andhae_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andhae_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}


## Andmar_New_REF_AndMar_VCF
#================================= vcf and sfs at each depth 1, 3, 5, 7, 10x for genome subset
Andmar_New_REF_AndMar_VCF_and_sfs(){
for depth in {80,240,400,560,800}
do
    ## output vcf at 1, 3, 5, 7, 10x
    echo -e "\n$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.new.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndMar_mask_region $Andmar_New_REF_AndMar_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndMar | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
    -Oz -o ./"$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## index vcf for next step
    echo "index vcf file"
    bcftools index -f "$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## output sfs at 1, 3, 5, 7, 10x
    bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' "$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.new.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andmar_New_REF_AndMar_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.new.sfs

done
}

## Andmar_New_REF_AndBic_VCF
Andmar_New_REF_AndBic_VCF_sfs(){
for depth in {240,400,560}
do
    ## output sfs data
    echo -e "\n$Andmar_New_REF_AndBic_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndBic_mask_region $Andmar_New_REF_AndBic_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndBic | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andmar_New_REF_AndBic_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andmar_New_REF_AndBic_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andmar_New_REF_AndBic_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andmar_New_REF_AndHae_VCF
Andmar_New_REF_AndHae_VCF_sfs(){
for depth in {240,400,560}
do
    ## output sfs data
    echo -e "\n$Andmar_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndHae_mask_region $Andmar_New_REF_AndHae_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndHae | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andmar_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andmar_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andmar_New_REF_AndHae_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andmar_New_REF_AndHat_VCF
Andmar_New_REF_AndHat_VCF_sfs(){
for depth in {240,400,560}
do
    ## output sfs data
    echo -e "\n$Andmar_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_AndHat_mask_region $Andmar_New_REF_AndHat_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndHat | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andmar_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andmar_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andmar_New_REF_AndHat_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andmar_New_REF_AndTri_VCF
Andmar_New_REF_AndTri_VCF(){
for depth in {240,400,560}
do
bcftools filter --soft-filter mask --mask-file $New_REF_AndTri_mask_region $Andmar_New_REF_AndTri_VCF | \
bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_AndTri | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
bcftools filter -e 'AC==0 || AC == AN' | \
bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
-Oz -o ./"$Andmar_New_REF_AndTri_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Andmar_New_REF_BomPas_VCF
Andmar_New_REF_BomPas_VCF_sfs(){
MAX_BIN=40
for depth in {560,}
#for depth in {240,400,560}
do
    ## output sfs data
    echo -e "\n$Andmar_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomPas_mask_region $Andmar_New_REF_BomPas_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomPas | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk -v max="$MAX_BIN" '{
        count[$2] = $1
        if ($2 > max) max = $2
    }
    END {
        for (i = 1; i <= max; i++)
        printf "%d ", (i in count ? count[i] : 0)
    }
    END { print "" }' > $output_SFS_dir/"$Andmar_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    
    ## not deal with sfs with 0
    ##awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Andmar_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Andmar_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Andmar_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

#output_SFS_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/SFS_data
#================================= vcf and sfs at each depth 1, 3, 5, 7, 10x for genome subset
Bompas_New_REF_BomPas_VCF_and_sfs(){
for depth in {68,204,340,476,680}
do
    ## output vcf at 1, 3, 5, 7, 10x
    echo -e "\n$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomPas_mask_region $Bompas_New_REF_BomPas_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomPas | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
    -Oz -o ./"$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
    ## Above to save corresponding vcf files (from "bcftools view -e")
    ## BomPas_New_REF_BomPas.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP68_1500x_noMS.vcf.gz

    ## index vcf for next step
    echo "index vcf file"
    bcftools index -f "$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## output sfs at 1, 3, 5, 7, 10x
    bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' "$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    ## $1=$1 is an assignment operation in awk. 
    ## It reassigns the first field ($1) to itself, which forces awk to rebuild the current line using the default output field separator (a single space).
    
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bompas_New_REF_BomPas_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    ## outout sfs file names
    ## BomPas_New_REF_BomPas.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP68_1500x_noMS.equal_self.sfs
done
}

## Bompas_New_REF_BomCon_VCF
Bompas_New_REF_BomCon_VCF_sfs(){
for depth in {204,340,476}
do
    ## output sfs data
    echo -e "\n$Bompas_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomCon_mask_region $Bompas_New_REF_BomCon_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomCon | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bompas_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bompas_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bompas_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Bompas_New_REF_BomHor_VCF
Bompas_New_REF_BomHor_VCF_sfs(){
for depth in {204,340,476}
do
    ## output sfs data
    echo -e "\n$Bompas_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomHor_mask_region $Bompas_New_REF_BomHor_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomHor | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bompas_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bompas_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bompas_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Bompas_New_alt_REF_BomMus_VCF
Bompas_New_alt_REF_BomMus_VCF_sfs(){
for depth in {204,340,476}
do
    ## output sfs data
    echo -e "\n$Bompas_New_alt_REF_BomMus_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_alt_BomMus_mask_region $Bompas_New_alt_REF_BomMus_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomMus | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bompas_New_alt_REF_BomMus_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bompas_New_alt_REF_BomMus_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bompas_New_alt_REF_BomMus_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Bompas_New_REF_BomHyp_VCF
Bompas_New_REF_BomHyp_VCF(){
for depth in {204,340,476}
do
bcftools filter --soft-filter mask --mask-file $New_REF_BomHyp_mask_region $Bompas_New_REF_BomHyp_VCF | \
bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomHyp | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
bcftools filter -e 'AC==0 || AC == AN' | \
bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
-Oz -o ./"$Bompas_New_REF_BomHyp_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Bompas_New_REF_ApisMel_VCF
Bompas_New_REF_ApisMel_VCF_sfs(){
MAX_BIN=34
#for depth in {204,340,476}
for depth in {476,}
do
    ## output sfs data
    echo -e "\n$Bompas_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_ApisMel_mask_region $Bompas_New_REF_ApisMel_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_ApisMel | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk -v max="$MAX_BIN" '{
        count[$2] = $1
        if ($2 > max) max = $2
    }
    END {
        for (i = 1; i <= max; i++)
        printf "%d ", (i in count ? count[i] : 0)
    }
    END { print "" }' > $output_SFS_dir/"$Bompas_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    
    ## not deal with sfs with 0
    #awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bompas_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bompas_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bompas_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

#================================= vcf and sfs at each depth 1, 3, 5, 7, 10x for genome subset
Bomvet_New_REF_BomVet_VCF_and_sfs(){
for depth in {58,174,290,416,580}
do
    ## output vcf at 1, 3, 5, 7, 10x
    echo -e "\n$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomVet_mask_region $Bomvet_New_REF_BomVet_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomVet | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
    -Oz -o ./"$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## index vcf for next step
    echo "index vcf file"
    bcftools index -f "$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz

    ## output sfs at 1, 3, 5, 7, 10x
    bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' "$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bomvet_New_REF_BomVet_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
done
}

## 17. Bomvet_New_REF_BomCon_VCF
Bomvet_New_REF_BomCon_VCF_sfs(){
for depth in {174,290,416}
do
    ## output sfs data
    echo -e "\n$Bomvet_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomCon_mask_region $Bomvet_New_REF_BomCon_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomCon | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bomvet_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bomvet_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bomvet_New_REF_BomCon_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## 19. Bomvet_New_REF_BomHor_VCF
Bomvet_New_REF_BomHor_VCF_sfs(){
for depth in {174,290,416}
do
    ## output sfs data
    echo -e "\n$Bomvet_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_BomHor_mask_region $Bomvet_New_REF_BomHor_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_BomHor | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bomvet_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bomvet_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bomvet_New_REF_BomHor_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

## Bomvet_New_REF_ApisMel_VCF ****************
Bomvet_New_REF_ApisMel_VCF_sfs(){
MAX_BIN=29
for depth in {174,290,416}
do
    ## output sfs data
    echo -e "\n$Bomvet_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs"\n"

    bcftools filter --soft-filter mask --mask-file $New_REF_ApisMel_mask_region $Bomvet_New_REF_ApisMel_VCF | \
    bcftools filter --SnpGap 5:indel | bcftools norm -a -f $REF_ApisMel | bcftools view -v snps -A -m 2 -M 2 -f 'PASS' | \
    bcftools filter -e 'AC==0 || AC == AN' | \
    bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP' | \
    awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
    awk -v max="$MAX_BIN" '{
        count[$2] = $1
        if ($2 > max) max = $2
    }
    END {
        for (i = 1; i <= max; i++)
        printf "%d ", (i in count ? count[i] : 0)
    }
    END { print "" }' > $output_SFS_dir/"$Bomvet_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
    
    ## not deal with sfs with 0
    #awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$Bomvet_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs
 
    ## print the sfs in terminal
    cat $output_SFS_dir/"$Bomvet_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.equal_self.sfs

#bcftools view -e "MEAN(FMT/DP) < $depth || MEAN(FMT/DP) > 1500" \
#-Oz -o ./"$Bomvet_New_REF_ApisMel_VCF_filter".SNP_softmask_genic_bi_FMT_DP"$depth"_1500x_noMS.vcf.gz
done
}

#*********  function list of array jobs   *************
function_list=(

"Andhae_New_REF_AndHae_VCF_and_sfs"
"Andhae_New_REF_AndHat_VCF_sfs"
"Andhae_New_REF_AndFul_VCF_sfs"
"Andhae_New_REF_BomPas_VCF_sfs"
"Andmar_New_REF_AndMar_VCF_and_sfs"
"Andmar_New_REF_AndBic_VCF_sfs"
"Andmar_New_REF_AndHat_VCF_sfs"
"Andmar_New_REF_AndHae_VCF_sfs"
"Andmar_New_REF_BomPas_VCF_sfs"
"Bompas_New_REF_BomPas_VCF_and_sfs"
"Bompas_New_alt_REF_BomMus_VCF_sfs"
"Bompas_New_REF_BomHor_VCF_sfs"
"Bompas_New_REF_BomCon_VCF_sfs"
"Bompas_New_REF_ApisMel_VCF_sfs"
"Bomvet_New_REF_BomVet_VCF_and_sfs"
"Bomvet_New_REF_BomPas_VCF_sfs"
"Bomvet_New_REF_BomHor_VCF_sfs"
"Bomvet_New_REF_BomCon_VCF_sfs"
"Bomvet_New_REF_ApisMel_VCF_sfs"

)

## make array jobs
function_array=$(echo ${function_list[*]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)

$function_array


exit 0
