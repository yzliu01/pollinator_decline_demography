#!/bin/sh
#SBATCH --account eDNA
##SBATCH --cpus-per-task 20
#SBATCH --mem 50g
#SBATCH --array=1%1
#SBATCH --time=02:00:00
#SBATCH --error=genome_subset_2025_DP_1x_3x_5x_7x_10x.%A_%a.e
#SBATCH --output=genome_subset_2025_DP_1x_3x_5x_7x_10x.%A_%a.o
#SBATCH --job-name=genome_subset_2025_DP_1x_3x_5x_7x_10x
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu2@gmail.com

## activate (env) tools of variant_calling_mapping
source /home/yzliu/miniforge3/etc/profile.d/conda.sh
conda activate variant_calling_mapping

## count line number from no header lines (bcftools view -H)
cd /home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/vcf/concated_vcf_each_species_REF

MAX_BIN=(
#"39"
"40"
#"34"
#"29"
)
species=(
#Hae_39
Mar_40
#Pas_34
#Vet_29
)
#for species in {Hae,Mar}
#for species in Hae
#for species in {Pas,Vet}

for index in ${!species[@]}
do
echo -e "ind_count:${MAX_BIN[index]}\t${species[index]}"

#done

## original vcf files need to be indexed
vcf_list=(
## genome subset with 100% read depth 1,3,5,7,10x
# DP 78,234,390,546,780
#And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.vcf.gz
#And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP234_1500x_noMS.vcf.gz
#And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP390_1500x_noMS.vcf.gz
#And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP546_1500x_noMS.vcf.gz
#And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP780_1500x_noMS.vcf.gz

#DP 80,240,400,560,800
And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP80_1500x_noMS.vcf.gz
And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP240_1500x_noMS.vcf.gz
And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP400_1500x_noMS.vcf.gz
And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP560_1500x_noMS.vcf.gz
And"$species"_New_REF_And"$species".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP800_1500x_noMS.vcf.gz

#DP 68,204,340,476,680
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP68_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP204_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP340_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP476_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP680_1500x_noMS.vcf.gz

#DP 58,174,290,416,580
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP58_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP174_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP290_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP416_1500x_noMS.vcf.gz
#Bom"${species[index]}"_New_REF_Bom"${species[index]}".individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP580_1500x_noMS.vcf.gz

)

#====================================== genome subset snp count =================================================
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


## output sfs directory
output_SFS_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/SFS_data

## test
#vcf_list=(
#BomMaj_REF_BomMaj.DP_1.5x.clip.vcf.gz
#)

## rm BomMaj_REF_BomMaj*g_prop*
genome_subset_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/sofwtare/stairway_plot_v2/stairway_plot_v2.1.2/stairway_plot_blueprint/pool_downsample_genome/templates
#rm $genome_subset_dir/genome_subset.txt

## create summary file before loop
echo -e "out_vcf\tprop_snp_no" > $genome_subset_dir/"${species[index]}"_genome_subset_snp_count.txt

for vcf_file in ${vcf_list[@]}
do

    #for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
    prop_value=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1)
    prop=(01 02 03 04 05 06 07 08 09 10)
    #prop_value=(0.1 0.2)
    #prop=(01 02)

    for i in ${!prop_value[@]}
    do
        #if [[ ${prop[i]} -ne 10 ]]; then
        if [[ ${prop[i]} != 10 ]]; then
        echo $i ${prop_value[i]} ${prop[i]}
        out_vcf_header=${vcf_file/vcf.gz/g_prop_"${prop[i]}".header.vcf}
        out_vcf=${vcf_file/vcf.gz/g_prop_"${prop[i]}".vcf}
        echo $vcf_file $out_vcf
        ##out_vcf
        ##AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.g_prop_01.vcf.gz
        
        ## avoid decimal: | cut -d '.' -f 1
        prop_snp_no=$(echo "$(bcftools view -H $vcf_file | wc -l) * ${prop_value[i]}" | bc | cut -d '.' -f 1)
        ## 1088680 * 0.1
        ## 108868
        prop_site_no=$(echo "$site_count * ${prop_value[i]}" | bc | cut -d '.' -f 1)

        ## tell what is doing
        echo -e "$out_vcf\t$prop_no"
        echo -e "$out_vcf\t$prop_no" >> $genome_subset_dir/"${species[index]}"_genome_subset_snp_count.txt

        ##done

        ##bcftools view $(shuf -n $(echo "$(bcftools view -H $vcf | wc -l) * 0.1 / 1" | bc) input.vcf)
        ## header lines
        bcftools view -h $vcf_file > $out_vcf_header
        ## without header lines
        bcftools view -H $vcf_file | shuf -n $prop_no > $out_vcf
        
        ## output vcf files
        #(cat $out_vcf_header $out_vcf) | bcftools sort -Oz -o $out_vcf.gz

        ## output sfs directly use newly combined vcf files
        ##( ... ) ensures the two cat operations feed one combined stream into bcftools sort
        ## use /tmp but not /tem - bcftools sort often needs a temporary directory, otherwise it may fail for large VCFs.
        ## creat temporary directory
        mkdir -p ./tmp
        (cat $out_vcf_header $out_vcf) | bcftools sort -T ./tmp | bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP\n' | \
        awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
        awk -v max="$MAX_BIN" '{
            count[$2] = $1
            if ($2 > max) max = $2
        }
        END {
            for (i = 1; i <= max; i++)
                printf "%d ", (i in count ? count[i] : 0)
        }
        END { print "" }' > $output_SFS_dir/"$out_vcf".equal_self.sfs

        ##awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$out_vcf".equal_self.sfs
        ##AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.g_prop_01.equal_self.sfs

        ## ddelete temporary vcf files
        rm *clip*g_prop*vcf

        else

        ## avoid decimal: | cut -d '.' -f 1
        prop_no=$(echo "$(bcftools view -H $vcf_file | wc -l) * ${prop_value[i]}" | bc | cut -d '.' -f 1)

        ## when prop=10, use the original vcf directly
        ## vcf_file
        #AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.vcf.gz
        out_vcf=${vcf_file/vcf.gz/g_prop_"${prop[i]}".vcf}

        ## tell what is doing
        echo -e "$out_vcf\t$prop_no"
        echo -e "$out_vcf\t$prop_no" >> $genome_subset_dir/"$species"_genome_subset_snp_count.txt
        
        bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%DP\n' $vcf_file | \
        awk '{if ($3 <= $4/2) print $3; if ($3 > $4/2) print $4-$3 }' | sort -V | uniq -c | \
        awk -v max="$MAX_BIN" '{
            count[$2] = $1
            if ($2 > max) max = $2
        }
        END {
            for (i = 1; i <= max; i++)
                printf "%d ", (i in count ? count[i] : 0)
        }
        END { print "" }' > $output_SFS_dir/"$out_vcf".equal_self.sfs
        ## lack values for 0
        ##awk '$1=$1'| cut -d ' ' -f 1 | tr '\n' ' ' > $output_SFS_dir/"$out_vcf".equal_self.sfs

        ##AndHae_New_REF_AndHae.individual_100kb_1500x_region.clip.SNP_softmask_genic_bi_FMT_DP78_1500x_noMS.g_prop_10.equal_self.sfs

        fi

    done

done

#mv genome_subset.txt Hae_genome_subset.txt

done

exit 0
