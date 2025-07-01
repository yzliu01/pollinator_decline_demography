
conda activate bedtools

ref_fai_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome
cd $ref_fai_dir

## fai file array variable
genome_fai=(
    "Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai"
    "Andrena_marginata_GCA_963932335.1-softmasked.fa.fai"
    "Bombus_pascuorum-GCA_905332965.1-softmasked.fa.fai"
    "Bombus_veteranus.hifi_asm_pl2.fa.fai"
)

for ref_fai in ${genome_fai[*]}
do
echo $ref_fai
#done
#Steps to Generate a Non-Overlapping BED File for a Given Genome Proportion

#    Get the Total Genome Length
#    Extract the total length from the .fai index of your reference genome:

#awk '{sum+=$2} END {print sum}' $ref_fai
# 330670691

#Compute Target Subset Length
#If you want, say, 10% of the genome:
prop=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
prop_name=(01 02 03 04 05 06 07 08 09)

# prepare full genome prop separately, win_whole
#rop=(1)
#prop_name=(10)

#Divide the Genome into Non-Overlapping Windows
#Use bedtools makewindows to generate evenly spaced non-overlapping regions:
##Adjust -w 100 to change window size to 100bp per fragment.
#bedtools makewindows -g $ref_fai -w 100 > ./random_prop_sample_genome/$ref_fai.win_100b.bed
## example
#1       0       100
#1       100     200
#1       200     300
#1       300     400

#Randomly Sample from Non-Overlapping Windows
#Ensure the selected windows sum up to your target length:
## shuf is produce non-overlapping regions
#shuf ./random_prop_sample_genome/$ref_fai.win_100b.bed | awk -v target="$target_length" '{

## example
#4       7761200 7761300
#4       4028500 4028600
#1       28858800        28858900
#6       6189500 6189600
#7       27870600        27870700

for i in ${!prop[*]}
do
echo "prop: ${prop[i]}"
#done
total_length=$(awk '{sum+=$2} END {print sum}' $ref_fai)
# calculate target genome length
target_length=$(echo "$total_length * ${prop[i]}" | bc)
echo -e "Target genome length: $ref_fai \t $target_length"
#Target genome length: 33067069

# in one step
## for whole genome random due to shuf
## check the sum of ranges, and take the range when in the required sizes
bedtools makewindows -g $ref_fai -w 100 | awk -v target="$target_length" '{
    sum += $3 - $2;
    print $0;
    if (sum >= target) exit;
}' | sort -V | bedtools merge > ./random_prop_sample_genome/$ref_fai.win_100b.shuf_subset_"${prop_name[i]}".sort.bed

## for whole genome linear sampling
#}' | sort -V > ./random_prop_sample_genome/$ref_fai.win_100b.subset_"${prop_name[i]}".bed

## each region
# | sort -V > ./random_prop_sample_genome/$ref_fai.win_100b.shuf_subset_"${prop_name[i]}".sort.bed
done
done




## then downsample reads using bam files



## optionally to check if overlapping exists
random_regions_0_9=Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed
bedtools intersect -a $random_regions_0_9 -b $random_regions_0_9 -wa -wb | awk '$1==$4 && $2<$5'

## merge overlapped regions
dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome/random_prop_sample_genome
bedtools merge -i Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed > prop_09.test.bed
bedtools merge -i Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed | wc -l
1143
wc -l prop_09.test.bed
1143
## original one
wc -l Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed
2976 Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed

## check the length of genome proportion
cat prop_09.test.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

cat random_prop_sample_genome/Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa.fai.prop_0.9.bed \
 | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# 297600000




