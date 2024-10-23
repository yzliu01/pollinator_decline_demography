#!/bin/sh
#SBATCH --account eDNA
#SBATCH --cpus-per-task 20
#SBATCH --mem 1000g
#SBATCH --array=1-4%4
#SBATCH --time=20:30:00
#SBATCH --error=1_bwa_samtools_REF_AndBic_BomSyl_2step_sambamba_markDup_stats.%A_%a.e
#SBATCH --output=1_bwa_samtools_REF_AndBic_BomSyl_2step_sambamba_markDup_stats.%A_%a.o
#SBATCH --job-name=1_bwa_samtools_REF_AndBic_BomSyl_2step_sambamba_markDup_stats

## read fastq files and Read group lines
FASTQ_CLEAN=/faststorage/project/eDNA/yzliu/DK_proj/data/bee_proj_data/fastq_clean
cd $FASTQ_CLEAN
# Andhae_fastq1.fq.clean.gz
# Andhae_fastq2.fq.clean.gz
# Andmar_fastq1.fq.clean.gz
# Andmar_fastq2.fq.clean.gz
seq1=$(ls *fastq1.fq.clean.gz | sed -n ${SLURM_ARRAY_TASK_ID}p) # forward sequence
seq2=$(ls *fastq2.fq.clean.gz | sed -n ${SLURM_ARRAY_TASK_ID}p) # reverse sequence

## Read group
Sorted_ReadGroup_FILE=/faststorage/project/eDNA/yzliu/DK_proj/data/bee_proj_data/bam/ReadGroup_pooled.tab.list
ReadGroup=$(cat $Sorted_ReadGroup_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p)
# @RG\tID:SN8522312001291\tLB:Andhae_pool1\tPL:BGI_DNBSEQ\tPU:gDNA\tSM:Andhae_pool1
# @RG\tID:SN8522312001292\tLB:Andmar_pool1\tPL:BGI_DNBSEQ\tPU:gDNA\tSM:Andmar_pool1
# @RG\tID:SN8522312001289\tLB:Bompas_pool1\tPL:BGI_DNBSEQ\tPU:gDNA\tSM:Bompas_pool1
# @RG\tID:SN8522312001290\tLB:Bomvet_pool1\tPL:BGI_DNBSEQ\tPU:gDNA\tSM:Bomvet_pool1

## mapping output
OUT_BAM=/faststorage/project/eDNA/yzliu/DK_proj/data/bee_proj_data/bam
## reference dir
REF_DIR=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/ref_genome
# Andrena_haemorrhoa-GCA_910592295.1-softmasked.fa
# Andrena_hattorfiana-GCA_944738655.1-softmasked.fa
# Apis_mellifera_HAv-GCF_003254395.2-softmasked.fa
# Bombus_hypnorum-GCA_911387925.1-softmasked.fa
# Bombus_pascuorum-GCA_905332965.1-softmasked.fa

#"Andrena_fulva-GCA_946251845.1-softmasked.fa"
#"Andrena_trimmerana-GCA_951215215.1-softmasked.fa"
#"Bombus_hortorum-GCA_905332935.1-softmasked.fa"
#"Bombus_confusus-GCA_014737475.1_ASM1473747v1-softmasked.fa"

REF1_list=(
            "Andrena_bicolor-GCA_960531205.1.fa"
            "Andrena_bicolor-GCA_960531205.1.fa"
            "Bombus_sylvestris-GCA_911622165.2-softmasked.fa"
            "Bombus_sylvestris-GCA_911622165.2-softmasked.fa"
            )
REF1=$(echo ${REF1_list[*]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)

OUT_NAME_list=(
            "New_REF_AndBic"
            "New_REF_AndBic"
            "New_REF_BomSyl"
            "New_REF_BomSyl"
            )
OUT_NAME=$(echo ${OUT_NAME_list[*]} | tr ' ' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
## out bam file name
File1=${seq1/_fastq1.fq.clean.gz/}

## activate (env) tools of variant_calling_mapping to use sambamba markdup
source /home/yzliu/miniforge3/etc/profile.d/conda.sh
conda activate variant_calling_mapping

## read mapping, convert sam file to bam file, and sort reads
cd $OUT_BAM

## attention to REF variable
bwa mem -t 20 -R $ReadGroup $REF_DIR/$REF1 $FASTQ_CLEAN/$seq1 $FASTQ_CLEAN/$seq2 | samtools sort  -@ 8 -m 100G -o $OUT_BAM/$File1.$OUT_NAME".sort.bam"

## mark duplicates
SORTED_BAM=$File1.$OUT_NAME".sort.bam"
MARKED_BAM=${SORTED_BAM/.sort.bam/.sort.marked_dups.bam}

#picard MarkDuplicates \
#    I=$SORTED_BAM \
#    O=$MARKED_BAM \
#    M=$MARKED_BAM".metrics.csv" >& $MARKED_BAM.log

## less mem intensive
sambamba markdup $SORTED_BAM $MARKED_BAM --nthreads=20 --tmpdir=$TEMPDIR

## index marked_dups
samtools index $MARKED_BAM
## stats
bamtools stats -in $MARKED_BAM > ./bam_stats/bamtools_stats/$MARKED_BAM
samtools stats -in $MARKED_BAM > ./bam_stats/samtools_stats/$MARKED_BAM

exit 0
