#!/bin/bash
# Pipeline: Recessive + De novo LOF variants (AF < 0.05)

set -euo pipefail

# ==================== CONFIGURATION (EDIT THESE) ====================
PROJECT_DIR=/path/to/your/project          
DATA_DIR=/path/to/data          
CHR_LIST=$(ls $PROJECT_DIR/data/ | grep tbi | awk -F. '{print $2}' | grep -v X | tail -n 28)

# Subdirectories (relative to PROJECT_DIR)
SAMPLE_LIST="$PROJECT_DIR/AC/data_sample/AC_sample.list"
COLNAME_FILE="$PROJECT_DIR/AC/data_sample/colname"
FAMILY_INPUT="$PROJECT_DIR/AC/data_sample/AC_family.input"
FILTER_RESULT_DIR="$PROJECT_DIR/randomf/filter_result"
GNOMAD_INPUT_DIR="$PROJECT_DIR/AC/get_AF/input"
DENOVO_REC_RESULT_DIR="$PROJECT_DIR/AC/denovo_rec/result"
PASS_VAR_DIR="$PROJECT_DIR/pass_var"
RAREVAR_DIR="$PROJECT_DIR/rarevar"
ANNO_DIR="$PROJECT_DIR/anno"
AC_DATA_VCF_DIR="$PROJECT_DIR/AC/data_vcf"
AC_DATA_GT_DIR="$PROJECT_DIR/AC/data_GT"
AC_PASS_VID_DIR="$PROJECT_DIR/AC/pass_vid"
AC_RAREVAR_DIR="$PROJECT_DIR/AC/rarevar"
PYTHON_SCRIPT="$PROJECT_DIR/AC/denovo_rec/GT_denovo_rec.py"
DNM_PASS_LIST="$PROJECT_DIR/AC/denovo/vid_dnm.pass.list"

# ==================== 1. Extract genotypes ====================
for CHR in $CHR_LIST; do
    bcftools view -S $SAMPLE_LIST \
        $DATA_DIR/vcf/ACH_CHD_GRCh37.${CHR}.dn.vcf.gz \
        -c1 -Oz -o $AC_DATA_VCF_DIR/AC_GRCh37.${CHR}.dn.nofilter.vcf.gz
    tabix -p vcf $AC_DATA_VCF_DIR/AC_GRCh37.${CHR}.dn.nofilter.vcf.gz

    bcftools query -f "%ID[\t%GT]\n" $AC_DATA_VCF_DIR/AC_GRCh37.${CHR}.dn.nofilter.vcf.gz | \
        grep -v '\*' | cat $COLNAME_FILE - > $AC_DATA_GT_DIR/AC_GRCh37.${CHR}.GT.nofilter.txt
done

# ==================== 2. Keep QC-passed variants ====================
for CHR in $CHR_LIST; do
    awk '{if($2==1) print $1}' $FILTER_RESULT_DIR/${CHR}.tsv > $PASS_VAR_DIR/${CHR}_pass_vid.tsv
done

for CHR in $CHR_LIST; do
    awk 'NR==FNR{a[$1]=$0;} NR>FNR && a[$1] {print a[$1]"\t"$0}' \
        $DENOVO_REC_RESULT_DIR/AC_GRCh37.${CHR}_rec_denovo.nofilter.txt \
        $PASS_VAR_DIR/${CHR}_pass_vid.tsv | \
        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $AC_PASS_VID_DIR/AC_GRCh37.${CHR}_rec_denovo.pass.txt
done

# ==================== 3. Keep rare variants (AF < 0.05) ====================
for CHR in $CHR_LIST; do
    awk '{if($5 < 0.05) print $0}' $GNOMAD_INPUT_DIR/ACH_CHD_GRCh37.${CHR}.gnomad.txt > $RAREVAR_DIR/AC_CHD_${CHR}_rarevar.txt
done

for CHR in $CHR_LIST; do
    awk 'NR==FNR{a[$1]=$0;} NR>FNR && a[$1] {print a[$1]"\t"$0}' \
        $AC_PASS_VID_DIR/AC_GRCh37.${CHR}_rec_denovo.pass.txt \
        $RAREVAR_DIR/AC_CHD_${CHR}_rarevar.txt > $AC_RAREVAR_DIR/AC_CHD_${CHR}_rec_denovo.pass_rarevar.txt
done

# ==================== 4. Annotate (exclude intergenic) ====================
for CHR in $CHR_LIST; do
    bcftools query -f '%ID\t%INFO/CSQ' $DATA_DIR/vcfqc/ACH_CHD_GRCh37.${CHR}.anno3.vcf.gz | \
        grep -v intergenic_variant | \
        awk -F'|' '{print $1"\t"$2"\t"$3"\t"$4}' > $ANNO_DIR/${CHR}.anno.txt
done

# Merge rare variants with annotations (adjust column numbers if needed)
for CHR in $CHR_LIST; do
    awk 'NR==FNR{a[$1]=$0;} NR>FNR && a[$1] {print a[$1]"\t"$0}' \
        $AC_RAREVAR_DIR/AC_CHD_${CHR}_rec_denovo.pass_rarevar.txt \
        $ANNO_DIR/${CHR}.anno.txt | \
        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$12"\t"$15"\t"$16"\t"$17}' >> $PROJECT_DIR/AC_all_rare_AF0.05_annotated.tsv
done

# ==================== 5. Identify recessive & de novo patterns (Python) ====================
for CHR in $CHR_LIST; do
    python $PYTHON_SCRIPT \
        $FAMILY_INPUT \
        $AC_DATA_GT_DIR/AC_GRCh37.${CHR}.GT.nofilter.txt \
        $DENOVO_REC_RESULT_DIR/AC_GRCh37.${CHR}_rec_denovo.nofilter.txt
done
# Merge de novo and recessive pass list with annotations
for CHR in $CHR_LIST; do
    awk 'NR==FNR{a[$1]=$0;} NR>FNR && a[$1] {print a[$1]"\t"$0}' \
        $DENOVO_REC_RESULT_DIR/AC_GRCh37.${CHR}_rec_denovo.nofilter.txt \
        $PROJECT_DIR/AC_all_rare_AF0.05_annotated.tsv >> $PROJECT_DIR/AC_rec.pass_coding_AF0.05.tsv
done

awk 'NR==FNR{a[$1]=$0;} NR>FNR && a[$1] {print a[$1]"\t"$0}' \
    $DNM_PASS_LIST \
    $PROJECT_DIR/AC_all_rare_AF0.05_annotated.tsv > $PROJECT_DIR/AC_dnm.pass_coding_AF0.05.tsv

# ==================== 6. Filter LOF variants (frameshift/stop_gain) ====================
LOF_OUTPUT="$PROJECT_DIR/AC_LOF_rare_AF0.05_rec_denovo.tsv"
# Column assumptions: $2 = rec (recessive), $9 = Consequence (adjust if needed)
awk 'BEGIN{OFS="\t"; print "#LOF_variants\tAF<0.05\trecessive_or_denovo"} 
     NR==1 {print $0, "LOF_type"} 
     NR>1 {
         lof=0; reason="";
         if($9 ~ /frameshift_variant/) {lof=1; reason="frameshift"}
         if($9 ~ /stop_gained/) {lof=1; reason="stop_gain"}
         # Recessive: $2 not empty and not ";"
         # De novo variants are already in this file (from DNM list)
         if(lof==1 && $2 !~ /^;$/ && $2 != "") print $0, reason;
     }' $PROJECT_DIR/AC_dnm.pass_coding_AF0.05.tsv > $LOF_OUTPUT

awk 'BEGIN{OFS="\t"; print "#LOF_variants\tAF<0.05\trecessive_or_denovo"} 
     NR==1 {print $0, "LOF_type"} 
     NR>1 {
         lof=0; reason="";
         if($9 ~ /frameshift_variant/) {lof=1; reason="frameshift"}
         if($9 ~ /stop_gained/) {lof=1; reason="stop_gain"}
         # Recessive: $2 not empty and not ";"
         # De novo variants are already in this file (from DNM list)
         if(lof==1 && $2 !~ /^;$/ && $2 != "") print $0, reason;
     }' $PROJECT_DIR/AC_rec.pass_coding_AF0.05.tsv >> $LOF_OUTPUT
