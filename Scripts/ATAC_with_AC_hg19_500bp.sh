list='MM011 MM029 MM031 MM047 MM057 MM074 MM087 MM099' 
for sample in $list; do


mkdir -p /staging/leuven/stg_00002/lcb/zkalender/melanoma_ATAC/60.MISC_ANALYSIS/${sample}_AC_ATAC

cd /staging/leuven/stg_00002/lcb/zkalender/melanoma_ATAC/60.MISC_ANALYSIS/${sample}_AC_ATAC

# find ATAC peaks without AC signal 
dir=/ddn1/vol1/staging/leuven/stg_00002/lcb/lvanaers/MM_ATAC/Paper_MMlines/20.MM_ATAC_H3K27ac/hg38/PIM_4kb

# select top 5000 peaks based on ATAC peak-score
cat $dir/${sample}_BL_non-overlap_ATAC-H3K27ac.bed | sort -k5 -rg | awk 'NR<=5000' | cut -f1-3 > ${sample}_ATAC_wo_AC_hg38.bed
module load bedtools/20170824-foss-2014a 

# the peaks are extended to 4kb by Linde, shrink it back to 500bp peaks
TARGET_LENGTH=500
 awk -vF=${TARGET_LENGTH} 'BEGIN{ OFS="\t"; } { len=$3-$2; diff=F-len; flank=int(diff/2); upflank=downflank=flank; if (diff%2==1) { downflank++; }; print $1, $2-upflank, $3+downflank; }' ${sample}_ATAC_wo_AC_hg38.bed | sortBed -i stdin > ${sample}_ATAC_wo_AC_hg38.bed_500bp.bed

# liftover from hg38 to hg19 for i-cisTarget analysis
/data/leuven/software/biomed/Kent/20170816/bin/liftOver ${sample}_ATAC_wo_AC_hg38.bed_500bp.bed ~/lcb/zkalender/resources/hg38ToHg19.over.chain.gz  ${sample}_ATAC_wo_AC_hg19.bed_500bp.bed unmapp

# find ATAC peaks with AC signal 
cat $dir/${sample}_PIM_ATAC-H3K27ac.bed | awk '$6>0' | sort -k5 -rg | awk 'NR<=5000' | cut -f1-3 > ${sample}_ATAC_with_AC_hg38.bed
module load bedtools/20170824-foss-2014a 

# the peaks are extended to 4kb by Linde, shrink it back to 500bp peaks
TARGET_LENGTH=500
 awk -vF=${TARGET_LENGTH} 'BEGIN{ OFS="\t"; } { len=$3-$2; diff=F-len; flank=int(diff/2); upflank=downflank=flank; if (diff%2==1) { downflank++; }; print $1, $2-upflank, $3+downflank; }' ${sample}_ATAC_with_AC_hg38.bed | sortBed -i stdin > ${sample}_ATAC_with_AC_hg38.bed_500bp.bed

# liftover from hg38 to hg19 for i-cisTarget analysis
/data/leuven/software/biomed/Kent/20170816/bin/liftOver ${sample}_ATAC_with_AC_hg38.bed_500bp.bed ~/lcb/zkalender/resources/hg38ToHg19.over.chain.gz  ${sample}_ATAC_with_AC_hg19.bed_500bp.bed unmapp 

done
