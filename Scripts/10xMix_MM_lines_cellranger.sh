#!/bin/bash
### cellranger mapping/counting pipeline
### for NextSeq run 20170420 10x experiment mix of MM lines
### update 20170508: rerun for combined sequencing runs NextSeq_20170420_10x & HiSeq4000_20170504


## cellranger prerequisites
#################################
# fastq files
IN="/staging/leuven/stg_00002/lcb/kspan/10xCombined/"
FQ="Mix_MM_lines"

# output folder
OUT="/staging/leuven/stg_00002/lcb/kspan/10xCombined/"
cd ${OUT}

# sample name (used in output files)
SNAME="Mix_MM_lines"

# modules
module load cellranger/1.3.1-foss-2014a
module load R/3.2.1-foss-2014a-x11-tcl

## run cellranger
#################################
# concatenate fastq files from this and previous run
declare -a a=(R1 R2)
for i in ${a[@]}; do
	echo $i;
	cat /ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/Runs/HiSeq4000_20170504/HJ5J5BBXX/outs/fastq_path/HJ5J5BBXX/CatFastq/Mix_MM_lines_S1_L007_${i}_001.fastq.gz /staging/leuven/stg_00002/lcb/kspan/Runs/NextSeq_20170420_10x/Mix_MM_lines_S1_L001_${i}_001.fastq.gz >${OUT}Mix_MM_lines_S1_L001_${i}_001.fastq.gz;
done
 
# cellranger
echo -e "$(date)\nrunning cellranger count ..."
cellranger count \
--id=${SNAME} \
--cells=5000 \
--fastqs=${OUT} \
--sample=${FQ} \
--transcriptome=/staging/leuven/stg_00002/lcb/kdavie/Resources/hg19/CellRanger/refdata-cellranger-hg19-1.2.0/ \
--jobmode=local \
--localcores=20

# convert gene expression matrix to tsv
echo -e "$(date)\nrunning\n\tRscript ~/src_kspan/10xCountMatrices.R ${OUT} ${SNAME}\n..."
Rscript ~/src_kspan/10xCountMatrices.R ${OUT} ${SNAME}

# gene body coverage
echo -e "$(date)\nrunning\n\t/user/leuven/306/vsc30697/src_kspan/calc_genebodycoverage_10x.sh ${OUT}${SNAME}\n..."
module purge
/user/leuven/306/vsc30697/src_kspan/calc_genebodycoverage_10x.sh ${OUT}${SNAME} >>${OUT}${SNAME}/gbc.log 2>&1
