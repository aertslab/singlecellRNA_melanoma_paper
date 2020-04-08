#!/bin/bash
### mapping of bulk RNA-Seq data
### count with ht-seq

# fastq files
IN="/ddn1/vol1/staging/leuven/stg_00002/lcb/kdavie/Runs/HiSeq4000_20170609/HiSeq4000_20170609/"
declare -a fastqs=(
MM029_S50_R1_001
MM031_S54_R1_001
MM032_S56_R1_001
MM034_S55_R1_001
MM043_S33_R1_001
MM047_S37_R1_001
)

# output folder
CAT="/staging/leuven/stg_00002/lcb/kspan/Runs/HiSeq4000_20170609/FastQC/"
MAP="/staging/leuven/stg_00002/lcb/kspan/Runs/HiSeq4000_20170609/STAR/"
COUNTS="/staging/leuven/stg_00002/lcb/kspan/Runs/HiSeq4000_20170609/Counts/"
mkdir -p ${CAT}
mkdir -p ${MAP}
mkdir -p ${COUNTS}
cd ${MAP}

# modules
module load STAR/2.5.3a-foss-2014a
module load SAMtools/1.2-foss-2014a
module load FastQC/0.11.5

# quality check
for f in ${fastqs[@]}; do 
	fastqc -t 20 ${IN}${f}.fastq.gz --outdir ${CAT}
done

# genome dictionary
#star indexed genome: StarGenomeIndex.sh
STARGENOME="/staging/leuven/stg_00002/lcb/kdavie/Resources/hg19/STAR_hg19_genome/"

# load mapping index
STAR \
	--runThreadN 20 \
	--genomeDir ${STARGENOME} \
	--genomeLoad LoadAndExit

# mapping & indexing
for f in ${fastqs[@]}; do
s=$(echo $f | perl -pe 's/_S\d+_R1_001//g');
STAR \
	--runThreadN 20 \
	--genomeDir ${STARGENOME} \
	--genomeLoad LoadAndKeep \
	--limitBAMsortRAM 50000000000 \
	--readFilesIn ${IN}${f}.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix ${MAP}${s} \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outWigType bedGraph \
	--outWigStrand Unstranded \
	--outWigNorm RPM \
	--outReadsUnmapped Fastx
samtools index ${MAP}${s}.bam

# filter for high quality reads
samtools view -h -q4 ${MAP}${s}Aligned.sortedByCoord.out.bam |\
       	samtools sort -@ 20 \
	-O bam \
	-T ${MAP}${s}_tmp_sort.bam \
	-o ${MAP}${s}_q4_sorted.bam
done

# count with ht-seq
module load HTSeq/0.6.1p1-foss-2014a-Python-2.7.6
for f in ${fastqs[@]}; do 
s=$(echo $f | perl -pe 's/_S\d+_R1_001//g');
	samtools view ${MAP}${s}_q4_sorted.bam \
		| htseq-count \
			--idattr="gene_name" \
			--minaqual=0 \
			--mode="union" \
			--stranded="no" \
			--type="exon" \
			- \
			/staging/leuven/stg_00002/lcb/resources/human/hg19/gencode.v18.annotation.gtf \
	>${COUNTS}${s}.counts"
done
