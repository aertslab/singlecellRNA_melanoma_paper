#!/bin/bash
### dropseq mapping/counting protocol from nemesh
### for NextSeq run 20170222 cell line MM057 siSOX10
        
#cookbook: http://mccarrolllab.com/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf 

## prerequisites
#################################
# number of beads
NB=37400
# estimated number of cells: 5% of # of beads
NC=$(echo $NB | awk '{printf("%.f\n",$1/20)}')
# estimated number of barcodes = 2x estimated # of cells
NBC=$(echo $NB | awk '{printf("%.f\n",$1/10)}')

# fastq files
IN="/staging/leuven/stg_00002/lcb/kdavie/Runs/NextSeq_20170222_Drop-seq/Demultiplexed/"
FQ1="DropSeq_MM057_Sox10KD_2_S2_R1_001"
FQ2="DropSeq_MM057_Sox10KD_2_S2_R2_001"

# output folder
OUT="/staging/leuven/stg_00002/lcb/kspan/Runs/NextSeq_20170222_Drop-seq/Nemesh_Protocol/"

# sample name (used in output files)
SNAME="DropSeq_MM057_SOX10"

# genome dictionary
#star indexed genome: StarGenomeIndex.sh
STARGENOME="/staging/leuven/stg_00002/lcb/kdavie/Resources/hg19/STAR_hg19_genome/"

#picard indexed genome: picard CreateSequenceDictionary REFERENCE=xy.fasta OUTPUT=xy.dict
GENOME="/staging/leuven/stg_00002/lcb/resources/human/hg19/Homo_sapiens_assembly19_sorted.fa"

# gene annotations
GENES="/staging/leuven/stg_00002/lcb/resources/human/hg19/Homo_sapiens.GRCh37.82.chr-added.gtf"

# modules
module load Picard/1.140-Java-1.7.0_51
DROPSEQTOOLS="/staging/leuven/stg_00002/lcb/kdavie/Drop-seq/Drop-seq_tools-1.12/"
module load R/3.2.1-foss-2014a-x11-tcl


## pre-filter for nextera primer and csp
########################################
module load ea-utils/r819-foss-2014a
fastq-mcf -H -l 20 /staging/leuven/stg_00002/lcb/kspan/fastqmcf.adapters ${IN}${FQ1}.fastq.gz ${IN}${FQ2}.fastq.gz -o ${OUT}${FQ1}.clean.fastq.gz -o ${OUT}${FQ2}.clean.fastq.gz

## run nemesh pipeline
#################################
# [1] combine fwd & rev read in interleaved bam file
#NemeshFastqToSam	

picard FastqToSam \
F1=${OUT}${FQ1}.clean.fastq.gz \
F2=${OUT}${FQ2}.clean.fastq.gz \
O=${OUT}${SNAME}_unmapped.bam \
SAMPLE_NAME=${SNAME}


# [2] extract cell (1st iteration) & molecular (2nd iteration) barcode from every read
# add bam tag to every read, remove low quality reads (NUM_BASES_BELOW_QUALITY)
#NemeshTagBamWithReadSequenceExtended

${DROPSEQTOOLS}TagBamWithReadSequenceExtended \
INPUT=${OUT}${SNAME}_unmapped.bam \
OUTPUT=${OUT}${SNAME}_unmapped_tagged_Cell.bam \
SUMMARY=${OUT}${SNAME}_unmapped_tagged_Cellular.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=2 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

${DROPSEQTOOLS}TagBamWithReadSequenceExtended \
INPUT=${OUT}${SNAME}_unmapped_tagged_Cell.bam \
OUTPUT=${OUT}${SNAME}_unmapped_tagged_CellMolecular.bam \
SUMMARY=${OUT}${SNAME}_unmapped_tagged_Molecular.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=2 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1


# [3] remove reads with low quality barcodes
#NemeshFilterBAM

${DROPSEQTOOLS}FilterBAM \
TAG_REJECT=XQ \
INPUT=${OUT}${SNAME}_unmapped_tagged_CellMolecular.bam \
OUTPUT=${OUT}${SNAME}_unmapped_tagged_filtered.bam


# [4] trim SMART adapter sequences
#NemeshTrimStartingSequence
#In our standard run, we look for at least 5 contiguous bases (NUM_BASES) of the SMART adapter (SEQUENCE) at the 5’ end of the read with no errors (MISMATCHES), and hard clip those bases off the read

${DROPSEQTOOLS}TrimStartingSequence \
INPUT=${OUT}${SNAME}_unmapped_tagged_filtered.bam \
OUTPUT=${OUT}${SNAME}_unmapped_tagged_trimmed_smart.bam \
OUTPUT_SUMMARY=${OUT}${SNAME}_adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5     

#plot trimming report
Rscript -e "trim <- read.table('${OUT}${SNAME}_adapter_trimming_report.txt', skip = 4, header = T); png('${OUT}${SNAME}_adapter_trimming_report.png'); barplot(trim[,2], names.arg=trim[,1]); dev.off()"


# [5] trim trailing polyA tails
#NemeshPolyATrimmer
#searches for at least 6 (NUM_BASES) contiguous A’s in the read with 0 mismatches (MISMATCHES), and hard clips the read to remove these bases and all bases 3’ of the polyA run

${DROPSEQTOOLS}PolyATrimmer \
INPUT=${OUT}${SNAME}_unmapped_tagged_trimmed_smart.bam \
OUTPUT=${OUT}${SNAME}_unmapped_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=${OUT}${SNAME}_polyA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6

#plot trimming report
Rscript -e "trim <- read.table('${OUT}${SNAME}_polyA_trimming_report.txt', skip = 4, header = T); png('${OUT}${SNAME}_polyA_trimming_report.png'); barplot(trim[,2], names.arg=trim[,1]); dev.off()"


# [6] convert to fastq
picard SamToFastq \
INPUT=${OUT}${SNAME}_unmapped_tagged_polyA_filtered.bam \
FASTQ=${OUT}${SNAME}_unmapped_tagged_polyA_filtered.fastq


# [7] align
#NemeshSTAR
module load STAR/2.5.1b-foss-2014a
STAR \
--runThreadN 10 \
--genomeDir $STARGENOME \
--readFilesIn ${OUT}${SNAME}_unmapped_tagged_polyA_filtered.fastq \
--outFileNamePrefix ${OUT}${SNAME}_star


# [8] sort by read name
picard SortSam \
I=${OUT}${SNAME}_starAligned.out.sam \
O=${OUT}${SNAME}_starAligned.srt.bam \
SO=queryname


# [9] merge aligned and unaligned (tagged) reads to recover bam tags
#using only primary alignments
#NemeshMergeBamAlignment

picard MergeBamAlignment \
REFERENCE_SEQUENCE=${GENOME} \
UNMAPPED_BAM=${OUT}${SNAME}_unmapped_tagged_polyA_filtered.bam \
ALIGNED_BAM=${OUT}${SNAME}_starAligned.srt.bam \
OUTPUT=${OUT}${SNAME}_merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false


# [10] tag reads that overlap with exons
#NemeshTagReadWithGeneExon

${DROPSEQTOOLS}TagReadWithGeneExon \
I=${OUT}${SNAME}_merged.bam \
O=${OUT}${SNAME}_star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${GENES} \
TAG=GE

# [11] detect and remove defective cel barcodes
#NemeshDetectBeadSynthesisErrors

${DROPSEQTOOLS}DetectBeadSynthesisErrors \
I=${OUT}${SNAME}_star_gene_exon_tagged.bam \
O=${OUT}${SNAME}_clean.bam \
OUTPUT_STATS=${OUT}${SNAME}_synthesis_stats.txt \
SUMMARY=${OUT}${SNAME}_synthesis_stats.summary.txt \
NUM_BARCODES=${NBC} \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

# [12] count expression
#NemeshDigitalExpression

${DROPSEQTOOLS}DigitalExpression \
I=${OUT}${SNAME}_clean.bam \
O=${OUT}${SNAME}_clean.dge.txt.gz \
SUMMARY=${OUT}${SNAME}_clean.dge.summary.txt \
NUM_CORE_BARCODES=${NC}

# [13] plot the cumulative distribution of reads
#NemeshBAMTagHistogram

${DROPSEQTOOLS}BAMTagHistogram \
I=${OUT}${SNAME}_clean.bam \
O=${OUT}${SNAME}_clean_readcounts.txt.gz \
TAG=XC

Rscript -e "counts <- read.table('${OUT}${SNAME}_clean_readcounts.txt.gz', head = F, stringsAsFactors=F); x <- cumsum(counts[,1]); x <- x/max(x); png('${OUT}${SNAME}_countsperbc.png'); plot(1:length(x), x, type='l', col='blue', xlab='cell barcodes sorted by number of reads [descending]', ylab='cumulative fraction of reads', xlim = c(1,500), lwd = 2); dev.off()"

# [14] delete intermediate files
rm ${OUT}${SNAME}_unmapped.bam \
${OUT}${FQ1}.clean.fastq.gz \
${OUT}${FQ2}.clean.fastq.gz \
${OUT}${SNAME}_unmapped_tagged_Cell.bam \
${OUT}${SNAME}_unmapped_tagged_CellMolecular.bam \
${OUT}${SNAME}_unmapped_tagged_filtered.bam \
${OUT}${SNAME}_unmapped_tagged_trimmed_smart.bam \
${OUT}${SNAME}_unmapped_tagged_polyA_filtered.bam \
${OUT}${SNAME}_unmapped_tagged_polyA_filtered.fastq \
${OUT}${SNAME}_starAligned.out.sam \
${OUT}${SNAME}_starAligned.srt.bam \
${OUT}${SNAME}_merged.bam \
${OUT}${SNAME}_star_gene_exon_tagged.bam \
###${OUT}${SNAME}_clean.bam	#final bam -> delete only after downstream analysis
