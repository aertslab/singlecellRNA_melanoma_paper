# task: link open (ATAC) and accessible (H3K27acetylation) chromatin regions to genes
# using the i-cistarget regions

#module load R/3.6.1-foss-2018a-X11-20180604
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#PATH=$PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/bin/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/lib/
#R
options(stringsAsFactors=FALSE)
.libPaths("/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/3.6")
staging <- "/ddn1/vol1/staging/leuven/stg_00002/lcb"
#############################################################################################
library(data.table)
library(rtracklayer)
library(GenomicRanges)
 
wd <- file.path(staging, 
	"kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/comb4subsets")
setwd(wd)

# function to intersect i-cistarget regions with custom regions
# queryRegions:  GRanges
getIctRegions <- function(queryRegions,  ictRegionsLocationFile, minOverlap=0.4, overlapType="any", returnRegionsOverlap=FALSE, verbose=TRUE) {
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")

  # Region ID location BED file
  ictRegions <- import.bed(con=ictRegionsLocationFile)
  if(verbose) message(paste("Regions in i-cisTarget database:", length(ictRegions)))

  dbRegionsOverlap <- findOverlaps(ictRegions, queryRegions,
                                   maxgap=0L, #minoverlap=1,
                                   type=overlapType, select="all", ignore.strand=TRUE)

  if(minOverlap>0)
  {
    # In i-cisTarget, the default is 40% minimum overlap. Both ways: It takes the maximum percentage (of the peak or the ict region)
    # To reproduce those results:
    overlaps <- pintersect(ictRegions[queryHits(dbRegionsOverlap)], queryRegions[subjectHits(dbRegionsOverlap)])
    percentOverlapHuman <- width(overlaps) / width(ictRegions[queryHits(dbRegionsOverlap)])
    percentOverlapPeaks <- width(overlaps) / width(queryRegions[subjectHits(dbRegionsOverlap)])
    maxOverlap <- apply(cbind(percentOverlapHuman, percentOverlapPeaks), 1, max)
    dbRegionsOverlap <- dbRegionsOverlap[maxOverlap > 0.4]
  }

  # Get regions names
  regionsSelected <- ictRegions[queryHits(dbRegionsOverlap)]@elementMetadata[,"name"]
  regionsSelected <- unique(regionsSelected)

  if(verbose) message(paste("Number of regions selected: ", length(regionsSelected)))

  ret <- regionsSelected
  if(returnRegionsOverlap) {
    regionsOverlap <- data.frame(inputRegion=as.character(queryRegions[subjectHits(dbRegionsOverlap)]),
               iCisTargetRegion=ictRegions[queryHits(dbRegionsOverlap)]@elementMetadata[,"name"])

    regionLocation <- ictRegionsLocation(regionsOverlap$iCisTargetRegion, ictRegionsLocationFile, locationAsString=TRUE)
    regionsOverlap <- data.frame(regionsOverlap, ictRegionLocation=regionLocation[regionsOverlap$iCisTargetRegion,])

    ret <- list(regionsSelected=regionsSelected, regionsOverlap=regionsOverlap)
  }
  return(ret)
}

#i-cistarget regions
ictRegionsLocationFile <- "/staging/leuven/stg_00002/lcb/icistarget/data/regions/homo_sapiens/hg19/refseq-r45/hg19__refseq_r45__ClusteredUniformDHS_all_merge_cleaned2_features_rm-insul_rm-exons2_extend.regionid-location.bed"
gene2regionFile <- "/staging/leuven/stg_00002/lcb/icistarget/data/regions/homo_sapiens/hg19/refseq-r45/hg19__refseq_r45__ClusteredUniformDHS_all_merge_cleaned2_features_rm-insul_rm-exons2_extend-refseq-hg19-ucsc-refgene-table_20kbaroundTSS_HumanRegions.bed_regions_closest_gene_names.regionid-to-geneid.tsv"

# get genes linked to ATAC peaks	
##################################### {
bedDir <-"/staging/leuven/stg_00002/lcb/zkalender/melanoma_ATAC/30.PEAKS/10.HG19/"
cellLines <- list.files(bedDir, pattern="_BL")
	#[1] "MM001_BL"      "MM011_BL"      "MM029_BL"      "MM031_BL"     
	#[5] "MM047_BL"      "MM057_BL"      "MM074_BL_REP1" "MM074_BL_REP2"
	#[9] "MM087_BL" 

peaks_thisDataset <- lapply(cellLines, function(x) import(con=file.path(bedDir, x, paste0(x, "_peaks.narrowPeak")), format="narrowPeak"))

bedDir <- "/staging/leuven/stg_00002/lcb/kspan/melanoma_ATAC/30.PEAKS/10.HG19/"
cellLines <- list.files(bedDir, pattern="_BL")
	#[1] "MM099_BL"
peaks_thisDataset <- append(peaks_thisDataset, import(con=file.path(bedDir, cellLines, paste0(cellLines, "_peaks.narrowPeak")), format="narrowPeak"))

peaks_thisDataset <- do.call(c, peaks_thisDataset)
length(peaks_thisDataset)
#[1] 1714769

peaks_thisDataset <- reduce(peaks_thisDataset)
length(peaks_thisDataset)
#[1] 518833

export(peaks_thisDataset, "all_atac_peaks.bed")

# convert to ictRegions 
thisDatasetRegions <- getIctRegions(peaks_thisDataset, ictRegionsLocationFile, minOverlap=0.4)
	#Regions in i-cisTarget database: 1223024
	#Number of regions selected:  349481

# Intersect with the ones mapping to genes...
gene2regions <- read.table(gene2regionFile, header = FALSE, sep = "\t", comment.char="")
allGenes <- unique(as.character(gene2regions$V2)); length(allGenes)
nrow(gene2regions)
# 220330
gene2regions <- gene2regions[which(gene2regions[,1] %in% thisDatasetRegions),]
nrow(gene2regions)
# 62159
head(gene2regions)
#               V1    V2
#2  chr10-reg68599 LOXL4
#11 chr10-reg68613 LOXL4
#12 chr10-reg68614 LOXL4

atacgenes <- unique(gene2regions$V2)
length(atacgenes)
#[1] 18567

write.table(atacgenes, file = "atac_genes.lst", quote = F, col.names = F, row.names = F)
#only those genes have open atac peak associated
################################################ }

# get genes linked to ATAC peaks intersected with H3K27ac
############################################################## {
bedDir <- file.path(staging, "zkalender/melanoma_ATAC/60.MISC_ANALYSIS")
cellLines <- list.files(bedDir, pattern="_AC_ATAC")
	#[1] "MM001_AC_ATAC" "MM011_AC_ATAC" "MM029_AC_ATAC" "MM031_AC_ATAC"
	#[5] "MM047_AC_ATAC" "MM057_AC_ATAC" "MM074_AC_ATAC" "MM087_AC_ATAC"
	#[9] "MM099_AC_ATAC"
cellLines <- vapply(strsplit(cellLines, "_"), '[', 1, FUN.VALUE=character(1))

atac_ac_peaks <- lapply(cellLines, function(x) import(
	con = file.path(bedDir, paste0(x, "_AC_ATAC"), 
			paste0(x, "_ATAC_with_AC_hg19.bed_500bp.bed")), 
	format="bed"))
atac_ac_peaks <- do.call(c, atac_ac_peaks)
length(atac_ac_peaks)
#[1] 44946
atac_ac_peaks <- reduce(atac_ac_peaks)
length(atac_ac_peaks)
#[1] 15761

# convert to ictRegions 
atac_ac_ict <- getIctRegions(atac_ac_peaks, ictRegionsLocationFile, minOverlap=0.4)

# Intersect with the ones mapping to genes...
gene2regions <- read.table(gene2regionFile, header = FALSE, sep = "\t", comment.char="")
allGenes <- unique(as.character(gene2regions$V2)); length(allGenes)
nrow(gene2regions)
# 220330
gene2regions <- gene2regions[which(gene2regions[,1] %in% atac_ac_ict),]
nrow(gene2regions)
# 5845
head(gene2regions)
#               V1    V2
#2  chr10-reg68599 LOXL4
#11 chr10-reg68613 LOXL4
#12 chr10-reg68614 LOXL4

atac_ac_genes <- unique(gene2regions$V2)
length(atac_ac_genes)
#[1] 5207

write.table(atac_ac_genes, file = "atac_ac_genes.lst", quote = F, col.names = F, row.names = F)
#only those genes have open atac peak associated & H3K27 acetylation
############################################################## }

sessionInfo() #{
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)
#
#Matrix products: default
#BLAS/LAPACK: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
#
#locale:
#[1] C
#
#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets
#[8] methods   base
#
#other attached packages:
#[1] rtracklayer_1.44.4   GenomicRanges_1.36.1 GenomeInfoDb_1.20.0
#[4] IRanges_2.18.3       S4Vectors_0.22.1     BiocGenerics_0.30.0
#[7] data.table_1.12.8
#
#loaded via a namespace (and not attached):
#[1] XVector_0.24.0              zlibbioc_1.30.0
#[3] GenomicAlignments_1.20.1    BiocParallel_1.18.1
#[5] lattice_0.20-38             tools_3.6.1
#[7] SummarizedExperiment_1.14.1 grid_3.6.1
#[9] Biobase_2.44.0              matrixStats_0.55.0
#[11] Matrix_1.2-17               GenomeInfoDbData_1.2.1
#[13] bitops_1.0-6                RCurl_1.95-4.12
#[15] DelayedArray_0.10.0         compiler_3.6.1
#[17] Biostrings_2.52.0           Rsamtools_2.0.1
#[19] XML_3.99-0.3
#}
