# task: aggregate single-cell counts per cell line / treatment / melanoma cell state
# library-size normalize & z-score
# for Fig4, Fig7d, S13 & S14, S22

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
library(robustbase)
library(DESeq2)

# row z-score scaling (center on mean = 'normal' zscore)
zscore <- function(x) {
        #first calculate row (1) means, then subtract mean ("-")
        y <- sweep(x, 1, apply(x, 1, mean), "-")
        #then calculate row standard deviation and devide ("/") by sd
        sweep(y, 1, apply(x, 1, sd), "/")}

zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k")
kdir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k")

load(file.path(zdir, "10.analysis_with_all_data/00.RData/meta_data.RData"))
raw.counts <- read.table(file.path(zdir, "10.analysis_with_all_data/10.TextData/raw_counts.tsv"))
	dim(raw.counts)
	#[1] 32738 43582
	#all cellranger genes, all 40k cells

# remove doublets
load(file.path(kdir, "10.analysis_with_all_data/00.RData/demuxlet_doublets.RData"))
table(demuxlet_doublets %in% colnames(raw.counts))
        #FALSE  TRUE 
        #    8   404
raw.counts <- raw.counts[,!colnames(raw.counts) %in% demuxlet_doublets]
meta_data <- meta_data[!rownames(meta_data) %in% demuxlet_doublets,]
raw.counts <- raw.counts[, rownames(meta_data)]

# [1] aggregate baselines of melanocytic/intermediate/mesenchymal-like w/o A375 (Fig4, S13 & S14) {
samples <- list(melanocytic = c("MM001_BL", "MM011_BL", "MM031_BL"),
	intermediate = c("MM074_BL", "MM087_BL", "MM057_BL"),
	mesenchymal_like = c("MM029_BL", "MM047_BL", "MM099_BL"))

aggr <- sapply(names(samples), function(set) { rowSums(raw.counts[,rownames(meta_data)[meta_data$CE_ID %in% samples[[set]]]]) })

#rm all-zero genes
aggr <- aggr[!rowSums(aggr)==0,]

#add pseudocount
aggr <- aggr + matrix(1, nrow = nrow(aggr), ncol = ncol(aggr))

#normalize
lib.size <- colSums(aggr)
aggr <- t(apply(aggr, 1, function(x) {x/lib.size*10^6}))

#log2
aggr <- log2(aggr)
	
#z-score
z <- zscore(aggr)

write.table(cbind(gene = rownames(z), z), file = "MMgroup_woA375__aggr_norm_scaled.mrna", 
		row.names = F, quote = F, sep = "\t")
#}

# [2] three intermediate lines KD treatments (Fig S22) {
lines <- c("MM074", "MM087", "MM057")
treats <- c("BL", "NTC", "SOX10_24h", "SOX10_48h", "SOX10_72h")
ceids <- sort(as.vector(outer(lines, treats, paste, sep = "_")))

aggr <- sapply(ceids, function(ceid) { 
	       rowSums(raw.counts[,rownames(meta_data)[meta_data$CE_ID == ceid]]) })

#rm all-zero genes
aggr <- aggr[!rowSums(aggr)==0,]

#add pseudocount
aggr <- aggr + matrix(1, nrow = nrow(aggr), ncol = ncol(aggr))

#normalize
lib.size <- colSums(aggr)
aggr <- t(apply(aggr, 1, function(x) {x/lib.size*10^6}))

#log2
aggr <- log2(aggr)
	
#z-score
z <- zscore(aggr)

z <- zscore(aggr[, ceids])
write.table(cbind(gene = rownames(z), z), file = "three_lines_aggr_norm_scaled.mrna",
		row.names = F, quote = F, sep = "\t")
#}

# [3] calculate log2FC for SOX10 treatments (Fig7d) {
lines <- c("MM074", "MM087", "MM057")
treats <- c("BL", "NTC", "SOX10_24h", "SOX10_48h", "SOX10_72h")
ceids <- sort(as.vector(outer(lines, treats, paste, sep = "_")))

aggr <- sapply(ceids, function(ceid) { 
	rowSums(raw.counts[,rownames(meta_data)[meta_data$CE_ID == ceid]])
	})

#rm all-zero genes
aggr <- aggr[!rowSums(aggr)==0,]

design <- data.frame(list(treatment = 
	c(rep("BL", 3), rep("NTC", 3), rep("SOX10_24h", 3), rep("SOX10_48h", 3), rep("SOX10_72h", 3)),
	cellline = c(rep(c("MM057", "MM074", "MM087"), 5))
	))
rownames(design) <- paste(design[, "cellline"], design[, "treatment"], sep = "_")
design <- design[colnames(aggr),]

dds <- DESeqDataSetFromMatrix(countData = aggr, 
			      colData = design,
			      design = formula(~ cellline + treatment))
dds <- DESeq(dds)

results <- list()
for (treat in c("SOX10_24h", "SOX10_48h", "SOX10_72h")) {
	print(treat)
	results[[paste0(treat, "_vs_BL")]] <- results(dds, contrast = c("treatment", treat, "BL"))
	results[[paste0(treat, "_vs_NTC")]] <- results(dds, contrast = c("treatment", treat, "NTC"))
}
saveRDS(results, file = "three_lines_aggr_deseq2.rds.gz")

#for supplementary data write full table w/o replacing padj by TRUE/FALSE
res_ntc <- results[c("SOX10_24h_vs_NTC", "SOX10_48h_vs_NTC","SOX10_72h_vs_NTC")]
names(res_ntc) <- c("SOX10-KD_24h", "SOX10-KD_48h","SOX10-KD_72h")
res_ntc <- lapply(names(res_ntc), function(x) {
	tab <- as.data.frame(res_ntc[[x]])
	colnames(tab) <- paste0(x, "_", colnames(tab))
	return(tab)
})
res_ntc <- do.call("cbind", res_ntc)

write.table(cbind(gene = rownames(res_ntc), res_ntc), 
	    file = "three_lines_aggr_deseq2_SOX10_vs_NTC.txt",
	    row.names = F, quote = F, sep = "\t")
#}

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
#[1] DESeq2_1.24.0               SummarizedExperiment_1.14.1
#[3] DelayedArray_0.10.0         BiocParallel_1.18.1
#[5] matrixStats_0.55.0          Biobase_2.44.0
#[7] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0
#[9] IRanges_2.18.3              S4Vectors_0.22.1
#[11] BiocGenerics_0.30.0         robustbase_0.93-5
#[13] data.table_1.12.8
#
#loaded via a namespace (and not attached):
#[1] bit64_0.9-7            splines_3.6.1          Formula_1.2-3
#[4] assertthat_0.2.1       BiocManager_1.30.4     latticeExtra_0.6-28
#[7] blob_1.2.0             GenomeInfoDbData_1.2.1 pillar_1.4.2
#[10] RSQLite_2.1.2          backports_1.1.5        lattice_0.20-38
#[13] glue_1.3.1             digest_0.6.21          RColorBrewer_1.1-2
#[16] XVector_0.24.0         checkmate_1.9.4        colorspace_1.4-1
#[19] htmltools_0.3.6        Matrix_1.2-17          XML_3.99-0.3
#[22] pkgconfig_2.0.3        genefilter_1.66.0      zlibbioc_1.30.0
#[25] purrr_0.3.2            xtable_1.8-4           scales_1.0.0
#[28] htmlTable_1.13.1       tibble_2.1.3           annotate_1.62.0
#[31] ggplot2_3.2.1          nnet_7.3-12            lazyeval_0.2.2
#[34] survival_2.44-1.1      magrittr_1.5           crayon_1.3.4
#[37] memoise_1.1.0          foreign_0.8-72         tools_3.6.1
#[40] stringr_1.4.0          locfit_1.5-9.1         munsell_0.5.0
#[43] cluster_2.1.0          AnnotationDbi_1.46.1   compiler_3.6.1
#[46] rlang_0.4.0            grid_3.6.1             RCurl_1.95-4.12
#[49] rstudioapi_0.10        htmlwidgets_1.3        bitops_1.0-6
#[52] base64enc_0.1-3        gtable_0.3.0           DBI_1.0.0
#[55] R6_2.4.1               gridExtra_2.3          knitr_1.25
#[58] dplyr_0.8.3            bit_1.1-15.2           zeallot_0.1.0
#[61] Hmisc_4.2-0            stringi_1.4.3          Rcpp_1.0.3
#[64] vctrs_0.2.0            geneplotter_1.62.0     rpart_4.1-15
#[67] acepack_1.4.1          DEoptimR_1.0-8         tidyselect_0.2.5
#[70] xfun_0.10
#}
