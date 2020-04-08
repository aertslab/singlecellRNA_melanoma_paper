# task: calculate AUCell scores on motif and track regulons 100 pySCENIC runs subset 10 baselines
# boxplots Fig S4

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
library(AUCell)
library(ggplot2)
library(parallel)

seed <- 42

wd <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL")
setwd(wd)

trkdir <- file.path(wd, "rerunSCENICtracks100x")
mtfdir <- file.path(wd, "100_scenic_runs")

# get AUCell rankings (from script make_AUCell_ranking_db.R)
aucellRankings <- readRDS("aucellRankings.rds.gz")
  #0c00ab9146f1d19227a309e974af008dfdaaee97 aucellRankings.rds.gz

# get meta data
  zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/")
  meta <- read.table(file = file.path(zdir, "20.analysis_with_sub_matrices/12.ten_lines_BL/10.TextData/meta_data.tsv"))
  #rm doublets
  meta <- meta[colnames(aucellRankings@assays$data$ranking), ]
 
# function to get median AUCell value per class
medianAUC <- function(AUC, classes) {
	#classes: vector of classes, names are cell names
	classlist <- split(names(classes), unname(classes))
	do.call(rbind, 
   	  lapply(classlist, function(cells) {
		apply(AUC[cells, ], 2, median)
   	  })
	)
}

# prepare colours for plots
cols <- setNames(c("#0066FFFF", "#800080FF", "#FF5B35FF"),
	c("melanocyte-like", "intermediate", "mesenchymal-like"))
	
# [1] make AUCell matrices of 100x runs {
# track dbs
  db <- "trk"
  mclapply(1:100, function(i) {
  	# get & format regulons
  	regs <- readRDS(file.path(trkdir, paste0("run_", i), "reg_trk.rds.gz"))
	#append "_neg/_pos" to tf name
	for (dir in names(regs)) { names(regs[[dir]]) <- paste0(names(regs[[dir]]), "_", dir) }
	
	#make one list
	regs <- c(regs[[1]], regs[[2]])

  	# calculate AUC
	regulonAUC <- AUCell_calcAUC(regs, aucellRankings, 
		aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 4)
	auc <- data.frame(t(getAUC(regulonAUC)))

	# save
	saveRDS(auc, 
		file = file.path(trkdir, paste0("run_", i), 
		paste0("regulons_extended_", db, "_aucMatrix01.rds.gz")))
  }, mc.cores = 8)

# mtf dbs
  db <- "mtf"
  mclapply(1:100, function(i) {
  	# get & format regulons
  	regs <- readRDS(file.path(mtfdir, paste0("run_", i), "regulons_extended_mtf.rds.gz"))
	#append "_pos" to tf name, for the sake of consistency
	names(regs) <- paste0(names(regs), "_pos")

  	# calculate AUC
	regulonAUC <- AUCell_calcAUC(regs, aucellRankings, 
		aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 4)
	auc <- data.frame(t(getAUC(regulonAUC)))

	# save
	saveRDS(auc, 
		file = file.path(mtfdir, paste0("run_", i), 
		paste0("regulons_extended_", db, "_aucMatrix01.rds.gz")))
  }, mc.cores = 8)
#}
	
# [2] calculate medianAUC per cell state / regulon {
# foreach run in 100 runs:

mclapply(1:100, function(i) {
	auc <- list()
	
	trk <- readRDS(file.path(trkdir, paste0("run_", i), "regulons_extended_trk_aucMatrix01.rds.gz"))
	colnames(trk) <- paste0(colnames(trk), "_trk")
	
	mtf <- readRDS(file.path(mtfdir, paste0("run_", i), "regulons_extended_mtf_aucMatrix01.rds.gz"))
	colnames(mtf) <- paste0(colnames(mtf), "_mtf")
	
	bothdbs <- cbind(trk, mtf)
	
	## compare three states
	classes <- setNames(meta$Cell_Line, rownames(meta))
	mmcells <- rownames(meta[meta$Cell_Line != "A375",])
	classes <- classes[mmcells]
	classes[classes %in% c("MM087", "MM057", "MM074")] <- "intermediate"
	classes[classes %in% c("MM001", "MM011", "MM031")] <- "melanocyte-like"
	classes[classes %in% c("MM029", "MM047", "MM099")] <- "mesenchymal-like"
	auc[["three_states"]] <- medianAUC(bothdbs, classes)
	
	saveRDS(auc, 
		file = file.path(trkdir, paste0("run_", i),
			"medianAUC_mtf_and_new_trk_regulons_aucell01.rds.gz"), 
		compress = "gzip")
}, mc.cores = 8)
#}

# [3] get auc of filtered regulons, calculate median per cell state {
auc_filtregs <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_auc01_auc.rds.gz")
auc_filtregs <- t(data.frame(auc_filtregs))
filtregs <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_auc01_regulons.rds.gz")

# calculate median AUC for three states
classes <- setNames(meta$Cell_Line, rownames(meta))
mmcells <- rownames(meta[meta$Cell_Line != "A375",])
classes <- classes[mmcells]
classes[classes %in% c("MM087", "MM057", "MM074")] <- "intermediate"
classes[classes %in% c("MM001", "MM011", "MM031")] <- "melanocyte-like"
classes[classes %in% c("MM029", "MM047", "MM099")] <- "mesenchymal-like"
mauc_filtregs <- medianAUC(auc_filtregs, classes)
saveRDS(mauc_filtregs, file = "all_pos_regulons_recurrent_100x_jw_new_trk_auc_median_three_states_aucell01.rds.gz")

mauc_filtregs <- reshape2::melt(as.matrix(mauc_filtregs))
setDT(mauc_filtregs)
colnames(mauc_filtregs) <- c("rn", "variable", "value")
#}

# [4] get median AUCell scores 100 runs, subset to regulons that recur >= 80x {
auc_100runs <- readRDS(file = "medianAUC_100runs_three_states_aucell01.rds.gz")
auc_100runs <- auc_100runs[grepl("pos", variable)]
auc_100runs[, variable := gsub("_pos", "", variable)]

# remove non-tf regulons
nontfs <- read.table(file.path(staging, "kspan/resources/hg19/hg19_nonTFs.lst"), header = F)$V1
auc_100runs[, tf := tstrsplit(variable, "_", fixed = TRUE)[1]]
auc_100runs <- auc_100runs[!tf %in% nontfs]

# subset auc to regulons that come up in >= 80 runs
#each regulon is listed x times for each pyscenic run (x: number of groups contrasted)
nlevels <- length(unique(auc_100runs$rn))
cutoff <- 80
recregs <- as.vector(auc_100runs[, .N, by = .(variable)][N >= cutoff*nlevels, variable])
print(paste0("regulons found >= ", cutoff, "x: ", length(recregs)))

auc_100 <- auc_100runs[variable %in% recregs]
#}

# [5] plot AUCell scores 100 runs as boxplots, filtered regulons as dots {
# match regulon names
table(unique(mauc_filtregs$variable) %in% unique(auc_100$variable))
	#FALSE  TRUE 
    	#1   363 
#[1] NKX3-1_mtf
auc_100[variable == "NKX3.1_mtf", variable := "NKX3-1_mtf"]

# combine both auc matrices to calculate zscore
mauc_filtregs[, run := 0]
auc <- rbind(auc_100[, .(rn, variable, value, run)], mauc_filtregs)
auc[, zscore := scale(value), by = variable]
saveRDS(auc, file = "medianAUC_100runs_three_states_regRecurrence80_zscore_plusFiltRegs_aucvalues_aucell01.rds.gz")

pdf(paste0("medianAUC_100runs_", comparison, "_regRecurrence80_zscore_plusFiltRegs_aucell01_meanOrdered.pdf"))
states <- as.vector(unique(auc$rn))
for (state in states) {

	#order by mean
	regorder <- as.vector(auc[run != 0 & rn == line, 
		mean(zscore, na.rm = T), by = variable][order(V1, decreasing = T), variable])
	auc$variable <- factor(auc$variable, levels = regorder, ordered = T)

	#plot top20 only
	p <- ggplot(auc[run != 0 & variable %in% regorder[1:20]], 
			aes(x = variable, y = zscore, fill = rn, colour = rn)) +
		geom_boxplot(na.rm = TRUE, outlier.shape=NA) +

		geom_point(data = auc[run == 0 & variable %in% regorder[1:20]],
			aes(x = variable, y = zscore, fill = rn), shape = 21, colour = "black") +

		scale_colour_manual(values = cols) +
		scale_fill_manual(values = cols) +
		labs(title = paste0("out of ", length(unique(auc$variable)), 
			" recurrent factors top20 (by mean) for ", line),
			x="", y = "zscore(median AUCell)") +
		theme_classic() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			axis.title.x=element_blank(),
			panel.background = element_blank())
	print(p)
}
dev.off()


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
#[1] parallel  stats     graphics  grDevices utils     datasets  methods
#[8] base
#
#other attached packages:
#[1] ggplot2_3.2.1     AUCell_1.6.1      data.table_1.12.8
#
#loaded via a namespace (and not attached):
# [1] SummarizedExperiment_1.14.1 tidyselect_0.2.5
# [3] purrr_0.3.2                 lattice_0.20-38
# [5] colorspace_1.4-1            vctrs_0.2.0
# [7] htmltools_0.3.6             stats4_3.6.1
# [9] blob_1.2.0                  XML_3.99-0.3
#[11] rlang_0.4.0                 R.oo_1.22.0
#[13] later_0.8.0                 pillar_1.4.2
#[15] withr_2.1.2                 glue_1.3.1
#[17] DBI_1.0.0                   R.utils_2.9.0
#[19] BiocParallel_1.18.1         BiocGenerics_0.30.0
#[21] bit64_0.9-7                 matrixStats_0.55.0
#[23] GenomeInfoDbData_1.2.1      zlibbioc_1.30.0
#[25] munsell_0.5.0               gtable_0.3.0
#[27] R.methodsS3_1.7.1           memoise_1.1.0
#[29] Biobase_2.44.0              IRanges_2.18.3
#[31] httpuv_1.5.2                GenomeInfoDb_1.20.0
#[33] AnnotationDbi_1.46.1        GSEABase_1.46.0
#[35] Rcpp_1.0.3                  xtable_1.8-4
#[37] promises_1.0.1              backports_1.1.5
#[39] scales_1.0.0                DelayedArray_0.10.0
#[41] S4Vectors_0.22.1            graph_1.62.0
#[43] annotate_1.62.0             XVector_0.24.0
#[45] mime_0.7                    bit_1.1-15.2
#[47] digest_0.6.21               dplyr_0.8.3
#[49] shiny_1.3.2                 GenomicRanges_1.36.1
#[51] grid_3.6.1                  bitops_1.0-6
#[53] tools_3.6.1                 magrittr_1.5
#[55] lazyeval_0.2.2              RCurl_1.95-4.12
#[57] tibble_2.1.3                RSQLite_2.1.2
#[59] crayon_1.3.4                pkgconfig_2.0.3
#[61] zeallot_0.1.0               Matrix_1.2-17
#[63] assertthat_0.2.1            R6_2.4.1
#[65] compiler_3.6.1
#}
