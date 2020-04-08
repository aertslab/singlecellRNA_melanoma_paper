# task: summarize regulons from 100x scenic runs

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
library(hdf5r)

samples <- c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")
for (sample in samples) {
	regulons <- list()
	scenicdir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets")
	
	wd <- file.path(scenicdir, sample)
	setwd(wd)
	
	# [1] get regulons from 100 runs
	###################################### {
	mtfdir <- file.path(wd, "100_scenic_runs")
	trkdir <- file.path(wd, "rerunSCENICtracks100x")

	# 100x pyscenic runs
  	run <- "mlt"
  	dir <- "pos"

  	db <- "mtf"
  	for (i in 1:100) {
		reg <- readRDS(file.path(mtfdir, paste0("run_", i), 
			paste0("regulons_extended_", db, ".rds.gz")))
		regulons[[paste(run, i, db, dir, sep = "_")]] <- reg
  	}
	
  	db <- "trk"
  	for (i in 1:100) {
		reg <- readRDS(file.path(trkdir, paste0("run_", i), "reg_trk.rds.gz"))
		# take only positive regulons because we have only positive mtf regulons anyway
		regulons[[paste(run, i, db, dir, sep = "_")]] <- reg[[dir]]
  	}

	#remove empty lists
	regulons <- regulons[!lapply(regulons, length)==0] 

	#remove non-tfs
	nontfs <- read.table(file.path(staging, "kspan/resources/hg19/hg19_nonTFs.lst"), header = F)$V1
	regulons <- lapply(regulons, function(r) {
		r[!names(r) %in% nontfs]
	})
	#}
	saveRDS(regulons, file = "regulons_extended_mlt_plus_sgl_list_new_trk.rds.gz", compress = "gzip")

	# [2]  make one data frame with presence/absence
	################################################ {
	reg_long <- lapply(regulons, function(regset) {
		stack(regset, drop = TRUE)
	})
	reg_long <- rbindlist(reg_long, idcol = names(reg_long))
	colnames(reg_long) <- c("psrun", "target", "tf")
	reg_long[, regulon := paste0(psrun, "_", tf)]
	incMat <- as.data.frame(with(reg_long, table(target, regulon)) > 0L) + 0L
		dim(incMat)
	#}

	# [3] summarize all runs
	######################## {
	## 100x scenic runs
	
	#how often is tf-target connection found in 100x runs
	mlt_runs_mtf <- sapply(1:100, function(i) {paste0("mlt_", i, "_mtf_pos")})
		#how often is regulator found in 100x runs
		tf_recurrence <- unique(reg_long[psrun %in% mlt_runs_mtf, .(psrun, tf)])[, .N, by = tf]
		colnames(tf_recurrence) <- c("tf", "tf_rec_100x_mtf")

	#join both
	mlt_runs_mtf <- reg_long[psrun %in% mlt_runs_mtf, .N, by = .(target, tf)]
	colnames(mlt_runs_mtf) <- c("target", "tf", "conn_rec_100x_mtf")
	
	mlt_runs_mtf[, dir := "pos"]
	mlt_runs_mtf <- mlt_runs_mtf[tf_recurrence, on = "tf"]
	
	mlt_runs_trk <- sapply(1:100, function(i) {paste0("mlt_", i, "_trk_pos")})
		#how often is regulator found in 100x runs
		tf_recurrence <- unique(reg_long[psrun %in% mlt_runs_trk, .(psrun, tf)])[, .N, by = tf]
		colnames(tf_recurrence) <- c("tf", "tf_rec_100x_trk")

	#join both
	mlt_runs_trk <- reg_long[psrun %in% mlt_runs_trk, .N, by = .(target, tf)]
	colnames(mlt_runs_trk) <- c("target", "tf", "conn_rec_100x_trk")
	mlt_runs_trk[, dir := "pos"]
	mlt_runs_trk <- mlt_runs_trk[tf_recurrence, on = "tf"]
	
	#full outer join
	allruns <- merge(mlt_runs_trk, mlt_runs_mtf, by = c("target", "tf", "dir"), all = TRUE)
	#}
	saveRDS(allruns, file = "all_runs_extended_summary_new_trk.rds.gz", compress = "gzip")

	# [4] save multiruns track connections & regulons as text files
	################################################################ {
	fwrite(allruns[dir == "pos" & !is.na(conn_rec_100x_trk), 
		.(tf, target, tf_rec_100x_trk, conn_rec_100x_trk)][
			order(tf_rec_100x_trk, tf, conn_rec_100x_trk, decreasing = T)], 
		file = "hundred_runs_extended_summary_new_trk.tsv", 
		quote = F, col.names = TRUE, row.names = FALSE, sep = "\t")
	
	dir.create("new_trk_regulons")
	for (TF in unique(allruns[dir == "pos" & !is.na(conn_rec_100x_trk), tf])) {
		fwrite(allruns[tf == TF & dir == "pos" & 
			!is.na(conn_rec_100x_trk), 
			.(target, conn_rec_100x_trk)][order(conn_rec_100x_trk, decreasing = T)], 
		file = file.path("new_trk_regulons", paste0(TF, ".tsv")),
		quote = F, sep = "\t", col.names = F, row.names = F)
	}
	#}

	# [5] filter regulons & calculate AUCell
	########################################################################### {
	allruns <- allruns[, .(tf, target, dir,
		conn_rec_100x_trk, tf_rec_100x_trk, conn_rec_100x_mtf, tf_rec_100x_mtf)]
	for (i in names(allruns)) allruns[is.na(get(i)), (i):=0]
	nontfs <- read.table(file.path(staging, "kspan/resources/hg19/hg19_nonTFs.lst"), header = F)$V1

	#filter rule:
	#keep only positive TF-target connections
	#remove non-TF regulons (histone deacetylases, dna methylases, rna polymerases)
	#for 100x recurrent TFs, take targets that come up 80x
	#for at least 80x recurrent TFs, take all targets
	filt_trk <- with(allruns,
  		dir == "pos" &
  		!tf %in% nontfs &
  		((tf_rec_100x_trk == 100 & conn_rec_100x_trk >= 80) |
  		(tf_rec_100x_trk < 100 & tf_rec_100x_trk >= 80 & conn_rec_100x_trk > 0)))

	filt_mtf <- with(allruns,
  		dir == "pos" &
  		!tf %in% nontfs &
  		((tf_rec_100x_mtf == 100 & conn_rec_100x_mtf >= 80) |
  		(tf_rec_100x_mtf < 100 & tf_rec_100x_mtf >= 80 & conn_rec_100x_mtf > 0)))

	regs_trk <- sapply(as.vector(unique(allruns[filt_trk, tf])), function(TF) {
        	as.vector(allruns[filt_trk & tf == TF, target])
		}, simplify = FALSE, USE.NAMES = TRUE)
	names(regs_trk) <- paste0(names(regs_trk), "_trk")

	regs_mtf <- sapply(as.vector(unique(allruns[filt_mtf, tf])), function(TF) {
        	as.vector(allruns[filt_mtf & tf == TF, target])
		}, simplify = FALSE, USE.NAMES = TRUE)
	names(regs_mtf) <- paste0(names(regs_mtf), "_mtf")

	regs <- c(regs_mtf, regs_trk)

	summary(unlist(lapply(regs, length)))

	# rm regulons <10 targets
	regs <- regs[!unlist(lapply(regs, length)) < 10]

	# calculate AUCell
	#get gene expression matrix from loom file 
	loom <- H5File$new("log_CPM_woDoublets.loom", mode = "r") 
	dgem <- t(loom[["matrix"]]$read()) 
	genes <- loom[["row_attrs"]][["Gene"]]$read() 
	rownames(dgem) <- genes 
	cells <- loom[["col_attrs"]][["CellID"]]$read() 
	colnames(dgem) <- cells 
	loom$close_all() 
 
	#calculate AUCell 
	#use ranking from make_AUCell_ranking_db.R
	aucellRankings <- readRDS("aucellRankings.rds.gz")
	regulonAUC <- AUCell_calcAUC(regs, aucellRankings,  
        	aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 20) 
	auc <- getAUC(regulonAUC) 
 	
	saveRDS(regulonAUC, file = "all_pos_regulons_recurrent_100x_jw_new_trk_regulonAUC01.rds.gz", compress = "gzip")
	saveRDS(auc, file = "all_pos_regulons_recurrent_100x_jw_new_trk_aucell01.rds.gz", compress = "gzip") 
	saveRDS(regs, file = "all_pos_regulons_recurrent_100x_jw_new_trk_regulons.rds.gz", compress = "gzip") 
	
	#save as text file
	rownames(auc) <- gsub("_mtf", "", rownames(auc))
	rownames(auc) <- gsub("_trk", "_track", rownames(auc))
	auc <- as.data.table(auc, keep.rownames = "rn")
	fwrite(auc, file = "all_pos_regulons_recurrent_100x_jw_new_trk_aucell01.tsv", quote = F, sep = "\t")
	#}
}

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
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
#[1] doRNG_1.7.1       rngtools_1.4      pkgmaker_0.27     registry_0.5-1
#[5] foreach_1.4.7     hdf5r_1.2.0       AUCell_1.6.1      data.table_1.12.4
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2                  lattice_0.20-38
# [3] zeallot_0.1.0               digest_0.6.21
# [5] mime_0.7                    R6_2.4.0
# [7] GenomeInfoDb_1.20.0         backports_1.1.5
# [9] stats4_3.6.1                RSQLite_2.1.2
#[11] pillar_1.4.2                zlibbioc_1.30.0
#[13] rlang_0.4.0                 annotate_1.62.0
#[15] blob_1.2.0                  S4Vectors_0.22.1
#[17] R.utils_2.9.0               R.oo_1.22.0
#[19] Matrix_1.2-17               BiocParallel_1.18.1
#[21] stringr_1.4.0               RCurl_1.95-4.12
#[23] bit_1.1-14                  shiny_1.3.2
#[25] DelayedArray_0.10.0         compiler_3.6.1
#[27] httpuv_1.5.2                pkgconfig_2.0.3
#[29] BiocGenerics_0.30.0         htmltools_0.3.6
#[31] SummarizedExperiment_1.14.1 tibble_2.1.3
#[33] GenomeInfoDbData_1.2.1      IRanges_2.18.3
#[35] codetools_0.2-16            matrixStats_0.55.0
#[37] doMC_1.3.6                  XML_3.98-1.20
#[39] crayon_1.3.4                withr_2.1.2
#[41] later_0.8.0                 bitops_1.0-6
#[43] R.methodsS3_1.7.1           grid_3.6.1
#[45] xtable_1.8-4                GSEABase_1.46.0
#[47] DBI_1.0.0                   magrittr_1.5
#[49] bibtex_0.4.2                graph_1.62.0
#[51] stringi_1.4.3               XVector_0.24.0
#[53] promises_1.0.1              doParallel_1.0.15
#[55] vctrs_0.2.0                 iterators_1.0.12
#[57] tools_3.6.1                 bit64_0.9-7
#[59] Biobase_2.44.0              parallel_3.6.1
#[61] AnnotationDbi_1.46.1        GenomicRanges_1.36.1
#[63] memoise_1.1.0
#}
