# task: combine results of 100x scenic runs of 4 subsets to cytoscape networks
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
library(hdf5r)
library(AUCell)

zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices")
indir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/")

wd <- file.path(staging,
  "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/comb4subsets")
dir.create(wd, recursive = TRUE)
setwd(wd)

### functions {
# replace NAs with zero in [[data.table]]
set_na_to_zero <- function(DT) {
 	for (i in names(DT)) DT[is.na(get(i)), (i):=0]; return(DT)}

# connection file stats
stc <- function(dt) {
	print(paste0(nrow(dt), " connections, ", length(unique(dt$tf)), " tfs, ", 
		length(unique(dt$target)), " targets"))
}
read.gmt <- function (fileName) {
	tmp <- unname(sapply(readLines(fileName), function(x) strsplit(as.character(x), "\t")))
	tmp <- tmp[which(lengths(tmp) > 0)]
	names(tmp) <- sapply(tmp, function(x) x[1])
	lapply(tmp, function(x) x[3:length(x)])
}
#}

# [1] summarize results of four subsets {
### get data from 4 subsets
###########################
subsets <- c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")

allruns_list <- sapply(subsets, function(subs) { readRDS(file.path(indir, subs, "all_runs_extended_summary_new_trk.rds.gz")) }, simplify = FALSE, USE.NAMES = TRUE)
names(allruns_list) <- subsets

allruns_list[["all_connections"]] <- unique(rbindlist(
	lapply(allruns_list, function(run) {run[, .(tf, target, dir)]})))
colnames(allruns_list[["12.ten_lines_BL"]])[4:9] <- paste0(colnames(allruns_list[["12.ten_lines_BL"]])[4:9], "_tenBL")
colnames(allruns_list[["5.MM074_SOX"]])[4:9] <- paste0(colnames(allruns_list[["5.MM074_SOX"]])[4:9], "_mm074")
colnames(allruns_list[["9.MM057_SOX_wo_TL"]])[4:9] <- paste0(colnames(allruns_list[["9.MM057_SOX_wo_TL"]])[4:9], "_mm057")
colnames(allruns_list[["10.M0087_SOX_wo_TL"]])[4:9] <- paste0(colnames(allruns_list[["10.M0087_SOX_wo_TL"]])[4:9], "_mm087")


### join
########
allruns <- Reduce(function(...) merge(..., on = c("target", "tf", "dir"), all = TRUE), allruns_list)
dim(allruns)
	#[1] 1090057      27
#replace all na's with zero
for (i in names(allruns)) {allruns[is.na(get(i)), (i):=0]}

#remove non-tfs
nontfs <- read.table(file.path(staging, "kspan/resources/hg19/hg19_nonTFs.lst"), header = F)$V1
allruns <- allruns[!tf %in% nontfs]
dim(allruns)
	#[1] 1090057      27
	
### sum up connection recurrences
#################################
#filter rule:
#for 100x recurrent TFs, take targets that come up 80x
#for at least 80x recurrent TFs, take all targets

allruns[, mlt_tenBL := ifelse(
	(tf_rec_100x_trk_tenBL >= 80 & conn_rec_100x_trk_tenBL >= tf_rec_100x_trk_tenBL*0.8) |
	(tf_rec_100x_mtf_tenBL >= 80 & conn_rec_100x_mtf_tenBL >= tf_rec_100x_mtf_tenBL*0.8)
  	, 1, 0), by = .(tf, target, dir)]
allruns[, mlt_mm074 := ifelse(
	(tf_rec_100x_trk_mm074 >= 80 & conn_rec_100x_trk_mm074 >= tf_rec_100x_trk_mm074*0.8) |
	(tf_rec_100x_mtf_mm074 >= 80 & conn_rec_100x_mtf_mm074 >= tf_rec_100x_mtf_mm074*0.8)
  	, 1, 0), by = .(tf, target, dir)]
allruns[, mlt_mm057 := ifelse(
	(tf_rec_100x_trk_mm057 >= 80 & conn_rec_100x_trk_mm057 >= tf_rec_100x_trk_mm057*0.8) |
	(tf_rec_100x_mtf_mm057 >= 80 & conn_rec_100x_mtf_mm057 >= tf_rec_100x_mtf_mm057*0.8)
  	, 1, 0), by = .(tf, target, dir)]
allruns[, mlt_mm087 := ifelse(
	(tf_rec_100x_trk_mm087 >= 80 & conn_rec_100x_trk_mm087 >= tf_rec_100x_trk_mm087*0.8) |
	(tf_rec_100x_mtf_mm087 >= 80 & conn_rec_100x_mtf_mm087 >= tf_rec_100x_mtf_mm087*0.8)
  	, 1, 0), by = .(tf, target, dir)]

allruns[, conn_rec_100x := sum(
	mlt_tenBL, mlt_mm074, mlt_mm057, mlt_mm087), by = .(tf, target, dir)]

#include biopsy data (Tirosh et al., 2016: MM-GSE72056_pySCENIC.txt)
####################################################################
biop <- readRDS(file.path(staging, "kspan/analyses/PublicData/Tirosh2016/all_extended_regulons_new.rds.gz"))
biop[, biop := .N, by = c("tf", "target", "dir")]	#sum connection recurrence
biop <- unique(biop[, .(tf, target, dir, biop)])
biop <- biop[!tf %in% nontfs]
allruns <- merge(allruns, biop, on = c("target", "tf", "dir"), all = TRUE)
allruns <- set_na_to_zero(allruns)

#include ATAC peaks & H3K27ac data
##################################
atacgenes <- fread("../comb4subsets/atac_genes.lst", header = F)$V1
allruns[, atac := ifelse(target %in% atacgenes, 1, 0)]
allruns[, .N, by = atac]
	#   atac       N
	#1:    1 1193378
	#2:    0   56298
	# 
atac_ac_genes <- fread("../comb4subsets/atac_ac_genes.lst", header = F)$V1
allruns[, atac_ac := ifelse(target %in% atac_ac_genes, 1, 0)]
allruns[, .N, by = atac_ac]
	#   atac_ac      N
	#1:       0 780789
	#2:       1 468887
#}
saveRDS(allruns, file = "allruns.rds.gz", compress = "gzip")
	#a99b12aefbd040ca2025bd13360ddc41d547521b  allruns.rds.gz
	
# [2] filter by link weight {
# get grnboost2 results of 100 runs & combine into one table (per subset)
for (sample in c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")) {
	grnbpath <- file.path(staging, "dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices", sample, "pySCENIC_100x/out/grnboost2withoutDask")
	grnrunlist <- list(
		c(1:10), c(11:20), c(21:30), c(31:40), c(41:50), 
		c(51:60), c(61:70), c(71:80), c(81:90), c(91:100))
	grnruns <- list()

	# read grnboost results in chunks of 10
	for (subs in 1:length(grnrunlist)) {
		print(paste0("subset ", subs))
	
		lls <- list()
		for (run in grnrunlist[[subs]]) {
			print(paste0("run ", run))
	
			file <- file.path(grnbpath, paste0("run_", run, "/run_", run, "__adj.tsv"))

			ll <- fread(cmd = paste0("awk '$3>=1' ", file), 
				col.names = c("TF", "Target", paste0("run", run)))

			ll[, tt := paste0(TF, "#", Target)]
			setkey(ll, tt)
			ll[, `:=`(TF = NULL, Target = NULL)]
			lls[[run]] <- ll
		}
		# merge 
		llt <- Reduce(function(...) merge(..., all.x = TRUE), lls[grnrunlist[[subs]]])
		# save
		grnruns[[subs]] <- llt
		# clean up
		rm(lls, llt)
	}
	
	# merge
	gbmerged <- Reduce(function(...) merge(..., all.x = TRUE), grnruns)

	# set NAs to 0
	set_na_to_zero(gbmerged)

	# calculate means of weights
	gbmerged[, meanweight := rowMeans(.SD), by = tt, .SDcols = c(2:ncol(gbmerged))]
	dim(gbmerged)
	#[1] 400661    101

	# bring into link-list format
	gbmerged <- gbmerged[, .(tt, meanweight)]
	gbmerged[, TF := vapply(strsplit(tt, "#"), '[', 1, FUN.VALUE=character(1))]
	gbmerged[, Target := vapply(strsplit(tt, "#"), '[', 2, FUN.VALUE=character(1))]
	gbmerged <- gbmerged[, .(TF, Target, meanweight)]
	colnames(gbmerged) <- c("TF", "target", "importance")
	setDF(gbmerged)
	
	# sort
	gbmerged <- gbmerged[order(gbmerged$importance, decreasing = T),]

	# save
	saveRDS(gbmerged, file = "100xGRNboost2_LinkList_sorted.rds.gz", compress = "gzip")
	write.table(gbmerged, file = "100xGRNboost2_LinkList_sorted.tsv", row.names = F, sep = "\t", quote = F)
}

# combine results of 4 subsets
gb <- list()
for (sample in c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")) {
	gb[[sample]] <- readRDS(file.path(indir, sample, "100xGRNboost2_LinkList_sorted.rds.gz"))
	#remove any scores below 90th percentile
	gb[[sample]] <- gb[[sample]][gb[[sample]]$importance >= quantile(gb[[sample]]$importance, prob = 0.9),]
	gb[[sample]]$subset <- strsplit(sample, "\\.")[[1]][2]
}
gb <- rbindlist(gb)
colnames(gb)[1] <- "tf"
setkeyv(gb, c("tf", "target"))
gb <- unique(gb[, .(tf, target)])
	#165177 connections

# connections with high weight in any grnboost2 run
allrunsf <- merge(allruns, gb, by = c("tf", "target"), all.x = FALSE)
	stc(allruns)
	#[1] "1248751 connections, 1511 tfs, 18842 targets" 
	
	stc(allrunsf)
	#[1] "33159 connections, 548 tfs, 8122 targets"
#}
saveRDS(allrunsf, file = "allruns_filtered.rds.gz", compress = "gzip")
	#48caf044a3a30dd750d0d50555159473ab8bac79  allruns_filtered.rds.gz
	
# [3] filter by connection recurrence in 100 runs, biopsies and atac (+ac) peaks -> Fig 4 {
prefix <- "ext_pos_atac_and_ac_and_conn_rec_mlt_ge4_or_mlt_ge3_and_biopsies_ge1_weigth_ge_90thperc"
filt <- with(allrunsf,
		dir == "pos" &
		atac_ac > 0 &
		(conn_rec_100x >= 4 |
		(conn_rec_100x >= 3 & biop >= 1)) )

stc(allrunsf[filt])
	#[1] "1345 connections, 55 tfs, 768 targets"
	
subs <- allrunsf[filt, .(tf, target, dir)]

	## save files for cytoscape
	write.table(subs[,.(tf, dir, target)],
		file = paste0(prefix, ".sif"),
		sep = "\t", row.names = F, col.names = F, quote = F)

	#mrna file
	mrna <- fread(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL/MM40k100xSCENIC_tenBL_stateSpecific_networks/MMgroup_woA375__aggr_norm_scaled.mrna"))
	mrna <- mrna[gene %in% unique(c(subs[, tf], subs[, target]))]
	mrna[, alltfs := ifelse(gene %in% unique(allruns[, tf]), TRUE, FALSE)]
	mrna[, tfs_this_network := ifelse(gene %in% unique(subs[, tf]), TRUE, FALSE)]

	write.table(mrna,
		file = paste0(prefix, ".mrna"),
		sep = "\t", row.names = F, col.names = T, quote = F)
		
	#connection strength
	mrna <- allrunsf[filt]
	mrna[, mlt_tenBL_mtf := ifelse(
		tf_rec_100x_mtf_tenBL >= 80 & conn_rec_100x_mtf_tenBL >= tf_rec_100x_mtf_tenBL*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_tenBL_trk := ifelse(
		tf_rec_100x_trk_tenBL >= 80 & conn_rec_100x_trk_tenBL >= tf_rec_100x_trk_tenBL*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_tenBL_sum := sum(mlt_tenBL_mtf, mlt_tenBL_trk), by = 1:NROW(mrna)]
		
	mrna[, mlt_mm074_mtf := ifelse(
		tf_rec_100x_mtf_mm074 >= 80 & conn_rec_100x_mtf_mm074 >= tf_rec_100x_mtf_mm074*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm074_trk := ifelse(
		tf_rec_100x_trk_mm074 >= 80 & conn_rec_100x_trk_mm074 >= tf_rec_100x_trk_mm074*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm074_sum := sum(mlt_mm074_mtf, mlt_mm074_trk), by = 1:NROW(mrna)]
	
	mrna[, mlt_mm057_mtf := ifelse(
		tf_rec_100x_mtf_mm057 >= 80 & conn_rec_100x_mtf_mm057 >= tf_rec_100x_mtf_mm057*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm057_trk := ifelse(
		tf_rec_100x_trk_mm057 >= 80 & conn_rec_100x_trk_mm057 >= tf_rec_100x_trk_mm057*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm057_sum := sum(mlt_mm057_mtf, mlt_mm057_trk), by = 1:NROW(mrna)]
	
	mrna[, mlt_mm087_mtf := ifelse(
		tf_rec_100x_mtf_mm087 >= 80 & conn_rec_100x_mtf_mm087 >= tf_rec_100x_mtf_mm087*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm087_trk := ifelse(
		tf_rec_100x_trk_mm087 >= 80 & conn_rec_100x_trk_mm087 >= tf_rec_100x_trk_mm087*0.8
  		, 1, 0), by = .(tf, target, dir)]
	mrna[, mlt_mm087_sum := sum(mlt_mm087_mtf, mlt_mm087_trk), by = 1:NROW(mrna)]

	
	# how often is connection found
	mrna[, mlt_sum := sum(mlt_tenBL_sum, mlt_mm074_sum, mlt_mm057_sum, mlt_mm087_sum, biop),
		by = 1:NROW(mrna)]

	mrna[, .N, by = mlt_sum][order(mlt_sum)]
   		#mlt_sum   N
   		#1:       4 572
   		#2:       5 254
   		#3:       6  99
   		#4:       7 124
   		#5:       8 146
   		#6:       9  54
   		#7:      10  96
		
	# is connection found it both mtf & trk
	mrna[, mlt_conf := ifelse(
		mlt_tenBL_sum > 1 | mlt_mm074_sum > 1 | mlt_mm057_sum > 1 | mlt_mm087_sum > 1 | biop > 1, 1, 0), by = 1:NROW(mrna)]
	mrna[, .N, by = mlt_conf]
	   #mlt_conf   N
	   #1:        1 578	#found by motif & track in at least one run
	   #2:        0 767	#found by motif or track
	   
	mrna[, shared_name := paste0(tf, " (pos) ", target)]
	write.table(mrna[, .(shared_name, mlt_sum, mlt_conf)],
		file = paste0(prefix, "_connection_strength.mrna"),
		sep = "\t", row.names = F, col.names = T, quote = F)

# mlt_sum:
# sum of connection recurrence in five data sets (four subsets this paper + Tirosh biopsy data)
# in each data set, a connection can be found once (in either the motif or the track databases)
# or twice (when found in both databases)
# => minimum connection recurrence is four (found once in all four subsets or in three subsets plus biopsies), maximum is ten (found in all five data sets with both motif and track databases)

# mlt_conf:
# is set to 0 when the connection is found using only one of the databases
# set to 1 when the connection is found in both motif and track databases in any of the five runs
#}

# [4] make supplementary networks Fig S13 + S14 from subset 10 baselines {
# [A] get filtered regulons for subset 10 baselines from filter_regulons.R
allruns <- readRDS(file.path(indir, "12.ten_lines_BL/all_runs_extended_summary_new_trk.rds.gz"))
filt <- with(allruns,
	dir == "pos" &
	!tf %in% nontfs &
	(((tf_rec_100x_trk == 100 & conn_rec_100x_trk >= 80) |
	  (tf_rec_100x_trk < 100 & tf_rec_100x_trk >= 80 & conn_rec_100x_trk > 0)) |
	 ((tf_rec_100x_mtf == 100 & conn_rec_100x_mtf >= 80) |
	  (tf_rec_100x_mtf < 100 & tf_rec_100x_mtf >= 80 & conn_rec_100x_mtf > 0))))
allruns <- allruns[filt]

#include biopsy data
biop <- readRDS(file.path(staging, "kspan/analyses/PublicData/Tirosh2016/all_extended_regulons_new.rds.gz"))
biop <- biop[!tf %in% nontfs & dir == "pos"]
biop[, biop_mtf := ifelse(db == "mtf", 1, 0)]
biop[, biop_trk := ifelse(db == "trk", 1, 0)]
biop[, db := NULL]
cols <- c("biop_mtf", "biop_trk")
biop <- unique(biop[, (cols) := lapply(.SD, function(col) sum(col)), .SDcols = cols, by = .(tf, target)])
allruns <- merge(allruns[, .(target, tf, dir, conn_rec_100x_trk, conn_rec_100x_mtf)],
		         biop, on = c("target", "tf", "dir"), all = TRUE)
allruns <- set_na_to_zero(allruns)

## intermediate state network {
tfs <- c("NFATC2", "EGR3", "ETV4", "ELF1", "SOX6")
stc(allruns[tf %in% tfs & dir == "pos" & (conn_rec_100x_trk > 0 | conn_rec_100x_mtf > 0)])
	#[1] "1711 connections, 5 tfs, 1427 targets" 
stc(allruns[tf %in% tfs & dir == "pos" & (conn_rec_100x_mtf == 100 | conn_rec_100x_trk == 100)])
        #[1] "707 connections, 5 tfs, 650 targets"

#for targets regulated by at least 2 TFs: keep all
#for targets regulated by only one TF: keep only those that have 100x recurrence
tmpnw <- allruns[tf %in% tfs & (conn_rec_100x_trk >= 80 | conn_rec_100x_mtf >= 80)]
stc(tmpnw[target %in% tmpnw[, .N, by = target][N > 1, target] |
                    (conn_rec_100x_mtf == 100 | conn_rec_100x_trk == 100)])
	#[1] "1010 connections, 5 tfs, 726 targets"
write.table(tmpnw[target %in% tmpnw[, .N, by = target][N > 1, target] |
	                            (conn_rec_100x_mtf == 100 | conn_rec_100x_trk == 100),
	file = "int_network_conn_rec100_orShared.sif",
        sep = "\t", row.names = F, col.names = F, quote = F)

# edge annotation: connection recurrence & confirmation with biopsy regulons
tmpnw[, shared_name := paste0(tf, " (activates) ", target)]
tmpnw[, max_conn_strength := max(conn_rec_100x_mtf, conn_rec_100x_trk), by = 1:NROW(tmpnw)]
tmpnw[, conn_strength := ifelse(max_conn_strength< 90, 1, 2)]
tmpnw[max_conn_strength == 100, conn_strength := 3]
tmpnw[, .N, by = conn_strength]
        #   conn_strength   N
        #1:             3 707
        #2:             2 680
        #3:             1 324

#note: only low-confidence track recurrence, max 71x
tmpnw[, motif := ifelse(conn_rec_100x_mtf >= 80, 1, 0)]
tmpnw[, track := ifelse(conn_rec_100x_trk >= 10, 1, 0)] #234 of 1477 motif connections confirmed

stc(tmpnw[conn_rec_100x_trk > 0])
        #[1] "528 connections, 3 tfs, 516 targets"
stc(tmpnw[conn_rec_100x_trk >= 10])
        #[1] "234 connections, 3 tfs, 231 targets"

write.table(tmpnw[, .(shared_name, conn_rec_100x_mtf, conn_rec_100x_trk, conn_strength,
	                      biop_mtf, biop_trk, motif, track)],
	file = "int_network_conn_rec_trk.attr",
	sep = "\t", row.names = F, quote = F)
#}

## SOX10 network {
#keep only connections found in our single-cell data (not the additional ones found in the biopsies)
sox10 <- allruns[tf == "SOX10" & (conn_rec_100x_trk > 0 | conn_rec_100x_mtf > 0)]
stc(sox10)
#[1] "594 connections, 1 tfs, 594 targets"
#76 found in track AND motif (of those, 30 confirmed by biopsies (30 mtf & 13 of those track, too))
#115 of all targets  confirmed by biopsies
#3 targets found by track only

sox10[, dir := "activates"]
write.table(sox10[, .(tf, dir, target)],
	file = "sox10_network_conn_rec_full.sif",
	sep = "\t", row.names = F, col.names = F, quote = F)

# edge annotation: connection recurrence & confirmation with biopsy motif regulon
sox10[, shared_name := paste0(tf, " (activates) ", target)]
sox10[, max_conn_strength := max(conn_rec_100x_mtf, conn_rec_100x_trk), by = 1:NROW(sox10)]
sox10[, conn_strength := ifelse(max_conn_strength< 90, 1, 2)]
sox10[max_conn_strength == 100, conn_strength := 3]

sox10[, motif := ifelse(conn_rec_100x_mtf >= 80, 1, 0)]
sox10[, track := ifelse(conn_rec_100x_trk >= 80, 1, 0)]

write.table(sox10[, .(shared_name, conn_rec_100x_mtf, conn_rec_100x_trk, conn_strength, 
		                      biop_mtf, biop_trk, motif, track)],
	file = "sox10/sox10_network_conn_rec_full.attr",
	sep = "\t", row.names = F, quote = F)
#}
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
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
#[1] AUCell_1.6.1      hdf5r_1.2.0       data.table_1.12.4
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2                  XVector_0.24.0
# [3] GenomeInfoDb_1.20.0         pillar_1.4.2
# [5] compiler_3.6.1              later_0.8.0
# [7] zlibbioc_1.30.0             tools_3.6.1
# [9] R.methodsS3_1.7.1           R.utils_2.9.0
#[11] bitops_1.0-6                zeallot_0.1.0
#[13] digest_0.6.21               bit_1.1-14
#[15] lattice_0.20-38             annotate_1.62.0
#[17] RSQLite_2.1.2               memoise_1.1.0
#[19] tibble_2.1.3                pkgconfig_2.0.3
#[21] rlang_0.4.0                 Matrix_1.2-17
#[23] graph_1.62.0                DelayedArray_0.10.0
#[25] shiny_1.3.2                 DBI_1.0.0
#[27] parallel_3.6.1              GenomeInfoDbData_1.2.1
#[29] vctrs_0.2.0                 S4Vectors_0.22.1
#[31] IRanges_2.18.3              grid_3.6.1
#[33] stats4_3.6.1                bit64_0.9-7
#[35] GSEABase_1.46.0             Biobase_2.44.0
#[37] R6_2.4.0                    AnnotationDbi_1.46.1
#[39] BiocParallel_1.18.1         XML_3.98-1.20
#[41] magrittr_1.5                blob_1.2.0
#[43] matrixStats_0.55.0          GenomicRanges_1.36.1
#[45] backports_1.1.5             promises_1.0.1
#[47] htmltools_0.3.6             BiocGenerics_0.30.0
#[49] SummarizedExperiment_1.14.1 mime_0.7
#[51] xtable_1.8-4                httpuv_1.5.2
#[53] RCurl_1.95-4.12             crayon_1.3.4
#[55] R.oo_1.22.0
#}
