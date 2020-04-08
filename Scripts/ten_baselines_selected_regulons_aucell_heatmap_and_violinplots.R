# task: make selection of regulons (subset 10 baselines)
# calculate AUCell
# plot as heatmap (Fig3b)
# and violin plots (Fig3a & FigS4b)

#module load R/3.6.1-foss-2018a-X11-20180604
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#PATH=$PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/bin/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/lib/
#R
options(bitmapType="cairo")
options(stringsAsFactors=FALSE)
.libPaths("/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/3.6")
staging <- "/ddn1/vol1/staging/leuven/stg_00002/lcb"
#############################################################################################
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggforce)


wd <- file.path(staging, 
  "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL/")
setwd(wd)
zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices")

# function to scale each row from 0-1
scale01 <- function(x, low = min(x), high = max(x)) {
       x <- (x - low)/(high - low)
       x
}

# colours
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40kColours/MM40kColours.RData"))

# [1] make selection of regulons {
# AUCell ranking from make_AUCell_ranking_db.R
aucellRankings <- readRDS("aucellRankings.rds.gz")
				    
# get filtered regulons from 100 runs
regs <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_regulons.rds.gz")

# add selected regulons
# ZEB1 motif & track, RXRG motif, SOX6 track, BHLHE40 track, 
# TFAP2A motif, TFAP2B motif, ELF1 track, ETV4 track, USF2 track
selection <- c("ZEB1_mtf", "RXRG_mtf", "SOX6_trk", "BHLHE40_trk", "ZEB1_trk", 
	"TFAP2A_mtf", "TFAP2B_mtf", "ELF1_trk", "ETV4_trk", "USF2_trk")

allruns <- readRDS("all_runs_extended_summary_new_trk.rds.gz")
dbfilter <- list(
	mtf = with(allruns, tf_rec_100x_mtf > 0),
	trk = with(allruns, tf_rec_100x_trk > 0))

for (add in selection) {
	thistf <- unlist(strsplit(add, "_"))[1]
	thisdb <- unlist(strsplit(add, "_"))[2]
	regs[[add]] <- allruns[tf == thistf & dir == "pos" & dbfilter[[thisdb]], target]
}

# calculate AUC
regulonAUC <- AUCell_calcAUC(regs, aucellRankings,
	aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 20)
auc <- getAUC(regulonAUC) 

saveRDS(regulonAUC, "all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_regulonAUC01.rds.gz", compress = "gzip")
saveRDS(auc, "all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_aucell01.rds.gz", compress = "gzip")
saveRDS(regs, "all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_regulons.rds.gz", compress = "gzip")

#}

# [2] plot heatmap (Fig3b) {
auc <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_aucell01.rds.gz")
regulons <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_regulons.rds.gz")

# for plot: change regulon names
#e.g. SOX10_mtf -> "SOX10"; SOX10_trk -> "SOX10_track"
rownames(auc) <- gsub("_mtf", "", rownames(auc))
rownames(auc) <- gsub("_trk", "_track", rownames(auc))

names(regulons) <- gsub("_trk", "_track", names(regulons))
names(regulons) <- gsub("_mtf", "", names(regulons))

# order AUC matrix by cell lines
meta <- read.table(file.path(zdir, "12.ten_lines_BL/10.TextData/meta_data.tsv"))
meta <- meta[colnames(auc),]
lines <- c("MM001", "MM011", "MM031", "MM074", "MM087", "MM057", "A375", "MM029", "MM047", "MM099")
meta$Cell_Line <- factor(meta$Cell_Line, levels = lines)
meta <- meta[order(meta$Cell_Line),]
auc <- auc[, rownames(meta)]

# manually order regulons & remove motif regulons NR3C1 & BHLHE41 
regorder <- c("SOX10", "SOX10_track", "TFAP2A", "MITF", "MITF_track", 
	"IRF4", "USF2", "USF2_track", "HES6", "EGR3", 
	"SOX6", "SOX6_track", "NFATC2", "ELF1", "ELF1_track", 
	"ETV4", "ETV4_track", "JUN", "JUN_track", "SOX9", 
	"NR3C1_track", "IRF1", "IRF1_track", "FOSL2", "FOSL2_track",
	"ATF5", "NFIB", "FOSL1", "FOSL1_track", "IRF3", 
	"STAT1", "STAT1_track", "ZEB1", "ZEB1_track", "ATF4", 
	"ATF4_track", "SOX4", "MYC", "MYC_track", "TFAP2B", 
	"SOX11", "MYLK", "GATA3", "GATA3_track", "E2F1", 
	"E2F1_track", "MYBL2", "MYBL2_track", "RAD21")
auc <- auc[regorder,]
regulons <- regulons[regorder]

# order regulons
auc <- auc[plotten,]
regulons <- regulons[plotten]

# scale data
sauc <- t(apply(auc, 1, scale01))

# colour gradient for 95 percentile of the data
f2 <- colorRamp2(c(min(sauc), quantile(unlist(sauc), 0.95)), c("#EEEEEE", "red"))

# get regulon sizes
regulonsize <- data.frame(t(as.data.frame(lapply(regulons, length), row.names = c("num_targets"))))
rscol <- colorRamp2(c(min(regulonsize$num_targets),
	max(regulonsize$num_targets)), c("#EEEEEE", "purple"))

# get recurrence (from supplementary file GRNs_100runs_tenBL.tsv)
tfs <- unique(vapply(strsplit(regorder, "_"), '[', 1, FUN.VALUE=character(1)))
grn <- fread("GRNs_100runs_tenBL.tsv")
recurrence <- list()
for (tf in tfs) {
	recurrence[[tf]] <- max(grn[TF == tf, TF_recurrence_100runs_motif_db_tenBL])
	recurrence[[paste0(tf, "_track")]] <- max(grn[TF == tf, TF_recurrence_100runs_track_db_tenBL])
}
recurrence <- data.frame(t(as.data.frame(recurrence, row.names = c("recurrence"))))
recurrence <- recurrence[rownames(sauc),,drop = F]
reccol <- colorRamp2(c(min(recurrence$recurrence),
	max(recurrence$recurrence)), c("#EEEEEE", "slateblue3"))

# make row annotation from regulon sizes & recurrence
row_annot <- rowAnnotation(
	heatm1 = anno_simple(regulonsize$num_targets, width = unit(7, "mm"), col = rscol),
	text1 = anno_text(regulonsize$num_targets, gp = gpar(fontsize = 8), just = "right", 
		location = 0.4),
	heatm2 = anno_simple(recurrence$recurrence, width = unit(5.5, "mm"), col = reccol),
	text = anno_text(recurrence$recurrence, gp = gpar(fontsize = 8), just = "right",
		location = 0.53),
	show_annotation_name = FALSE, show_legend = TRUE,
	gap = unit(-3, "mm"))

# make regulon size & recurrence legends
lgd_list = list(
	Legend(col_fun = rscol, title = "size"),
	Legend(col_fun = reccol, title = "recur."))

# split cell state groups in columns
groups <- list(melanocytic = c("MM001", "MM031", "MM011"),
        intermediate = c("MM074", "MM087", "MM057"),
        A375 = c("A375"),
        mesenchymal_like = c("MM029", "MM047", "MM099"))

# plot cell state groups separately
hs <- list()
df <- as.data.frame(meta[, "Cell_Line", drop = F])
set <- names(groups[1])
tmpdf <- df[which(df$Cell_Line %in% groups[[set]]),, drop=F]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = MM40kColours["Cell_Line"], 
	show_annotation_name = FALSE)
hs[[set]] <- Heatmap(sauc[, rownames(tmpdf)],
		cluster_columns = FALSE, cluster_rows = FALSE, 
		col = f2, top_annotation = ha,
		show_heatmap_legend = F, show_column_names = FALSE, show_row_names = FALSE,
		left_annotation = row_annot)
set <- names(groups[2])
tmpdf <- df[which(df$Cell_Line %in% groups[[set]]),, drop=F]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = MM40kColours["Cell_Line"],
	show_annotation_name = FALSE)
hs[[set]] <- Heatmap(sauc[, rownames(tmpdf)],
		cluster_columns = FALSE, cluster_rows = FALSE, 
		col = f2, top_annotation = ha,
		show_heatmap_legend = F, show_column_names = FALSE, show_row_names = FALSE)
set <- names(groups[3])
tmpdf <- df[which(df$Cell_Line %in% groups[[set]]),, drop=F]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = MM40kColours["Cell_Line"],
	show_annotation_name = FALSE)
hs[[set]] <- Heatmap(sauc[, rownames(tmpdf)],
		cluster_columns = FALSE, cluster_rows = FALSE, 
		col = f2, top_annotation = ha,
		show_heatmap_legend = F, show_column_names = FALSE, show_row_names = FALSE)
set <- names(groups[4])
tmpdf <- df[which(df$Cell_Line %in% groups[[set]]),, drop=F]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = MM40kColours["Cell_Line"],
	show_annotation_name = FALSE)
hs[[set]] <- Heatmap(sauc[, rownames(tmpdf)],
		cluster_columns = FALSE, cluster_rows = FALSE, 
		col = f2, top_annotation = ha,
		show_heatmap_legend = T, name = "AUCell",
		show_column_names = FALSE, 
		row_names_gp = gpar(fontsize = 8))
hslist <- hs[[1]] + hs[[2]] +  hs[[3]] + hs[[4]]
draw(hslist, gap = unit(2, "mm"), annotation_legend_list = lgd_list)
#}

# [3] violin plots (Fig3a) {
aucs <- data.table(reshape::melt(auc))
setnames(aucs, c("regulon", "Cell_ID", "AUCell"))
meta$Cell_ID <- rownames(meta)
aucs <- merge(aucs, meta, by = "Cell_ID")[,.(Cell_ID, Cell_Line, Experiment, nGene, regulon, AUCell)]

selregs <- c("SOX10", "SOX6", "NFATC2", "EGR3", "ELF1", "ETV4")
aucs <- aucs[regulon %in% selregs]
aucs[, regulon := droplevels(regulon)]

p <- ggplot(aucs[regulon %in% selregs], aes(x = Cell_Line, y = AUCell, fill = Cell_Line)) +
	geom_violin(scale = "width") +
	scale_fill_manual(values = MM40kColours$Cell_Line) +

   	geom_jitter(aes(x = Cell_Line, y = AUCell),
		shape=16, size = 0.3, width = 0.2) +

	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf("test.pdf")
for (i in seq_along(selregs)) {
  print(p + facet_wrap_paginate(regulon~ ., nrow = 1, ncol = 1, page = i, scales = "free", drop = F))
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
#[1] grid      stats     graphics  grDevices utils     datasets  methods
#[8] base
#
#other attached packages:
#[1] ggforce_0.3.1        ggplot2_3.2.1        circlize_0.4.6
#[4] ComplexHeatmap_2.0.0 data.table_1.12.8
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3          pillar_1.4.2        compiler_3.6.1
# [4] RColorBrewer_1.1-2  plyr_1.8.4          digest_0.6.21
# [7] tibble_2.1.3        gtable_0.3.0        clue_0.3-57
#[10] pkgconfig_2.0.3     png_0.1-7           rlang_0.4.0
#[13] parallel_3.6.1      withr_2.1.2         dplyr_0.8.3
#[16] cluster_2.1.0       GlobalOptions_0.1.0 tidyselect_0.2.5
#[19] reshape_0.8.8       glue_1.3.1          R6_2.4.1
#[22] GetoptLong_0.1.7    polyclip_1.10-0     purrr_0.3.2
#[25] tweenr_1.0.1        farver_2.0.1        magrittr_1.5
#[28] scales_1.0.0        MASS_7.3-51.4       assertthat_0.2.1
#[31] shape_1.4.4         colorspace_1.4-1    labeling_0.3
#[34] lazyeval_0.2.2      munsell_0.5.0       crayon_1.3.4
#[37] rjson_0.2.20
#}
