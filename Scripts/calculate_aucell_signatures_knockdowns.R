# task: calculate AUCell on literature signatures on three lines SOX10-KD (10x & Drop-Seq)
# make new violin plots (Fig6h) and heatmaps (Fig S17d)

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
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggforce)
library(hdf5r)

# functions
read.gmt <- function (fileName) {
	tmp <- unname(sapply(readLines(fileName), function(x) strsplit(as.character(x), "\t")))
	tmp <- tmp[which(lengths(tmp) > 0)]
	names(tmp) <- sapply(tmp, function(x) x[1])
	lapply(tmp, function(x) x[3:length(x)])
}
# scale each row from 0-1
scale01 <- function(x, low = min(x), high = max(x)) {
       x <- (x - low)/(high - low)
       x
}
# colours 
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40kColours/MM40kColours.RData"))

wd <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/40.three_MM_lines_SOX_wo_TL/20.Plots/MM40k3linesAUCsignVioplots")
dir.create(wd)
setwd(wd)

samples <- c("5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL", "53.MM057_DropSeq")
lines <- c("MM074", "MM087", "MM057")

# [1] for dropseq data, first calculate aucell ranking
loomfile <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/53.MM057_DropSeq/30.SCENIC/pySCENIC/logCPM.loom")
loom <- H5File$new(loomfile, mode = "r") 

dgem <- t(loom[["matrix"]]$read())
genes <- loom[["row_attrs"]][["Gene"]]$read()
rownames(dgem) <- genes
cells <- loom[["col_attrs"]][["CellID"]]$read()
colnames(dgem) <- cells

set.seed(123) 
aucellRankings <- AUCell_buildRankings(dgem, nCores = 20, plotStats = F) 
saveRDS(aucellRankings, file = file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/53.MM057_DropSeq/30.SCENIC/pySCENIC/aucellRankings.rds.gz"), compress = "gzip") 
aucellRankings@nGenesDetected
    #min      1%      5%     10%     50%    100% 
 #874.00  955.77 1053.00 1178.70 2063.50 6024.00 

rankings <- list()
rankings[["53.MM057_DropSeq"]] <- aucellRankings

# [2] get AUCell rankings of 10x data
#from [[MM40k100xSCENIC_getRegulons_switchers.R]]
for (sample in samples[1:3]) {
	rankfile <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets", sample, "aucellRankings.rds.gz")
	aucellRankings <- readRDS(rankfile)
	rankings[[sample]] <- aucellRankings
}

# [3] calculate AUCell
# new signatures [[GeneSets]]
signatures <- read.gmt(file.path(staging, "kspan/resources/genesets/melanoma_signatures_190618.gmt"))

aucells <- list()
for (sample in samples) {
	auc <- AUCell_calcAUC(signatures, rankings[[sample]], 
		aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 20)
	saveRDS(auc, file = file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets", sample, "auc_signatures.rds.gz"))
	saveRDS(auc, file = paste0("auc_signatures_", sample, ".rds.gz"), compress = "gzip")
	aucells[[sample]] <- auc
}

# [4] meta data
# 10x
meta <- read.table(file.path(staging, 
	"zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/40.three_MM_lines_SOX_wo_TL/10.TextData/meta_data.tsv"))
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/10.analysis_with_all_data/00.RData/demuxlet_doublets.RData"), verbose = TRUE)
meta <- meta[!rownames(meta) %in% demuxlet_doublets,]
setDT(meta, keep.rownames = "Cell_ID")
meta[, technique := "tenX"]

# dropseq metadata
dm <- read.table(file.path(staging,
	"zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/53.MM057_DropSeq/10.TextData/meta_data.tsv"))
setDT(dm, keep.rownames = "Cell_ID")
dm[Experiment == "SOX10_KD", Experiment := "SOX10_72h"]
setDT(dm, keep.rownames = "Cell_ID")
dm[, technique := "dropseq"]
dm[, CE_ID := paste0(CE_ID, "_dropseq")]
	
meta <- rbind(meta[, .(Cell_ID, Cell_Line, Experiment, CE_ID, technique)],
	dm[, .(Cell_ID, Cell_Line, Experiment, CE_ID, technique)])

selsigs <- setNames(c(
	"VERF_INV", "VERF_PRO", "HOEK_INV", "HOEK_PRO",
	"raj_resistance", "WU_CELL_MIGRATION", "hugo_IPRES", "AXL_SIGNATURE_GODING",
	"GO_PIGMENTATION", "ANASTASSIOU_CANCER_MESENCHYMAL_TRANSITION_SIGNATURE",
	"VECCHI_GASTRIC_CANCER_ADVANCED_VS_EARLY_UP", "LU_TUMOR_VASCULATURE_UP",
	
	"HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
	"HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_ANGIOGENESIS",
	
	"KEGG_MELANOGENESIS", "KEGG_DNA_REPLICATION", "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
	"KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION", "KEGG_ALLOGRAFT_REJECTION",
	"KEGG_GRAFT_VERSUS_HOST_DISEASE", "RIESENBERG_MELANOMA_TNF_RESPONSE_GODING"),
	
	c("Verfaillie invasive", "Verfaillie proliferative", "Hoek invasive", "Hoek proliferative",
	"Raj resistance", "Wu cell migration", "Hugo IPRES", "Tirosh AXL signature",
	"GO pigmentation", "Cheng Slug-based EMT",
	"Vecchi gastric cancer advanced vs early up", "Lu tumor vasculature up",
	
	"Hallmark E2F targets", "Hallmark G2M checkpoint", "Hallmark EMT",
	"Hallmark TGFbeta signaling", "Hallmark angiogenesis",
	
	"KEGG melanogenesis", "KEGG DNA replication", "KEGG regulation of actin cytoskeleton",
	"KEGG leukocyte transendothelial migration", "KEGG allograft rejection",
	"KEGG graft vs host disease", "Riesenberg Melanoma TNF response"))

auc <- rbindlist(lapply(aucells, function(auc) {
	auc <- SummarizedExperiment::assays(auc)[["AUC"]]
	auc <- data.table(reshape::melt(auc[as.vector(selsigs),]))
	setnames(auc, c("signature", "Cell_ID", "AUCell"))
	auc
}))

# make long table
auc <- merge(auc, meta, by = "Cell_ID")
saveRDS(auc, file = "auc_signatures_combined.rds.gz", compress = "gzip")

###### heatmap
############## {
## note: dropseq has higher aucell values than 10x (much fewer genes -> different rank threshold)
## scaling from 0-1 across tenX and dropseq results in dark red dropseq, tenX colours faint
## solution: scale separately!

ausels <- lapply(aucells, function(a) {
	getAUC(a[as.vector(selsigs),])
})
ausels <- do.call(cbind, ausels)

df <- data.frame(meta[, .(Cell_Line, Experiment, technique)], row.names = meta[, Cell_ID])
plotcolours <- MM40kColours[c("Cell_Line", "Experiment")]
ps_tenX <- t(apply(ausels[, meta[technique == "tenX", Cell_ID]], 1, scale01))
#dropseq data: rm baseline!
ps_drops <- t(apply(ausels[, meta[Experiment != "BL" & technique == "dropseq", Cell_ID]], 1, scale01))

#cluster tenx data
dend <- as.dendrogram(hclust(dist(ps_tenX)))

### plot tenX data with white space between lines
hs <- list()
line <- lines[1]
#sort by experiment
tmpdf <- df[df$Cell_Line==line & df$technique == "tenX", c("Cell_Line", "Experiment")]
tmpdf$Experiment <- factor(tmpdf$Experiment, levels = names(MM40kColours$Experiment))
tmpdf <- tmpdf[with(tmpdf, order(Experiment)),]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = plotcolours, show_annotation_name = FALSE)
hs[[line]] <- Heatmap(ps_tenX[,rownames(tmpdf)], 
	cluster_rows = dend, cluster_columns = FALSE,
	col = f1, show_heatmap_legend = F,
	show_column_names = FALSE, show_row_names = FALSE, top_annotation = ha)

line <- lines[2]
#sort by experiment
tmpdf <- df[df$Cell_Line==line & df$technique == "tenX", c("Cell_Line", "Experiment")]
tmpdf$Experiment <- factor(tmpdf$Experiment, levels = names(MM40kColours$Experiment))
tmpdf <- tmpdf[with(tmpdf, order(Experiment)),]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = plotcolours, show_annotation_name = FALSE)
hs[[line]] <- Heatmap(ps_tenX[,rownames(tmpdf)], 
	cluster_rows = FALSE, cluster_columns = FALSE,
	col = f1, show_heatmap_legend = F,
	show_column_names = FALSE, show_row_names = FALSE, top_annotation = ha)
	
line <- lines[3]
#sort by experiment
tmpdf <- df[df$Cell_Line==line & df$technique == "tenX", c("Cell_Line", "Experiment")]
tmpdf$Experiment <- factor(tmpdf$Experiment, levels = names(MM40kColours$Experiment))
tmpdf <- tmpdf[with(tmpdf, order(Experiment)),]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = plotcolours, show_annotation_name = FALSE)
hs[[line]] <- Heatmap(ps_tenX[,rownames(tmpdf)], 
	cluster_rows = FALSE, cluster_columns = FALSE,
	col = f1, show_heatmap_legend = T, name = "AUCell\nscaled",
	show_column_names = FALSE, row_names_gp = gpar(fontsize = 6), top_annotation = ha)
hslist <- hs[[1]] + hs[[2]] + hs[[3]]

pdf("auc_selectedSigs_heatmaps_tenX.pdf")
draw(hslist, gap = unit(0.5, "mm"))
dev.off()

### plot dropseq data
line <- "MM057"
#sort by experiment
tmpdf <- df[df$Experiment != "BL" & 
	df$Cell_Line==line & 
	df$technique == "dropseq", c("Cell_Line", "Experiment")]
tmpdf$Experiment <- factor(tmpdf$Experiment, levels = names(MM40kColours$Experiment))
tmpdf <- tmpdf[with(tmpdf, order(Experiment)),]
ha <- HeatmapAnnotation(df = tmpdf, show_legend = F, col = plotcolours, show_annotation_name = FALSE)
hs[["dropseq"]] <- Heatmap(ps_drops[,rownames(tmpdf)], 
	cluster_rows = dend, cluster_columns = FALSE,
	col = f1, show_heatmap_legend = T, name = "AUCell\nscaled",
	show_column_names = FALSE, row_names_gp = gpar(fontsize = 6), top_annotation = ha)
pdf("auc_selectedSigs_heatmaps_dropseq.pdf", width = 4)
png("auc_selectedSigs_heatmaps_dropseq.png", width = 400, heigh = 800)
draw(hs[["dropseq"]])
dev.off()
#}

###### violin plots
################### {
auc[technique == "dropseq", Cell_Line := "MM057_DropSeq"]
auc$Cell_Line <- factor(auc[, Cell_Line], levels =  c("MM074", "MM087", "MM057", "MM057_DropSeq"))

#set colours of non-existing dropseq treatments to white
colours <- MM40kColours[["Experiment"]][unique(auc[,Experiment])]
col2 <- setNames(rep(as.vector(colours), length(unique(auc$Cell_Line))),
	paste(rep(names(colours), length(unique(auc$Cell_Line))), unique(auc$Cell_Line), sep = "_"))
col2[c("BL_MM057_DropSeq", "SOX10_24h_MM057_DropSeq", "SOX10_48h_MM057_DropSeq")] <- "white"
auc[,fillcol := paste(Experiment, Cell_Line, sep = "_")]

p <- ggplot(auc,
	aes(x = Cell_Line, y = AUCell, colour = fillcol , fill = fillcol)) +
	geom_violin() +
	scale_colour_manual(values = col2, guide = F) +
	scale_fill_manual(values = col2, guide = F) +
	
	#change x-axis tick labels
	scale_x_discrete(breaks=c("MM074", "MM087", "MM057", "MM057_DropSeq"),
                      labels=c("MM074 10x", "MM087 10x", "MM057 10x", "MM057 DropSeq")) +


	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf("auc_selectedSigs_vioplot.pdf")
for (i in c(1:(length(levels(auc$signature))/3))) {
  print(paste0("plotting page ", i))
  print(p + facet_grid_paginate(signature ~ ., nrow = 3, ncol = 1, page = i, scales = "free", drop = F))
}
dev.off()
#}

sessionInfo() #{
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)
#
#Matrix products: default
#BLAS/LAPACK: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/lib
#openblas_haswellp-r0.2.20.so
#
#locale:
#[1] C
#
#attached base packages:
#[1] grid      stats     graphics  grDevices utils     datasets  methods
#[8] base
#
#other attached packages:
#[1] hdf5r_1.2.0          ggforce_0.3.1        ggplot2_3.2.1
#[4] circlize_0.4.6       ComplexHeatmap_2.0.0 AUCell_1.6.1
#[7] data.table_1.12.8
#
#loaded via a namespace (and not attached):
# [1] Biobase_2.44.0              bit64_0.9-7
# [3] R.utils_2.9.0               shiny_1.3.2
# [5] assertthat_0.2.1            stats4_3.6.1
# [7] blob_1.2.0                  GenomeInfoDbData_1.2.1
# [9] pillar_1.4.2                RSQLite_2.1.2
#[11] backports_1.1.5             lattice_0.20-38
#[13] glue_1.3.1                  digest_0.6.21
#[15] polyclip_1.10-0             GenomicRanges_1.36.1
#[17] RColorBrewer_1.1-2          promises_1.0.1
#[19] XVector_0.24.0              colorspace_1.4-1
#[21] plyr_1.8.4                  htmltools_0.3.6
#[23] httpuv_1.5.2                Matrix_1.2-17
#[25] R.oo_1.22.0                 GSEABase_1.46.0
#[27] XML_3.99-0.3                pkgconfig_2.0.3
#[29] GetoptLong_0.1.7            zlibbioc_1.30.0
#[31] purrr_0.3.2                 xtable_1.8-4
#[33] scales_1.0.0                tweenr_1.0.1
#[35] later_0.8.0                 BiocParallel_1.18.1
#[37] tibble_2.1.3                annotate_1.62.0
#[39] farver_2.0.1                IRanges_2.18.3
#[41] withr_2.1.2                 SummarizedExperiment_1.14.1
#[43] BiocGenerics_0.30.0         lazyeval_0.2.2
#[45] magrittr_1.5                crayon_1.3.4
#[47] mime_0.7                    memoise_1.1.0
#[49] R.methodsS3_1.7.1           MASS_7.3-51.4
#[51] graph_1.62.0                tools_3.6.1
#[53] GlobalOptions_0.1.0         matrixStats_0.55.0
#[55] stringr_1.4.0               S4Vectors_0.22.1
#[57] munsell_0.5.0               cluster_2.1.0
#[59] DelayedArray_0.10.0         AnnotationDbi_1.46.1
#[61] compiler_3.6.1              GenomeInfoDb_1.20.0
#[63] rlang_0.4.0                 RCurl_1.95-4.12
#[65] rjson_0.2.20                labeling_0.3
#[67] tcltk_3.6.1                 bitops_1.0-6
#[69] gtable_0.3.0                reshape_0.8.8
#[71] DBI_1.0.0                   reshape2_1.4.3
#[73] R6_2.4.1                    dplyr_0.8.3
#[75] bit_1.1-15.2                zeallot_0.1.0
#[77] clue_0.3-57                 shape_1.4.4
#[79] stringi_1.4.3               parallel_3.6.1
#[81] Rcpp_1.0.3                  vctrs_0.2.0
#[83] png_0.1-7                   tidyselect_0.2.5
#}
