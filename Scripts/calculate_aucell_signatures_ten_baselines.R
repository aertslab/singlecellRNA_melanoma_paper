# task: calculate AUCell on gene signatures subset ten baselines & plot
#violin plots: Fig1e,f,h,i,l; FigS1a
#boxplot: Fig2c
#plot expression on tSNE/diffusion map: FigS1b

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

# functions
read.gmt <- function (fileName) {
	tmp <- unname(sapply(readLines(fileName), function(x) strsplit(as.character(x), "\t")))
	tmp <- tmp[which(lengths(tmp) > 0)]
	names(tmp) <- sapply(tmp, function(x) x[1])
	lapply(tmp, function(x) x[3:length(x)])
}
# scale each row from 0-1 (for heatmap)
scale01 <- function(x, low = min(x), high = max(x)) {
       x <- (x - low)/(high - low)
       x
}
# colours [[MM40kColours]]
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40kColours/MM40kColours.RData"))

wd <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/12.ten_lines_BL/20.Plots/")
dir.create(wd)
setwd(wd)

# calculate AUCell
################## {
# re-use ranking db from make_AUCell_ranking_db.R
rankdir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL")
aucellRankings <- readRDS(file.path(rankdir, "aucellRankings.rds.gz"))
#0c00ab9146f1d19227a309e974af008dfdaaee97	aucellRankings.rds.gz

# new signatures
signatures <- read.gmt(file.path(staging, "kspan/resources/genesets/melanoma_signatures_190618.gmt"))

auc <- AUCell_calcAUC(signatures, aucellRankings, 
		aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 20)
saveRDS(auc, file = file.path(rankdir, "auc_signatures.rds.gz"), compress = "gzip")

# sort aucell values according to cell lines
auc <- SummarizedExperiment::assays(auc)[["AUC"]]
load(file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/12.ten_lines_BL/00.RData/meta_data.RData"))
meta_data <- meta_data[rownames(meta_data) %in% colnames(auc),]
lines <- c("MM001", "MM011", "MM031", "MM074", "MM087", "MM057", "A375", "MM029", "MM047", "MM099")

meta_data$Cell_Line <- factor(meta_data$Cell_Line, levels = lines)
meta <- meta_data[order(meta_data$Cell_Line),]
meta$Cell_ID <- rownames(meta)
df <- meta[,"Cell_Line", drop = F]
#}

# plot heatmap FigS1a
##################### {
#select signatures
selsigs <- setNames(c(
	"VERF_INV", "VERF_PRO", "HOEK_INV", "HOEK_PRO", 
	"rambow_neuro", "rambow_immune", "rambow_pigmentation",
	"raj_resistance", "hugo_IPRES", "WU_CELL_MIGRATION"),
	
	c("Verfaillie invasive", "Verfaillie proliferative", 
	"Hoek invasive", "Hoek proliferative",
	"Rambow neuro", "Rambow immune", "Rambow pigmentation",
	"Raj resistance", "Hugo IPRES", "Wu cell migration"))

ps <- auc[as.vector(selsigs), rownames(df)]
rownames(ps) <- names(selsigs)
pss <- t(apply(ps, 1, scale01))

plotcolours <- MM40kColours["Cell_Line"]
ha <- HeatmapAnnotation(df = df[,"Cell_Line", drop = F], 
	col = plotcolours)
f1 <- colorRamp2(seq(min(pss), quantile(unlist(pss), 0.95), length = 2), 
		c("#EEEEEE", "red"))
h1 <- Heatmap(pss, col = f1, name = "AUCell\nscaled", cluster_columns = FALSE,
	show_column_names = FALSE, row_names_gp = gpar(fontsize = 8), top_annotation = ha)

pdf("heatmap_aucell_sc_selected_sigs.pdf", height = 2.5)
draw(h1)
dev.off()
#}

# violin plots Fig1e,f,h,i,l
############################ {
selsigs <- c("HOEK_INV", "HOEK_PRO", "GO_PIGMENTATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
	"AXL_SIGNATURE_GODING", "KARAKAS_TGFB1_SIGNALING", "KEGG_GRAFT_VERSUS_HOST_DISEASE",
	"rambow_immune", "rambow_neuro", "LEE_NEURAL_CREST_STEM_CELL_UP")

aucs <- data.table(reshape::melt(auc[selsigs,]))
setnames(aucs, c("signature", "Cell_ID", "AUCell"))
aucs <- merge(aucs, meta, by = "Cell_ID")[,.(Cell_ID, Cell_Line, Experiment, nGene, signature, AUCell)]

p <- ggplot(aucs, aes(x = Cell_Line, y = AUCell, fill = Cell_Line)) +
	geom_violin(scale = "width") +
	scale_fill_manual(values = MM40kColours$Cell_Line) +

   	geom_jitter(aes(x = Cell_Line, y = AUCell),
		shape=16, size = 0.3, width = 0.2) +

	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf("vioplots_aucell_sc_selected_sigs.pdf")
for (i in c(1:(length(levels(aucs$signature))))) {
  print(paste0("plotting page ", i))
  print(p + facet_wrap_paginate(signature ~ ., nrow = 1, ncol = 1, page = i, scales = "free", drop = F))
}
dev.off()
#}

# plot heatmaps of signature AUCells along PC & CCA axes Fig1h,j
################################################################ {
# get cell order
pca <- read.table(file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/12.ten_lines_BL/10.TextData/seurat_pc1_pc2.tsv"))
pca <- pca[colnames(auc),]

#[[MM40k_12baselines_JasperCCA.Rmd]]
cca <- read.table(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/12.ten_lines_BL/10.TextData/CCA_jasper_CCA1_CCA2.tsv"))
cca <- cca[colnames(auc),]

#same colour gradient across all values
f1 <- colorRamp2(c(0, 0.8), c("#EEEEEE", "red"))

##### PC1 
scaled_auc <- t(apply(auc[c("HOEK_PRO", "HOEK_INV"),], 1, scale01))

#sort by PC1
cells <- rownames(pca)[order(pca$PC1)]
sauc <- scaled_auc[,cells]

h1 <- Heatmap(sauc, col = f1, 
	name = "AUCell\nscaled", 
	cluster_columns = FALSE, 
	cluster_rows = F, 
	show_column_names = F)

png("heatmap_aucell_sc_hoek_PC1.png", width = 1000, height = 100)
draw(h1)
dev.off()

pdf("heatmap_aucell_sc_hoek_PC1.pdf", width = 10, height = 2)
draw(h1)
dev.off()


#### CC1
cells <- rownames(cca)[order(cca$CC1)]
sauc <- scaled_auc[,cells]

h1 <- Heatmap(sauc, col = f1, 
	name = "AUCell\nscaled", 
	cluster_columns = FALSE, 
	cluster_rows = F, 
	show_column_names = F)

png("heatmap_aucell_sc_hoek_CC1.png", width = 1000, height = 100)
draw(h1)
dev.off()

pdf("heatmap_aucell_sc_hoek_CC1.pdf", width = 10, height = 2)
draw(h1)
dev.off()

#### PC2
scaled_auc <- scale01(data.frame(auc["KEGG_ALLOGRAFT_REJECTION",,drop = F]))
rownames(scaled_auc) <- c("immII")
cells <- rownames(pca)[order(pca$PC2)]
sauc <- as.matrix(scaled_auc[,cells])

h1 <- Heatmap(sauc, col = f1, 
	name = "AUCell\nscaled", 
	cluster_columns = FALSE, 
	cluster_rows = F, 
	show_column_names = F)

png("heatmap_aucell_sc_immune_PC2.png", width = 1000, height = 100)
draw(h1)
dev.off()

pdf("heatmap_aucell_sc_immune_PC2.pdf", width = 10, height = 2)
draw(h1)
dev.off()

#### CC2
scaled_auc <- scale01(data.frame(auc["rambow_mitosis",,drop = F]))
rownames(scaled_auc) <- c("mitosis")
cells <- rownames(cca)[order(cca$CC2)]
sauc <- as.matrix(scaled_auc[,cells])

h1 <- Heatmap(sauc, col = f1, 
	name = "AUCell\nscaled", 
	cluster_columns = FALSE, 
	cluster_rows = F, 
	show_column_names = F)

png("heatmap_aucell_sc_mitosis_CC2.png", width = 1000, height = 100)
draw(h1)
dev.off()

pdf("heatmap_aucell_sc_mitosis_CC2.pdf", width = 10, height = 2)
draw(h1)
dev.off()
#}

# boxplot migration signatures Wu et al Fig2c
############################################# {
sign <- "WU_CELL_MIGRATION"

aucs <- data.table(reshape::melt(auc[sign,,drop = F]))
setnames(aucs, c("signature", "Cell_ID", "AUCell"))
aucs <- merge(aucs, meta, by = "Cell_ID")[,.(Cell_ID, Cell_Line, Experiment, nGene, signature, AUCell)]

p <- ggplot(aucs, aes(x = Cell_Line, y = AUCell, fill = Cell_Line)) +
	scale_fill_manual(values = MM40kColours$Cell_Line) +
	scale_colour_manual(values = MM40kColours$Cell_Line) +
	theme_classic() +
	guides(fill=FALSE) +    #rm legend
	labs(x = "cell line", y = "AUCell", title = "Migratory potential Wu et al.)") +
	geom_boxplot() 

pdf("aucell_sc_boxpl_WU_CELL_MIGRATION.pdf")
print(p)
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
#[4] ComplexHeatmap_2.0.0 AUCell_1.6.1         data.table_1.12.8
#
#loaded via a namespace (and not attached):
#[1] Biobase_2.44.0              bit64_0.9-7
#[3] R.utils_2.9.0               shiny_1.3.2
#[5] assertthat_0.2.1            stats4_3.6.1
#[7] blob_1.2.0                  GenomeInfoDbData_1.2.1
#[9] pillar_1.4.2                RSQLite_2.1.2
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
#[55] S4Vectors_0.22.1            munsell_0.5.0
#[57] cluster_2.1.0               DelayedArray_0.10.0
#[59] AnnotationDbi_1.46.1        compiler_3.6.1
#[61] GenomeInfoDb_1.20.0         rlang_0.4.0
#[63] RCurl_1.95-4.12             rjson_0.2.20
#[65] labeling_0.3                bitops_1.0-6
#[67] gtable_0.3.0                reshape_0.8.8
#[69] DBI_1.0.0                   R6_2.4.1
#[71] dplyr_0.8.3                 bit_1.1-15.2
#[73] zeallot_0.1.0               clue_0.3-57
#[75] shape_1.4.4                 parallel_3.6.1
#[77] Rcpp_1.0.3                  vctrs_0.2.0
#[79] png_0.1-7                   tidyselect_0.2.5
#}
