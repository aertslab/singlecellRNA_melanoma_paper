# question:
#are gene expression values significantly different between cell lines belonging 
#to the melanocytic, the intermediate and the mesenchymal-like cell state?

# task: run seurat differential expression tests (wilcoxon & MAST) on subset ten baselines
#results added as asterisks to Fig1

#module load R/3.5.0-iomkl-2018a-X11-20180131
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#module load binutils/2.28-GCCcore-6.4.0
#R
options(stringsAsFactors=FALSE)
.libPaths("/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/3.5")
staging <- "/ddn1/vol1/staging/leuven/stg_00002/lcb"
#############################################################################################
library(data.table)
library(Seurat)
seed <- 123

sample <- "12.ten_lines_BL"
zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices")
kdir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices")

wd <- file.path(kdir, sample)
setwd(wd)
dir.create("10.TextData/markers")

# get raw counts
indir <- file.path(zdir, sample)
counts <- read.table(file.path(indir, "10.TextData/raw_counts.tsv"))
meta <- read.table(file.path(indir, "10.TextData/meta_data.tsv"))

# remove doublets
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/10.analysis_with_all_data/00.RData/demuxlet_doublets.RData"), verbose = TRUE)
counts <- counts[, !colnames(counts) %in% demuxlet_doublets]
meta <- meta[colnames(counts),]

groups <- setNames(c(rep("mel", 3), rep("int", 3), "A375", rep("mes", 3)), 
	c("MM001", "MM011", "MM033", "MM057","MM074","MM087", "A375","MM029","MM047","MM099"))
meta[, "group"] <- groups[meta$Cell_Line]
table(meta$group)
#A375  int  mes  mel
# 566 1309 1242  745 

so <- CreateSeuratObject(counts = counts,
	min.cells = 0, min.features = 1000,
	meta.data = meta)
DefaultAssay(so)	#RNA
so <- SCTransform(so, seed.use = seed)
DefaultAssay(so)	#SCT
saveRDS(so, file = "00.RData/seurat_obj_diffex.rds.gz", compress = "gzip")

# define contrasts
################## {
contrasts <- list()
  meta$contr <- 0
  contr <- meta[, "contr", drop = F]

# test1: mesenchymal-like vs (melanocytic + intermediate)
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% 
	c(names(groups)[groups == "mes"], names(groups)[groups == "int"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mes"]), "contr"] <- 2

  contrasts[["mesenchymal_like_vs_melanocytic_and_intermediate"]] <- contr

# test2: intermediate vs (melanocytic & mesenchymal-like)
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% 
	c(names(groups)[groups == "mel"], names(groups)[groups == "mes"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "int"]), "contr"] <- 2

  contrasts[["intermediate_vs_mesenchymal_like_and_melanocytic"]] <- contr

# test3: melanocytic vs (intermediate + mesenchymal-like)
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% 
	c(names(groups)[groups == "mes"], names(groups)[groups == "int"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mel"]), "contr"] <- 2

  contrasts[["melanocytic_vs_intermediate_and_mesenchymal_like"]] <- contr

# test4: intermediate vs melanocytic 
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mel"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "int"]), "contr"] <- 2

  contrasts[["intermediate_vs_melanocytic"]] <- contr

# test5: mesenchymal-like vs intermediate 
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% c(names(groups)[groups == "int"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mes"]), "contr"] <- 2

  contrasts[["mesenchymal_like_vs_intermediate"]] <- contr

# test6: mesenchymal-like vs melanocytic 
  contr[, "contr"] <- 0
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mel"]), "contr"] <- 1
  contr[meta$Cell_Line %in% c(names(groups)[groups == "mes"]), "contr"] <- 2

  contrasts[["mesenchymal_like_vs_melanocytic"]] <- contr
#}
 
# wilcox
######## {
markers_wilcox <- sapply(contrasts, function(contr) {
	so <- AddMetaData(so, metadata = contr$contr, col.name = "contr")
	Idents(so) <- "contr"

	m <- FindMarkers(so, ident.1 = 2, ident.2 = 1,
		only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
	m
}, simplify = FALSE, USE.NAMES = TRUE)
saveRDS(markers_wilcox, file = "10.TextData/markers/markers_wilcox.rds.gz", compress = "gzip")

genes <- c("MITF", "TYR", "TGFBI", "SERPINE1", "FN1", "S100A16", "IFITM3", "HLA-B", "NES", "MIA")
lapply(markers_wilcox, function(m) { 
	print(m[genes,])
})

lapply(names(markers_wilcox), function(name) {
	write.table(markers_wilcox[[name]], 
		file = paste0("10.TextData/markers/markers_", name, "_wilcox.tsv"),
		quote = F, sep = "\t")
})
#}

# MAST
###### {
markers_mast <- sapply(contrasts, function(contr) {
	so <- AddMetaData(so, metadata = contr$contr, col.name = "contr")
	Idents(so) <- "contr"

	m <- FindMarkers(so, ident.1 = 2, ident.2 = 1,
		only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
		test.use = "MAST")
	m
}, simplify = FALSE, USE.NAMES = TRUE)
saveRDS(markers_mast, file.path("10.TextData/markers/markers_mast.rds.gz"), compress = "gzip")

lapply(names(markers_mast), function(name) {
	write.table(markers_mast[[name]], 
		file = paste0("10.TextData/markers/markers_", name, "_mast.tsv"),
		quote = F, sep = "\t")
})

genes <- c("MITF", "TYR", "TGFBI", "SERPINE1", "FN1", "S100A16", "IFITM3", "HLA-B", "NES", "MIA")
lapply(markers_mast, function(m) { 
	print(m[genes,])
})
#}

## genes used in fig1
genes <- c("MITF","TYR","TGFBI","SERPINE1", "S100A16", "FN1", "IFITM3", "HLA-B", "NES", "MIA")
lapply(markers_wilcox, function(m) { 
	print(m[genes,])
})
#{
#$mesenchymal_like_vs_melanocytic_and_intermediate
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      0.000000e+00 -1.4164035 0.312 0.935  0.000000e+00
#TYR       0.000000e+00 -1.7867011 0.024 0.904  0.000000e+00
#TGFBI    2.289276e-277  1.0725137 0.798 0.095 4.270872e-273
#SERPINE1 5.737141e-231  1.2021843 0.622 0.021 1.070321e-226
#S100A16  4.463208e-294  1.3368108 0.947 0.406 8.326562e-290
#FN1       7.282611e-14 -0.2858948 0.880 0.909  1.358644e-09
#IFITM3   1.316495e-108  0.7497035 0.933 0.698 2.456052e-104
#HLA-B     4.845715e-22 -0.5090694 0.850 0.788  9.040165e-18
#NES      5.643973e-126 -0.6007241 0.065 0.481 1.052940e-121
#MIA      7.840516e-265 -1.1346252 0.028 0.704 1.462727e-260
#
#$intermediate_vs_mesenchymal_like_and_melanocytic
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      5.998592e-68  0.2743267 0.935 0.565  1.119097e-63
#TYR       0.000000e+00  1.2947778 0.904 0.277  0.000000e+00
#TGFBI    9.211865e-135 -0.7688894 0.095 0.503 1.718565e-130
#SERPINE1 2.740292e-131 -0.8968356 0.021 0.394 5.112288e-127
#S100A16   4.635591e-87 -0.9606149 0.406 0.620  8.648158e-83
#FN1       2.921450e-97  0.6652736 0.909 0.641  5.450258e-93
#IFITM3    1.722634e-03 -0.3768134 0.698 0.603  1.000000e+00
#HLA-B     7.325953e-98  0.8462176 0.788 0.565  1.366730e-93
#NES      1.810688e-165  0.5820009 0.481 0.076 3.378020e-161
#MIA       0.000000e+00  1.1251787 0.704 0.035  0.000000e+00
#
#$melanocytic_vs_intermediate_and_mesenchymal_like
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     1.566120e-192  0.9442321 0.985 0.632 2.921754e-188
*#TYR       2.738308e-02 -0.2756440 0.698 0.475  1.000000e+00*
#TGFBI     3.502696e-99 -0.7868719 0.012 0.437  6.534629e-95
#SERPINE1  1.904767e-62 -0.7654194 0.013 0.314  3.553534e-58
#S100A16  2.557218e-168 -1.3339038 0.074 0.670 4.770747e-164
#FN1      9.434850e-266 -2.0022926 0.243 0.895 1.760166e-261
#IFITM3   3.968855e-251 -1.4572845 0.052 0.813 7.404295e-247
#HLA-B    3.416737e-245 -1.7341535 0.090 0.819 6.374265e-241
#NES       1.689389e-26 -0.3032610 0.095 0.279  3.151723e-22
#MIA       8.019649e-68 -0.7086720 0.047 0.375  1.496146e-63
#
#$intermediate_vs_melanocytic
#                p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     7.744261e-69 -0.4841612 0.935 0.985  1.444769e-64
#TYR      2.144037e-97  0.7953619 0.904 0.698  3.999915e-93
#NA                 NA         NA    NA    NA            NA
#NA.1               NA         NA    NA    NA            NA
#S100A16  3.093745e-59  0.4724574 0.406 0.074  5.771691e-55
#FN1     9.167957e-236  2.1313356 0.909 0.243 1.710374e-231
#IFITM3  6.583633e-167  1.0232151 0.698 0.052 1.228243e-162
#HLA-B   1.711595e-200  1.9501150 0.788 0.090 3.193152e-196
#NES      2.896504e-72  0.5515462 0.481 0.095  5.403717e-68
#MIA     2.889411e-171  1.1096260 0.704 0.047 5.390485e-167
#
#$mesenchymal_like_vs_intermediate
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      0.000000e+00 -1.4164035 0.312 0.935  0.000000e+00
#TYR       0.000000e+00 -1.7867011 0.024 0.904  0.000000e+00
#TGFBI    2.289276e-277  1.0725137 0.798 0.095 4.270872e-273
#SERPINE1 5.737141e-231  1.2021843 0.622 0.021 1.070321e-226
#S100A16  4.463208e-294  1.3368108 0.947 0.406 8.326562e-290
#FN1       7.282611e-14 -0.2858948 0.880 0.909  1.358644e-09
#IFITM3   1.316495e-108  0.7497035 0.933 0.698 2.456052e-104
#HLA-B     4.845715e-22 -0.5090694 0.850 0.788  9.040165e-18
#NES      5.643973e-126 -0.6007241 0.065 0.481 1.052940e-121
#MIA      7.840516e-265 -1.1346252 0.028 0.704 1.462727e-260
#
#$mesenchymal_like_vs_melanocytic
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     1.063684e-303 -1.9005647 0.312 0.985 1.984409e-299
#TYR      6.721751e-231 -0.9913392 0.024 0.698 1.254010e-226
#TGFBI    5.729068e-221  1.1987052 0.798 0.012 1.068815e-216
#SERPINE1 4.133654e-150  1.2100156 0.622 0.013 7.711746e-146
#S100A16  3.963438e-282  1.8092682 0.947 0.074 7.394191e-278
#FN1      1.122833e-209  1.8454407 0.880 0.243 2.094757e-205
#IFITM3   3.168118e-278  1.7729186 0.933 0.052 5.910441e-274
#HLA-B    1.208818e-222  1.4410456 0.850 0.090 2.255171e-218
#NA                  NA         NA    NA    NA            NA
#NA.1                NA         NA    NA    NA            NA
##}

lapply(markers_mast, function(m) { 
	print(m[genes,])
})
#{
#$mesenchymal_like_vs_melanocytic_and_intermediate
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      0.000000e+00 -1.4164035 0.312 0.935  0.000000e+00
#TYR       0.000000e+00 -1.7867011 0.024 0.904  0.000000e+00
#TGFBI    6.916919e-323  1.0725137 0.798 0.095 1.290420e-318
#SERPINE1 8.603083e-280  1.2021843 0.622 0.021 1.604991e-275
#S100A16   0.000000e+00  1.3368108 0.947 0.406  0.000000e+00
#FN1       2.323037e-13 -0.2858948 0.880 0.909  4.333859e-09
#IFITM3   1.411515e-117  0.7497035 0.933 0.698 2.633322e-113
#HLA-B     4.961597e-70 -0.5090694 0.850 0.788  9.256356e-66
#NES      3.459284e-143 -0.6007241 0.065 0.481 6.453639e-139
#MIA       0.000000e+00 -1.1346252 0.028 0.704  0.000000e+00
#
#$intermediate_vs_mesenchymal_like_and_melanocytic
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     5.368672e-133  0.2743267 0.935 0.565 1.001579e-128
#TYR       0.000000e+00  1.2947778 0.904 0.277  0.000000e+00
#TGFBI    1.846750e-155 -0.7688894 0.095 0.503 3.445297e-151
#SERPINE1 9.454050e-168 -0.8968356 0.021 0.394 1.763748e-163
#S100A16  6.259322e-162 -0.9606149 0.406 0.620 1.167739e-157
#FN1      1.780351e-102  0.6652736 0.909 0.641  3.321423e-98
#IFITM3    5.248978e-62 -0.3768134 0.698 0.603  9.792493e-58
#HLA-B    1.783666e-117  0.8462176 0.788 0.565 3.327607e-113
#NES      1.405041e-172  0.5820009 0.481 0.076 2.621244e-168
#MIA       0.000000e+00  1.1251787 0.704 0.035  0.000000e+00
#
#$melanocytic_vs_intermediate_and_mesenchymal_like
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     1.072286e-219  0.9442321 0.985 0.632 2.000456e-215
*#TYR       2.203311e-94 -0.2756440 0.698 0.475  4.110497e-90*	???
#TGFBI    2.561174e-141 -0.7868719 0.012 0.437 4.778126e-137
#SERPINE1  1.454274e-88 -0.7654194 0.013 0.314  2.713094e-84
#S100A16  1.021875e-216 -1.3339038 0.074 0.670 1.906411e-212
#FN1       0.000000e+00 -2.0022926 0.243 0.895  0.000000e+00
#IFITM3    0.000000e+00 -1.4572845 0.052 0.813  0.000000e+00
#HLA-B     0.000000e+00 -1.7341535 0.090 0.819  0.000000e+00
#NES       2.822798e-31 -0.3032610 0.095 0.279  5.266213e-27
#MIA       1.424564e-90 -0.7086720 0.047 0.375  2.657666e-86
#
#$intermediate_vs_melanocytic
#                p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF     4.133334e-68 -0.4841612 0.935 0.985  7.711148e-64
#TYR     2.759923e-106  0.7953619 0.904 0.698 5.148913e-102
#NA                 NA         NA    NA    NA            NA
#NA.1               NA         NA    NA    NA            NA
#S100A16  3.199687e-68  0.4724574 0.406 0.074  5.969337e-64
#FN1     2.137121e-293  2.1313356 0.909 0.243 3.987012e-289
#IFITM3  5.901115e-211  1.0232151 0.698 0.052 1.100912e-206
#HLA-B   3.166253e-262  1.9501150 0.788 0.090 5.906962e-258
#NES      9.817477e-83  0.5515462 0.481 0.095  1.831548e-78
#MIA     5.654955e-220  1.1096260 0.704 0.047 1.054988e-215
#
#$mesenchymal_like_vs_intermediate
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      0.000000e+00 -1.4164035 0.312 0.935  0.000000e+00
#TYR       0.000000e+00 -1.7867011 0.024 0.904  0.000000e+00
#TGFBI    6.916919e-323  1.0725137 0.798 0.095 1.290420e-318
#SERPINE1 8.603083e-280  1.2021843 0.622 0.021 1.604991e-275
#S100A16   0.000000e+00  1.3368108 0.947 0.406  0.000000e+00
#FN1       2.323037e-13 -0.2858948 0.880 0.909  4.333859e-09
#IFITM3   1.411515e-117  0.7497035 0.933 0.698 2.633322e-113
#HLA-B     4.961597e-70 -0.5090694 0.850 0.788  9.256356e-66
#NES      3.459284e-143 -0.6007241 0.065 0.481 6.453639e-139
#MIA       0.000000e+00 -1.1346252 0.028 0.704  0.000000e+00
#
#$mesenchymal_like_vs_melanocytic
#                 p_val  avg_logFC pct.1 pct.2     p_val_adj
#MITF      0.000000e+00 -1.9005647 0.312 0.985  0.000000e+00
#TYR      1.932112e-256 -0.9913392 0.024 0.698 3.604548e-252
#TGFBI    1.160616e-308  1.1987052 0.798 0.012 2.165246e-304
#SERPINE1 1.356094e-200  1.2100156 0.622 0.013 2.529929e-196
#S100A16   0.000000e+00  1.8092682 0.947 0.074  0.000000e+00
#FN1      9.250795e-246  1.8454407 0.880 0.243 1.725828e-241
#IFITM3    0.000000e+00  1.7729186 0.933 0.052  0.000000e+00
#HLA-B    1.135414e-283  1.4410456 0.850 0.090 2.118228e-279
#NA                  NA         NA    NA    NA            NA
#NA.1                NA         NA    NA    NA            NA
#}

sessionInfo() #{
#R version 3.5.0 (2018-04-23)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)
#
#Matrix products: default
#BLAS: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/lib/libR.so
#LAPACK: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/modules/lapack.so
#
#locale:
#[1] C
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
#[1] Seurat_3.0.1.9010 data.table_1.12.2
#
#loaded via a namespace (and not attached):
# [1] httr_1.4.0          tidyr_0.8.3         jsonlite_1.6
# [4] viridisLite_0.3.0   splines_3.5.0       lsei_1.2-0
# [7] R.utils_2.9.0       gtools_3.8.1        Rdpack_0.11-0
#[10] assertthat_0.2.1    ggrepel_0.8.1       globals_0.12.4
#[13] pillar_1.4.1        lattice_0.20-38     glue_1.3.1
#[16] reticulate_1.12     digest_0.6.21       RColorBrewer_1.1-2
#[19] SDMTools_1.1-221.1  colorspace_1.4-1    cowplot_0.9.4
#[22] htmltools_0.3.6     Matrix_1.2-17       R.oo_1.22.0
#[25] plyr_1.8.4          pkgconfig_2.0.2     bibtex_0.4.2
#[28] tsne_0.1-3          listenv_0.7.0       purrr_0.3.2
#[31] scales_1.0.0        RANN_2.6.1          gdata_2.18.0
#[34] Rtsne_0.15          tibble_2.1.3        ggplot2_3.1.1
#[37] ROCR_1.0-7          pbapply_1.4-0       lazyeval_0.2.2
#[40] survival_2.44-1.1   magrittr_1.5        crayon_1.3.4
#[43] R.methodsS3_1.7.1   future_1.14.0       nlme_3.1-140
#[46] MASS_7.3-51.4       gplots_3.0.1.1      ica_1.0-2
#[49] tools_3.5.0         fitdistrplus_1.0-14 gbRd_0.4-11
#[52] stringr_1.4.0       plotly_4.9.0        munsell_0.5.0
#[55] cluster_2.0.9       irlba_2.3.3         compiler_3.5.0
#[58] rsvd_1.0.1          caTools_1.17.1.2    rlang_0.3.4
#[61] grid_3.5.0          ggridges_0.5.1      htmlwidgets_1.3
#[64] igraph_1.2.4.1      bitops_1.0-6        npsurv_0.4-0
#[67] gtable_0.3.0        codetools_0.2-16    reshape2_1.4.3
#[70] R6_2.4.0            gridExtra_2.3       zoo_1.8-6
#[73] dplyr_0.8.1         future.apply_1.2.0  KernSmooth_2.23-15
#[76] metap_1.1           ape_5.1             stringi_1.4.3
#[79] parallel_3.5.0      Rcpp_1.0.2          sctransform_0.2.0
#[82] png_0.1-7           tidyselect_0.2.5    lmtest_0.9-37
#}
