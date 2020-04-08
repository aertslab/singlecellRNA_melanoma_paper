# question:
#are AUCell values calculated on the filtered regulons of the ten baselines
#significantly different between cell lines belonging 
#to the melanocytic, the intermediate and the mesenchymal-like cell state?

# task: apply wilcoxon rank sum test on AUCell values & calculate bonferroni FDR

#module load R/3.6.1-foss-2018a-X11-20180604
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#PATH=$PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/bin/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/lib/
#R
options(stringsAsFactors=FALSE)
options(bitmapType="cairo")
.libPaths("/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/3.6")
staging <- "/ddn1/vol1/staging/leuven/stg_00002/lcb"
#############################################################################################
library(reshape2)
library(data.table)

wd <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL")
setwd(wd)

# aucell values from ten_baselines_selected_regulons_aucell_heatmap_and_violinplots.R
auc <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_aucell01.rds.gz")

# meta data
zdir <- file.path(staging, "zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices")
meta <- read.table(file.path(zdir, "12.ten_lines_BL/10.TextData/meta_data.tsv"))
load(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/10.analysis_with_all_data/00.RData/demuxlet_doublets.RData"), verbose = TRUE)
meta <- meta[!rownames(meta) %in% demuxlet_doublets,]

# remove A375 completely
meta <- meta[meta$Cell_Line != "A375",]
auc <- auc[, rownames(meta)]

# melt auc data
tt <- reshape2::melt(auc)
tt <- merge(tt, meta, by.x = "cells", by.y = "row.names")
tt <- subset(tt, select=c("gene sets","value","Cell_Line"))
setDT(tt)
colnames(tt) <- c("regulon", "auc", "cellline")

# define cell state groups
groups <- list(mel = c("MM001", "MM011", "MM033"),
	int = c("MM057","MM074","MM087"),
	mes = c("MM029","MM047","MM099"))

# one group vs all others
################################ {
# test1: mesenchymal-like vs all others
  tt[, phenotype := ifelse(cellline %in% groups[["mes"]], "mes", "rest")]
  tt[, mes_vs_all_mel := wilcox.test(auc ~ phenotype)$p.value, by = regulon]
 
# test2: intermediate vs rest
  tt[, phenotype := ifelse(cellline %in% groups[["int"]], "int", "rest")]
  tt[, int_vs_mes_and_mel_only := wilcox.test(auc ~ phenotype)$p.value, by = regulon]
 
# test3: melanocytic vs rest
  tt[, phenotype := ifelse(cellline %in% groups[["mel"]], "mel", "rest")]
  tt[, mel_only_vs_intermediate_and_mes := wilcox.test(auc ~ phenotype)$p.value, by = regulon]
#}

# contrast groups
######################## {
# test1: melanocytic vs intermediate
  tt[, phenotype := ifelse(cellline %in% groups[["mel"]], "mel", "rest")]
  tt[cellline %in% c(groups[["mel"]], groups[["int"]]), 
	mel_only_vs_intermediate := wilcox.test(auc ~ phenotype)$p.value, by = regulon]

# test2: melanocytic vs mesenchymal-like
  tt[, phenotype := ifelse(cellline %in% groups[["mel"]], "mel", "rest")]
  tt[cellline %in% c(groups[["mel"]], groups[["mes"]]), 
	mel_only_vs_mes := wilcox.test(auc ~ phenotype)$p.value, by = regulon]

# test3: intermediate vs mesenchymal-like
  tt[, phenotype := ifelse(cellline %in% groups[["int"]], "int", "rest")]
  tt[cellline %in% c(groups[["int"]], groups[["mes"]]), 
	intermediate_vs_mes := wilcox.test(auc ~ phenotype)$p.value, by = regulon]

res <- unique(tt[, .(regulon, 
	mes_vs_all_mel, int_vs_mes_and_mel_only, mel_only_vs_intermediate_and_mes,
		mel_only_vs_intermediate, mel_only_vs_mes, intermediate_vs_mes)])
#}

#remove NA values
res <- res[, .(
	mes_vs_all_mel = mean(mes_vs_all_mel, na.rm = TRUE),
	int_vs_mes_and_mel_only = mean(int_vs_mes_and_mel_only, na.rm = TRUE),
	mel_only_vs_intermediate_and_mes = mean(mel_only_vs_intermediate_and_mes, na.rm = TRUE),

	mel_only_vs_intermediate = mean(mel_only_vs_intermediate, na.rm = TRUE), 
	mel_only_vs_mes = mean(mel_only_vs_mes, na.rm = TRUE), 
	intermediate_vs_mes = mean(intermediate_vs_mes, na.rm = TRUE)), .(regulon)]

#multiple testing correction
res[, mes_vs_all_mel_padj := p.adjust(mes_vs_all_mel, method = "bonferroni")]
res[, int_vs_mes_and_mel_only_padj := p.adjust(int_vs_mes_and_mel_only, method = "bonferroni")]
res[, mel_only_vs_intermediate_and_mes_padj := p.adjust(mel_only_vs_intermediate_and_mes, method = "bonferroni")]
res[, mel_only_vs_intermediate_padj := p.adjust(mel_only_vs_intermediate, method = "bonferroni")]
res[, mel_only_vs_mes_padj := p.adjust(mel_only_vs_mes, method = "bonferroni")]
res[, intermediate_vs_mes_padj := p.adjust(intermediate_vs_mes, method = "bonferroni")]

res <- res[, .(regulon, 
	mes_vs_all_mel, mes_vs_all_mel_padj,
	int_vs_mes_and_mel_only, int_vs_mes_and_mel_only_padj,
	mel_only_vs_intermediate_and_mes, mel_only_vs_intermediate_and_mes_padj,

	mel_only_vs_intermediate, mel_only_vs_intermediate_padj,
	mel_only_vs_mes, mel_only_vs_mes_padj,
	intermediate_vs_mes, intermediate_vs_mes_padj)]
#}
fwrite(res, file = "wilcoxon_rank_sum_test_selectedRegulons_tenBaselines.tsv",
	quote = F, sep = "\t")

#regulons used in figS4
regs <- c("MITF_mtf", "SOX10_mtf", "MITF_trk", "SOX10_trk")
res[regulon %in% regs, .(regulon, 
	mel_only_vs_intermediate_padj, mel_only_vs_mes_padj, intermediate_vs_mes_padj,
	mes_vs_all_mel_padj, int_vs_mes_and_mel_only_padj, mel_only_vs_intermediate_and_mes_padj)]

#     regulon mel_only_vs_intermediate_padj mel_only_vs_mes_padj
#	1:  MITF_mtf                 1.545426e-192        4.979380e-303
#2: SOX10_mtf                  3.338674e-50        4.979313e-303
#3: SOX10_trk                  4.315096e-73        5.534206e-303
#4:  MITF_trk                 1.456675e-206        4.979340e-303
#   intermediate_vs_mes_padj mes_vs_all_mel_padj int_vs_mes_and_mel_only_padj
#1:                        0                   0                 2.108662e-08
#2:                        0                   0                 0.000000e+00
#3:                        0                   0                 0.000000e+00
#4:                        0                   0                 5.348077e-08
#   mel_only_vs_intermediate_and_mes_padj
#1:                         3.540938e-246
#2:                          8.915883e-43
#3:                          1.292764e-29
#4:                         1.543612e-259

#regulons used in fig3
regs <- c("SOX10_mtf", "NFATC2_mtf", "EGR3_mtf", "ELF1_mtf", "ETV4_mtf", "SOX6_mtf")
res[regulon %in% regs, .(regulon, 
	mel_only_vs_intermediate_padj,
	intermediate_vs_mes_padj,
	int_vs_mes_and_mel_only_padj)]

#      regulon mel_only_vs_intermediate_padj intermediate_vs_mes_padj
#1:   EGR3_mtf                 6.957415e-259             0.000000e+00
#2:   ELF1_mtf                 1.536243e-275             0.000000e+00
#3:  SOX10_mtf                  3.338674e-50             0.000000e+00
#4:   ETV4_mtf                  0.000000e+00             0.000000e+00
#5: NFATC2_mtf                 2.516774e-282            2.691582e-258
#6:   SOX6_mtf                 4.421458e-220             0.000000e+00
#   int_vs_mes_and_mel_only_padj
#1:                 0.000000e+00
#2:                9.980695e-290
#3:                 0.000000e+00
#4:                7.722521e-279
#5:                 0.000000e+00
#6:                 0.000000e+00

sessionInfo() #{
# version 3.6.1 (2019-07-05)
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
#[1] AUCell_1.6.1      data.table_1.12.8 reshape2_1.4.3
#
#loaded via a namespace (and not attached):
# [1] SummarizedExperiment_1.14.1 lattice_0.20-38
# [3] vctrs_0.2.0                 htmltools_0.3.6
# [5] stats4_3.6.1                blob_1.2.0
# [7] XML_3.99-0.3                rlang_0.4.0
# [9] R.oo_1.22.0                 pillar_1.4.2
#[11] later_0.8.0                 DBI_1.0.0
#[13] R.utils_2.9.0               BiocParallel_1.18.1
#[15] BiocGenerics_0.30.0         bit64_0.9-7
#[17] matrixStats_0.55.0          GenomeInfoDbData_1.2.1
#[19] plyr_1.8.4                  stringr_1.4.0
#[21] zlibbioc_1.30.0             R.methodsS3_1.7.1
#[23] memoise_1.1.0               Biobase_2.44.0
#[25] IRanges_2.18.3              httpuv_1.5.2
#[27] GenomeInfoDb_1.20.0         parallel_3.6.1
#[29] AnnotationDbi_1.46.1        GSEABase_1.46.0
#[31] Rcpp_1.0.3                  xtable_1.8-4
#[33] backports_1.1.5             promises_1.0.1
#[35] DelayedArray_0.10.0         S4Vectors_0.22.1
#[37] graph_1.62.0                annotate_1.62.0
#[39] XVector_0.24.0              mime_0.7
#[41] bit_1.1-15.2                digest_0.6.21
#[43] stringi_1.4.3               shiny_1.3.2
#[45] GenomicRanges_1.36.1        grid_3.6.1
#[47] tools_3.6.1                 bitops_1.0-6
#[49] magrittr_1.5                RCurl_1.95-4.12
#[51] tibble_2.1.3                RSQLite_2.1.2
#[53] crayon_1.3.4                pkgconfig_2.0.3
#[55] zeallot_0.1.0               Matrix_1.2-17
#[57] R6_2.4.1                    compiler_3.6.1
#}
