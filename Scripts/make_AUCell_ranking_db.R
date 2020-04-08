# task: make data base with AUCell rankings for each 10x subset

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
library(hdf5r)
library(AUCell)
seed <- 123

indir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/")

# for each subset, create rankings for AUCell {
for (ss in c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")) {
	wd <- file.path(indir, ss) 
	setwd(wd)
	
	loomfile <- "log_CPM_woDoublets.loom"
	loom <- H5File$new(loomfile, mode = "r")

	### retrieve gene expression matrix (genes x cells)
	dgem <- t(loom[["matrix"]]$read())

	### get genes and cells
	genes <- loom[["row_attrs"]][["Gene"]]$read()
	rownames(dgem) <- genes
	cells <- loom[["col_attrs"]][["CellID"]]$read()
	colnames(dgem) <- cells

	### close loom file
	loom$close_all()

	# create rankings 
	set.seed(seed)
	aucellRankings <- AUCell_buildRankings(dgem, nCores = 20, plotStats = F)
	saveRDS(aucellRankings, file = "aucellRankings.rds.gz", compress = "gzip")
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
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
#[1] hdf5r_1.2.0  AUCell_1.6.1
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3                  XVector_0.24.0
# [3] GenomeInfoDb_1.20.0         pillar_1.4.2
# [5] compiler_3.6.1              later_0.8.0
# [7] zlibbioc_1.30.0             tools_3.6.1
# [9] R.methodsS3_1.7.1           R.utils_2.9.0
#[11] bitops_1.0-6                zeallot_0.1.0
#[13] digest_0.6.21               bit_1.1-15.2
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
#[37] data.table_1.12.8           R6_2.4.1
#[39] AnnotationDbi_1.46.1        BiocParallel_1.18.1
#[41] XML_3.99-0.3                blob_1.2.0
#[43] magrittr_1.5                matrixStats_0.55.0
#[45] GenomicRanges_1.36.1        backports_1.1.5
#[47] promises_1.0.1              htmltools_0.3.6
#[49] BiocGenerics_0.30.0         SummarizedExperiment_1.14.1
#[51] mime_0.7                    xtable_1.8-4
#[53] httpuv_1.5.2                RCurl_1.95-4.12
#[55] crayon_1.3.4                R.oo_1.22.0
#}
