#<[[MM40k100xSCENIC]]

# task: score tumor biopsy data from Jerby-Arnon et al. 2018 (GSE115978)
# with regulons from subset '10 baselines'
 
# [1] take SKCM loom file
# [2] score with combined regulons baselines
# [3] make new loom files with AUCell values for SCope 

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
library(rjson)
library(SCopeLoomR)	#v0.5.0
seed <- 123

wd <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40kCancerCompRescore")
dir.create(wd, recursive = TRUE)
setwd(wd)

# selected regulons from 10 melanoma lines baselines
regfile <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets/12.ten_lines_BL/all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_regulons.rds.gz")
#776c4d992c34da0369e8f129d5449e9ce077b69c all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_regulons.rds.gz
regs <- readRDS(regfile)

# jerby-arnon data
gset <- "SKCM_GSE115978"

# get dgem from loom file
loomfile <- paste0(gset, ".loom")
loom <- H5File$new(loomfile, mode = "r")
genes <- loom[["row_attrs"]][["Gene"]]$read()
cells <- loom[["col_attrs"]][["CellID"]]$read()
dgem <- t(loom[["matrix"]]$read())
rownames(dgem) <- genes
colnames(dgem) <- cells
emb <- loom[["col_attrs"]][["Embedding"]]$read()
rownames(emb) <- cells

# get meta data
meta <- data.frame(
	row.names = loom[["col_attrs"]][["CellID"]]$read()
	, ClusterID = loom[["col_attrs"]][["ClusterID"]]$read()
	, nGenes = loom[["col_attrs"]][["nGene"]]$read())
meta$clusters <- paste0("cluster_", meta$ClusterID)

#make new md w/o regulon thresholds
md <- get_global_meta_data(loom)
md[[4]] <- NULL
loom$close_all()

#make annotation table to add as meta data
annot <- rbindlist(lapply(md[["clusterings"]][[1]][["clusters"]], function(cl) { return(data.frame(cl)) }) , fill = TRUE)
annot[, description := gsub(" ", "_", description)]
annot[, description := gsub("\\?", "Unknown", description)]
annot <- setNames(annot[,description], paste0("cluster_", annot[, id]))
meta[, "celltype"] <- as.vector(annot[meta$clusters])

# add more metadata
annotfile <- file.path(staging, "zkalender/scMelanoma_paper/Jerby_Arnon_scRNAseq/GSE115978_cell.annotations.csv")
annot2 <- read.table(annotfile, sep = ",", header = T, row.names = 1)
meta[, "sample"] <- annot2[rownames(meta), "samples"]
meta[, "treatment_group"] <- annot2[rownames(meta), "treatment.group"]
meta[, "cohort"] <- annot2[rownames(meta), "Cohort"]
meta[, "no_of_reads"] <- annot2[rownames(meta), "no.of.reads"]

# and make plot with cell type annotations (Fig5b)
library(wesanderson)
colours <- setNames(c(wes_palette("Darjeeling1"),
	"#00C957", wes_palette("Royal1")[2], wes_palette("Royal2")[3], 
	wes_palette("Cavalcanti1")[2], "slateblue3"),
	unique(meta$celltype))

pdf(paste0(gset, "_tsne.pdf"))
	plot(emb, col = colours[meta[rownames(emb), "celltype"]],
		pch = 16, cex = 0.9, main = "cell line")
	plot.new()
	legend("top", names(colours), col = colours,
		pch = 20, cex = 1, pt.cex = 3)

	#CD4+: T helper
	#CD8+: cytotoxic
	#CD3+: T cell
	#malignant cells (S100+, MITF+)

# plot sample annotation (32)
colours <- setNames(c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"),
	wes_palette("BottleRocket2"), wes_palette("Zissou1"),
	wes_palette("Royal2"), wes_palette("Cavalcanti1"), wes_palette("IsleofDogs2")[2:3]),
	unique(meta$sample))
plot(emb, col = colours[meta[rownames(emb), "sample"]],
		pch = 16, cex = 0.9, main = "sample")
	plot.new()
	legend("top", names(colours), col = colours,
		pch = 20, cex = 0.7, pt.cex = 2)
dev.off()

# calculate AUCell
set.seed(seed)
aucellRankings <- AUCell_buildRankings(dgem, nCores = 20, plotStats = F)
regulonAUC <- AUCell_calcAUC(regs, aucellRankings, 
		aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 20)
auc <- getAUC(regulonAUC)

# build new loom
fileName <- paste0(gset, "_MMregulons.loom")
newloom <- build_loom(file.name = fileName,
	title = gset,
	genome = "human",
	dgem = dgem,
	default.embedding = emb
)
add_scenic_regulons_auc_matrix(loom = newloom, regulons.AUC = auc)
add_scenic_regulons(loom = newloom, regulons = regs)

add_col_attr(newloom, key = "celltype", value = meta[, "celltype"], as.annotation = T)
add_col_attr(newloom, key = "sample", value = meta[, "sample"], as.annotation = T)
add_col_attr(newloom, key = "treatment_group", value = meta[, "treatment_group"], as.annotation = T)
add_col_attr(newloom, key = "cohort", value = meta[, "cohort"], as.annotation = T)
add_col_attr(newloom, key = "no_of_reads", value = meta[, "no_of_reads"], as.metric = T)

## add hierarchy
tree <- c("Cancer Compendium", "", "")
for (treen in seq_along(tree)) {
	dtype <- hdf5r::guess_dtype(tree[treen])
	dtype <- SCopeLoomR:::hdf5_utf8_encode(value = tree[treen], dtype = dtype)
	newloom$create_attr(attr_name = paste0("SCopeTreeL", treen), robj = tree[treen], dtype = dtype,
		space = H5S$new(type = "scalar", dims = NULL, maxdims = NULL))
}
list.attributes(newloom)

finalize(loom = newloom)

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
# [1] SCopeLoomR_0.5.0  rjson_0.2.20      doRNG_1.7.1       rngtools_1.4
# [5] pkgmaker_0.27     registry_0.5-1    foreach_1.4.7     hdf5r_1.2.0
# [9] AUCell_1.6.1      data.table_1.12.4
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2                  lattice_0.20-38
# [3] zeallot_0.1.0               digest_0.6.21
# [5] mime_0.7                    plyr_1.8.4
# [7] R6_2.4.0                    GenomeInfoDb_1.20.0
# [9] backports_1.1.5             stats4_3.6.1
#[11] RSQLite_2.1.2               pillar_1.4.2
#[13] zlibbioc_1.30.0             rlang_0.4.0
#[15] annotate_1.62.0             blob_1.2.0
#[17] S4Vectors_0.22.1            R.utils_2.9.0
#[19] R.oo_1.22.0                 Matrix_1.2-17
#[21] BiocParallel_1.18.1         stringr_1.4.0
#[23] RCurl_1.95-4.12             bit_1.1-14
#[25] shiny_1.3.2                 DelayedArray_0.10.0
#[27] compiler_3.6.1              httpuv_1.5.2
#[29] base64enc_0.1-3             pkgconfig_2.0.3
#[31] BiocGenerics_0.30.0         htmltools_0.3.6
#[33] SummarizedExperiment_1.14.1 tibble_2.1.3
#[35] GenomeInfoDbData_1.2.1      IRanges_2.18.3
#[37] codetools_0.2-16            matrixStats_0.55.0
#[39] doMC_1.3.6                  XML_3.98-1.20
#[41] crayon_1.3.4                withr_2.1.2
#[43] later_0.8.0                 bitops_1.0-6
#[45] R.methodsS3_1.7.1           grid_3.6.1
#[47] xtable_1.8-4                GSEABase_1.46.0
#[49] DBI_1.0.0                   magrittr_1.5
#[51] bibtex_0.4.2                graph_1.62.0
#[53] stringi_1.4.3               reshape2_1.4.3
#[55] XVector_0.24.0              promises_1.0.1
#[57] doParallel_1.0.15           vctrs_0.2.0
#[59] iterators_1.0.12            tools_3.6.1
#[61] bit64_0.9-7                 Biobase_2.44.0
#[63] parallel_3.6.1              AnnotationDbi_1.46.1
#[65] GenomicRanges_1.36.1        memoise_1.1.0
#}
