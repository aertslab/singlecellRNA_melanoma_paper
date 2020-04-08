# task: UMAP on AUCell values of filtered regulons
# Fig3c & loom files

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
library(reticulate)
library(irlba)
library(Rtsne)
seed <- 42

# set UMAP parameters
umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
reticulate::py_set_seed(seed)

max.dim = 2L
metric = "correlation"
n_neighbours = 100
min_dist = 0.75
pcs = 5 

umap <- umap_import$UMAP(
         n_neighbors = as.integer(x = n_neighbours),
         n_components = as.integer(x = max.dim),
         metric = metric, 
         min_dist = min_dist)

samples <- c("12.ten_lines_BL", "5.MM074_SOX", "9.MM057_SOX_wo_TL", "10.M0087_SOX_wo_TL")
for (sample in samples) {
	print(sample)
        scenicdir <- file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets")
	wd <- file.path(scenicdir, sample)
	setwd(wd)

	# get AUCell values of filtered regulons
	auc <- readRDS("all_pos_regulons_recurrent_100x_jw_new_trk_aucell01.rds.gz")

	# pca
	npcs = 20
	npcs <- min(npcs, nrow(x = auc) - 1)
	pca.results <- irlba(A = t(x = auc), nv = npcs)
	ce <- pca.results$u %*% diag(pca.results$d) 
	rownames(ce) <- colnames(auc) 

	# umap
	set.seed(seed)
	umap_output <- umap$fit_transform(as.matrix(x = ce[,1:pcs]))
	rownames(x = umap_output) <- rownames(ce) 

	write.table(umap_output, file = "auc_100runs_umap1_umap2.tsv", quote = F, sep = "\t")
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
#[1] Rtsne_0.15        irlba_2.3.3       Matrix_1.2-17     reticulate_1.13  
#[5] data.table_1.12.8
#
#loaded via a namespace (and not attached):
#[1] compiler_3.6.1  Rcpp_1.0.3      grid_3.6.1      jsonlite_1.6.1 
#[5] lattice_0.20-38
#}
