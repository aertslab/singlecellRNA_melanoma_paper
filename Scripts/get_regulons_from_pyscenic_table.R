# script to get regulons from pyscenic ctx output table
# and top motifs for each regulon (max NES)

#module load R/3.6.1-foss-2018a-X11-20180604 
#module load GCC/6.4.0-2.28 
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28 
#PATH=$PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/bin/ 
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VSC_HOME/progs/hdf5-1.10.4/bin/lib/ 
options(stringsAsFactors=FALSE) 
.libPaths("/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/3.6") 
staging <- "/ddn1/vol1/staging/leuven/stg_00002/lcb" 
############################################################################################# 
library(data.table) 
library(parallel)

sample <- "5.MM074_SOX"
scenicdir <- file.path(staging, "dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices", sample, "pySCENIC_100x/out/cistarget")
setwd(scenicdir)

# function to get all target genes for a specific TF from reg table  
getRegulon <- function(reg, tf, d) {  
   #d: "neg"/"pos" 
   digitpattern <- "[[:digit:]]+\\.[[:digit:]e-]+"
   alltargets <- lapply(1:reg[TF == tf & dir == d, .N], function(i) {  
        targets <- strsplit(  
         gsub(digitpattern, "",  
          gsub(" ", "",  
           gsub("'", "",   
            gsub("]", "",  
             gsub("[", "",   
              gsub("[()]", "", reg[TF == tf & dir == d, TargetGenes][i]),   
             fixed = TRUE),   
            fixed = TRUE),  
           )  
          )  
         )  
        , ",")[[1]]  
   targets[sapply(targets, function(p) {p != ""})]  
   })  
   Reduce(union, alltargets)  
}

mclapply(1:100, function(runno) {
	motiffile <- paste0("run_", runno, "/run_", runno, "__reg_mtf.csv")
        tab <- fread(motiffile, skip = 3, sep = ",")
        colnames(tab) <- c("TF", "MotifID", "AUC", "Annotation", "Context", "MotifSimilarityQvalue",
			                  "NES", "OrthologousIdentity", "RankAtMax", "TargetGenes")

	tab[, dir := ifelse(grepl("activating", Context), "pos", "neg")]

	regulons <- list()
	d <- "neg"
	regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
		getRegulon(tab, tf, d)})  
	d <- "pos"
	regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
		getRegulon(tab, tf, d)})
	
	print("final regulons:")
	print(sapply(regulons, length))
	
	saveRDS(regulons, file = paste0("run_", runno, "/reg_mtf.rds.gz"), compress = "gzip")

        tab <- tab[, .(TF, MotifID, NES)]
        tab <- tab[tab[, .I[NES == max(NES)], by = TF]$V1]

        saveRDS(tab, file = paste0("run_", runno, "/top_mtf_per_TF.rds.gz"), compress = "gzip")
}, mc.cores = 30)

#combine top motifs into one table
motifs <- rbindlist(mclapply(1:100, function(runno) {
	tab <- readRDS(file.path("100_scenic_runs", paste0("run_", runno, "/top_mtf_per_TF.rds.gz")))
	tab[, run := runno]
	tab
}, mc.cores = 30))
saveRDS(motifs, file = "top_mtf_per_TF_100runs.rds.gz", compress = "gzip")

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
#[1] parallel  stats     graphics  grDevices utils     datasets  methods  
#[8] base     
#
#other attached packages:
#[1] data.table_1.12.8
#
#loaded via a namespace (and not attached):
#[1] compiler_3.6.1 tools_3.6.1   
#}
