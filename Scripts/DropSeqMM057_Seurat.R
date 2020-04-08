# make count matrices with cells that have >=1000 genes expresses
# combined runs [[NextSeq20170321_DropSeq]] and [[NextSeq20170406_DropSeq]]

library(Seurat)	
### functions {
# function to merge treatments of same cell line
mergetreats <- function(m1, m2) {
  #m1, m2: count matrices
  m <- merge(m1, m2, by = "row.names", all = T)
  m[is.na(m)] <- 0
  rownames(m) <- m[,1]
  m <- m[,-1]
  print(dim(m))
  m
}

# function to filter out cells with <=1000 genes expressed
cellswith1kgenes <- function(countmatrix) {
    num_genes <- colSums(apply(countmatrix, 2, function(x) (x!=0))*1)
    print(paste("keeping ", as.vector(table(num_genes>=1000)[2]), " cells", sep =""))
    subset(countmatrix, select = num_genes>=1000)
}

# function to make list object containing read/gene count etc.
cellstats <- function(m, statsfile, treat) {
  #m: count matrix
  #statsfile: *_genecounts_vs_cells.txt
  #treat: treatment (e.g. "BL", "siSOX10")
  cells <- colnames(m)
  stats <- read.table(statsfile)
  
  #add column with treatment information
  cellstats <- cbind(stats[cells,], treatment= c(rep(treat, length(cells))))
}
########################################################################################################}

treatment <- "MM057"
samplename1 <- "DropSeq_MM057_BL"
treatment1 <- "BL"

datafile <- "MM057_1k.RData"

## read count matrix treat1
treat1 <- read.table(paste(treatment, "/", samplename1, "_clean.dge.txt.gz", sep = ""),
	header = T, row.names = 1)

## filter out cells with <=1k genes
treat1 <- cellswith1kgenes(treat1)
#774 cells

## append treatment name to cell barcode
colnames(treat1) <- as.vector(sapply(colnames(treat1), function(x) x <- paste(x, "_", treatment1, sep = "")))

## stats (#genes & treatment)
stats1 <- as.data.frame(colSums(apply(treat1, 2, function(x) (x!=0))*1))
stats1 <- cbind(stats1, rep(treatment1, dim(stats1)[1]))
colnames(stats1) <- c("genes", "treatment")
rownames(stats1) <- as.vector(sapply(rownames(stats1), function(x) x <- paste(x, "_", treatment1, sep = "")))

## join with old data run NextSeq_20170222_Drop-seq
treatment2 <- "NTC"
treat2 <- cellswith1kgenes(read.table("/home/katina/uni/data/DropSeqMM057/DropSeq_MM057_NTC_clean.dge.txt.gz", header = T, row.names = 1))
colnames(treat2) <- as.vector(sapply(colnames(treat2), function(x) x <- paste(x, "_", treatment2, sep = "")))
stats2 <- as.data.frame(colSums(apply(treat2, 2, function(x) (x!=0))*1))
stats2 <- cbind(stats2, rep(treatment2, dim(stats2)[1]))
colnames(stats2) <- c("genes", "treatment")
rownames(stats2) <- as.vector(sapply(rownames(stats2), function(x) x <- paste(x, "_", treatment2, sep = "")))

treatment3 <- "SOX10"
treat3 <- cellswith1kgenes(read.table("/home/katina/uni/data/DropSeqMM057/DropSeq_MM057_SOX10_clean.dge.txt.gz", header = T, row.names = 1))
colnames(treat3) <- as.vector(sapply(colnames(treat3), function(x) x <- paste(x, "_", treatment3, sep = "")))
stats3 <- as.data.frame(colSums(apply(treat3, 2, function(x) (x!=0))*1))
stats3 <- cbind(stats3, rep(treatment3, dim(stats3)[1]))
colnames(stats3) <- c("genes", "treatment")
rownames(stats3) <- as.vector(sapply(rownames(stats3), function(x) x <- paste(x, "_", treatment3, sep = "")))

## merge treatments
bothtreats <- mergetreats(treat1, treat2)
alltreats <- list(exprMatrix = mergetreats(bothtreats, treat3),
	stats = rbind(stats1, stats2, stats3))
	
save(alltreats, file = datafile)

# Initialize the Seurat object with the raw (non-normalized data)
# filter simlarly to 10x data: [[HiSeq4000_20170504CountMatrices.R]]
MM057 <- new("seurat", raw.data = alltreats$exprMatrix)
MM057 <- Setup(MM057, 
		"MM057", 						#project name
		meta.data = alltreats$stats[,"treatment",drop = F],	#nUMI & nGene calculated automatically
		min.genes = 999,
		total.expr = 1e4)
		
#add % of mitochondrial genes
mito.genes <- grep("^MT-", rownames(MM057@data), value = T)
percent.mito <- colSums(as.matrix(expm1(MM057@data[mito.genes, ])))*100/colSums(as.matrix(expm1(MM057@data)))
MM057 <- AddMetaData(MM057, percent.mito, "percent.mito")

png("MM057/seurat_MM057_vlnplot.png", width = 1200, height = 600)
VlnPlot(MM057, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "treatment", cols.use = c("cadetblue4", "lightblue", "orange"))
dev.off()

# choose cutoffs
cut_umi <- 50000
png("MM057/seurat_MM057_nUMI-cutoff.png", width = 800, height = 600)
par(cex = 1.5)
t <- MM057@data.info$treatment
plot(MM057@data.info$nUMI, 
	MM057@data.info$nGene, 
	col=c("cadetblue4", "lightblue", "orange")[as.numeric(t)],
	pch = 18, cex = 1,
	xlab = "nUMI", ylab = "nGene")
legend("topleft",legend = unique(t),col=c("cadetblue4", "lightblue", "orange")[unique(t)],pch=18, title = "MM057")
abline(v=cut_umi)
b <- dim(MM057@data.info)[1]
a <- length(which(MM057@data.info$nUMI <cut_umi))
text(30000, 3000, labels = print(paste(a, " of ", b, " cells kept", sep = "")))
dev.off()

cut_mito <- 10
png("MM057/seurat_MM057_mito-cutoff.png", width = 800, height = 600)
par(cex = 1.5)
t <- MM057@data.info$treatment
plot(MM057@data.info$nUMI, 
	MM057@data.info$percent.mito, 
	col=c("cadetblue4", "lightblue", "orange")[as.numeric(t)],
	pch = 18, cex = 1,
	xlab = "nUMI", ylab = "percent.mito")
legend("topright",legend = unique(t),col=c("cadetblue4", "lightblue", "orange")[unique(t)],pch=18, title = "MM057")
abline(h=cut_mito)
b <- dim(MM057@data.info)[1]
a <- length(which(MM057@data.info$percent.mito <cut_mito))
text(30000, 13, labels = print(paste(a, " of ", b, " cells kept", sep = "")))
dev.off()

# filter
MM057 <- SubsetData(MM057, subset.name = "nUMI", accept.high = cut_umi)
MM057 <- SubsetData(MM057, subset.name = "percent.mito", accept.high = cut_mito)
#remove genes with 0 expression across this set of cells
MM057@data <- MM057@data[which(rowSums(as.matrix(MM057@data))!=0),]
dim(MM057@data)
#[1] 19188  2278

#subset raw data, too; save
MM057@raw.data <- MM057@raw.data[rownames(MM057@data), colnames(MM057@data)]

MM057 <- list(exprMatrix = MM057@raw.data,
                  info = MM057@data.info)
save(MM057, file="MM057_seurat.RData")

#how many genes are expressed in less than 1% of all cells?
num_cells <- rowSums(apply(MM057$exprMatrix, 2, function(x) (x!=0))*1)
summary(num_cells)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.0    12.0    85.0   296.1   383.0  2278.0 
length(which(num_cells < round(dim(MM057$exprMatrix)[2]*0.01)))
#[1] 6409

summary(MM057$info)
#     nGene           nUMI       treatment    orig.ident    percent.mito   
# Min.   :1000   Min.   : 1450   BL   : 773   MM057:2278   Min.   :0.7027  
# 1st Qu.:1621   1st Qu.: 3261   NTC  : 379                1st Qu.:3.4855  
# Median :2245   Median : 5175   SOX10:1126                Median :4.3711  
# Mean   :2498   Mean   : 6789                             Mean   :4.4853  
# 3rd Qu.:3121   3rd Qu.: 8462                             3rd Qu.:5.3297  
# Max.   :7636   Max.   :47897                             Max.   :9.9880  

sessionInfo() #{
#R version 3.4.4 (2018-03-15)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Debian GNU/Linux 8 (jessie)
#
#Matrix products: default
#BLAS: /usr/lib/openblas-base/libblas.so.3
#LAPACK: /usr/lib/libopenblasp-r0.2.12.so
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
#[1] Seurat_2.3.4  Matrix_1.2-14 cowplot_0.9.2 ggplot2_3.2.1
#
#loaded via a namespace (and not attached):
#  [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15       
#  [4] modeltools_0.2-22   ggridges_0.5.1      mclust_5.4.5       
#  [7] htmlTable_1.13.2    base64enc_0.1-3     rstudioapi_0.10    
# [10] proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-14     
# [13] bit64_0.9-7         codetools_0.2-16    splines_3.4.4      
# [16] R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-2  
# [19] knitr_1.26          zeallot_0.1.0       Formula_1.2-3      
# [22] jsonlite_1.6        ica_1.0-2           cluster_2.1.0      
# [25] kernlab_0.9-27      png_0.1-7           R.oo_1.23.0        
# [28] compiler_3.4.4      httr_1.4.1          backports_1.1.5    
# [31] assertthat_0.2.0    lazyeval_0.2.2      lars_1.2           
# [34] acepack_1.4.1       htmltools_0.4.0     tools_3.4.4        
# [37] igraph_1.2.2        gtable_0.3.0        glue_1.3.1         
# [40] RANN_2.6            reshape2_1.4.3      dplyr_0.8.3        
# [43] Rcpp_1.0.3          vctrs_0.2.0         gdata_2.18.0       
# [46] ape_5.1             nlme_3.1-142        iterators_1.0.12   
# [49] fpc_2.2-3           gbRd_0.4-11         lmtest_0.9-37      
# [52] xfun_0.3            stringr_1.4.0       lifecycle_0.1.0    
# [55] irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8     
# [58] MASS_7.3-51.4       zoo_1.8-6           scales_1.0.0       
# [61] doSNOW_1.0.18       parallel_3.4.4      RColorBrewer_1.1-2 
# [64] reticulate_1.10     pbapply_1.4-2       gridExtra_2.3      
# [67] rpart_4.1-15        segmented_1.0-0     latticeExtra_0.6-28
# [70] stringi_1.2.4       foreach_1.4.7       checkmate_1.9.4    
# [73] caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.9-0       
# [76] SDMTools_1.1-221.1  rlang_0.4.2         pkgconfig_2.0.2    
# [79] dtw_1.21-3          prabclus_2.3-1      bitops_1.0-6       
# [82] lattice_0.20-38     ROCR_1.0-7          purrr_0.3.3        
# [85] htmlwidgets_1.5.1   bit_1.1-14          tidyselect_0.2.5   
# [88] plyr_1.8.4          magrittr_1.5        R6_2.2.2           
# [91] snow_0.4-3          gplots_3.0.1        Hmisc_4.3-0        
# [94] pillar_1.4.2        foreign_0.8-72      withr_2.1.2        
# [97] fitdistrplus_1.0-14 mixtools_1.1.0      survival_3.1-7     
#[100] nnet_7.3-12         tsne_0.1-3          tibble_2.1.3       
#[103] crayon_1.3.4        hdf5r_1.3.0         KernSmooth_2.23-15 
#[106] grid_3.4.4          data.table_1.12.6   metap_1.1          
#[109] digest_0.6.23       diptest_0.75-7      tidyr_1.0.0        
#[112] R.utils_2.9.0       stats4_3.4.4        munsell_0.5.0      
#}
