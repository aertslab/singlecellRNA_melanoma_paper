# module load R/3.3.2-foss-2014a-noX
# module load cairo/1.14.6-foss-2014a

.libPaths('/user/leuven/303/vsc30324/lcb/zkalender/R/x86_64-pc-linux-gnu-library/3.3')
library(Seurat)
# setwd("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/Seurat_analysis_on_40k_set")
# load("MMlines_scRNA_seq_40k_seurat_obj.RData")
setwd("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data")
load("all_data_combined.RData")
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MMlines_scRNA_seq_40k_seurat_obj@data), value = TRUE)
percent.mito <- Matrix::colSums(MMlines_scRNA_seq_40k_seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(MMlines_scRNA_seq_40k_seurat_obj@raw.data)
# stash QC stats
MMlines_scRNA_seq_40k_seurat_obj <- AddMetaData(object = MMlines_scRNA_seq_40k_seurat_obj, metadata = percent.mito, col.name = "percent.mito")
# VLN plot
setwd("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/10.analysis_with_all_data")
png(file="20.Plots/MMlines_scRNA_seq_40k_VLN_plot.png", width=2400, height=1400,res=96)
VlnPlot(object = MMlines_scRNA_seq_40k_seurat_obj, group.by="CE_ID", features.plot = c("nGene", "nUMI", "percent.mito"), nCol=1, x.lab.rot=T, point.size.use=0.7)
dev.off()
# Gene Plots UMI vs nGenes
png(file="20.Plots/MMlines_scRNA_seq_40k_UMI_vs_nGene_plot.png", width=800, height=800,res=96)
GenePlot(object = MMlines_scRNA_seq_40k_seurat_obj, gene1 = "nUMI", gene2 = "nGene",cex.use=0.7,col.use="black")
points(MMlines_scRNA_seq_40k_seurat_obj@meta.data[MMlines_scRNA_seq_40k_seurat_obj@meta.data$percent.mito>=0.2,"nUMI"],MMlines_scRNA_seq_40k_seurat_obj@meta.data[MMlines_scRNA_seq_40k_seurat_obj@meta.data$percent.mito>=0.2,"nGene"],col="red",pch=20,cex=0.7)
dev.off()

# Gene Plots UMI vs percent mitochondria
png(file="20.Plots/MMlines_scRNA_seq_40k_UMI_vs_Mito_plot.png", width=800, height=800,res=96)
GenePlot(object = MMlines_scRNA_seq_40k_seurat_obj, gene1 = "nUMI", gene2 = "percent.mito",cex.use=0.7,col.use="black")
dev.off()

# UMI histogram
png(file="20.Plots/UMI_histogram.png",width=800, height=800,res=96)
hist(MMlines_scRNA_seq_40k_seurat_obj@meta.data$nUMI, breaks=100,main="UMI counts",xlab="")
dev.off()

# Gene histogram
png(file="20.Plots/nGene_histogram.png",width=800, height=800,res=96)
hist(MMlines_scRNA_seq_40k_seurat_obj@meta.data$nGene, breaks=100,main="Gene counts",xlab="")
dev.off()

# png("variable_genes_wo_CPM.png",width=800, height=800,res=150)
# mixed10x_seurat_obj <- FindVariableGenes(object = MMlines_scRNA_seq_40k_seurat_obj, do.plot = T)
# dev.off()

# Filter cells
# dim(MMlines_scRNA_seq_40k_seurat_obj@meta.data[MMlines_scRNA_seq_40k_seurat_obj@meta.data$percent.mito>=0.2,])
# 530
MMlines_scRNA_seq_40k_seurat_obj <- FilterCells(object = MMlines_scRNA_seq_40k_seurat_obj, subset.names = c("percent.mito"),
    low.thresholds = c(-Inf), high.thresholds = c(0.2))

MMlines_scRNA_seq_40k_seurat_obj@ident<-factor(MMlines_scRNA_seq_40k_seurat_obj@meta.data$CE_ID)
names(MMlines_scRNA_seq_40k_seurat_obj@ident)<-rownames(MMlines_scRNA_seq_40k_seurat_obj@meta.data)

raw_counts<-MMlines_scRNA_seq_40k_seurat_obj@raw.data
raw_counts<-raw_counts[,rownames(MMlines_scRNA_seq_40k_seurat_obj@meta.data)]
raw_counts_matrix<-as.matrix(raw_counts)
write.table(file="10.TextData/raw_counts.tsv",raw_counts_matrix, sep="\t",quote=F)
save(file="00.RData/seurat_obj_with_raw_counts.RData",MMlines_scRNA_seq_40k_seurat_obj)
#dim(raw_counts_matrix)
#[1] 32738 43582



filt_counts_cpm<-edgeR::cpm.default(raw_counts,log = F)
save(file="00.RData/CPM_counts.RData",filt_counts_cpm)

# Convert to Seurat object
full_meta_data<-MMlines_scRNA_seq_40k_seurat_obj@meta.data
full_meta_data[sapply(full_meta_data, is.character)] <- lapply(full_meta_data[sapply(full_meta_data, is.character)],
                                       as.factor)
MMlines_scRNA_seq_40k_CPM_seurat_obj <- CreateSeuratObject(raw.data = filt_counts_cpm, meta.data = full_meta_data)
MMlines_scRNA_seq_40k_CPM_seurat_obj@ident<-MMlines_scRNA_seq_40k_CPM_seurat_obj@meta.data$CE_ID
names(MMlines_scRNA_seq_40k_CPM_seurat_obj@ident)<-rownames(MMlines_scRNA_seq_40k_CPM_seurat_obj@meta.data)

save(file="00.RData/CPM_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_seurat_obj)

## Create subsets for further analysis

# if creating new subsets: load("~/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/10.analysis_with_all_data/00.RData/CPM_seurat_obj.RData")

####  Write raw counts to text files for bokeh
# SOX perturbations across three lines
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM057_SOX10_24h", "MM087_SOX10_24h", "MM087_SOX10_48h", "MM074_SOX10_48h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM074_SOX10_72h", "MM087_SOX10_72h", "MM057_NTC", "MM074_NTC", "MM087_NTC", "MM087_BL", "MM057_BL", "MM074_BL", "MM057_TL", "MM087_TL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/1.three_MM_lines_SOX/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# SOX & TGF perturbations across three lines
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM057_SOX10_24h", "MM087_SOX10_24h", "MM087_SOX10_48h", "MM074_SOX10_48h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM074_SOX10_72h", "MM087_SOX10_72h", "MM057_NTC", "MM074_NTC", "MM087_NTC", "MM087_BL", "MM057_BL", "MM074_BL", "MM057_TL", "MM087_TL" , "MM057_TGF_TNF", "MM057_TGF_only", "MM074_TGF_TNF", "MM074_TGF_only", "MM087_TGF_TNF", "MM087_TGF_only"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/2.three_MM_lines_SOX_and_TGF/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM057 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL", "MM057_TL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/3.MM057_SOX/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM074 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM074_SOX10_48h", "MM074_SOX10_72h", "MM074_NTC", "MM074_BL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM087 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL", "MM087_TL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/7.M0087_SOX/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM057 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL", "MM057_TL", "MM057_TGF_TNF", "MM057_TGF_only" ))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/4.MM057_SOX_TGF/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM074 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM074_SOX10_48h", "MM074_SOX10_72h", "MM074_NTC", "MM074_BL","MM074_TGF_TNF", "MM074_TGF_only"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/6.MM074_SOX_TGF/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM087 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL", "MM087_TL","MM087_TGF_TNF", "MM087_TGF_only"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/8.M0087_SOX_TGF/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM057 SOX10  without timeline
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)

# MM087 SOX10 without timeline
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL"))
seurat_raw_counts<-seurat_obj@raw.data
seurat_raw_counts<-seurat_raw_counts[,rownames(seurat_obj@meta.data)]
seurat_raw_counts<-as.matrix(seurat_raw_counts)
write.table(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/10.TextData/raw_counts.tsv",seurat_raw_counts, sep="\t",quote=F)
rm(seurat_obj)


# SOX perturbations across three lines
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM057_SOX10_24h", "MM087_SOX10_24h", "MM087_SOX10_48h", "MM074_SOX10_48h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM074_SOX10_72h", "MM087_SOX10_72h", "MM057_NTC", "MM074_NTC", "MM087_NTC", "MM087_BL", "MM057_BL", "MM074_BL", "MM057_TL", "MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/1.three_MM_lines_SOX/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# SOX & TGF perturbations across three lines
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM057_SOX10_24h", "MM087_SOX10_24h", "MM087_SOX10_48h", "MM074_SOX10_48h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM074_SOX10_72h", "MM087_SOX10_72h", "MM057_NTC", "MM074_NTC", "MM087_NTC", "MM087_BL", "MM057_BL", "MM074_BL", "MM057_TL", "MM087_TL" , "MM057_TGF_TNF", "MM057_TGF_only", "MM074_TGF_TNF", "MM074_TGF_only", "MM087_TGF_TNF", "MM087_TGF_only"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/2.three_MM_lines_SOX_and_TGF/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM057 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL", "MM057_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/3.MM057_SOX/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM074 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM074_SOX10_48h", "MM074_SOX10_72h", "MM074_NTC", "MM074_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL", "MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/7.M0087_SOX/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM057 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL", "MM057_TL", "MM057_TGF_TNF", "MM057_TGF_only" ))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/4.MM057_SOX_TGF/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM074 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_SOX10_24h", "MM074_SOX10_48h", "MM074_SOX10_72h", "MM074_NTC", "MM074_BL","MM074_TGF_TNF", "MM074_TGF_only"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/6.MM074_SOX_TGF/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 & TGF
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL", "MM087_TL","MM087_TGF_TNF", "MM087_TGF_only"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/8.M0087_SOX_TGF/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM057 SOX10  without timeline
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM057_SOX10_24h", "MM057_SOX10_48h", "MM057_SOX10_72h", "MM057_NTC", "MM057_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 without timeline
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 without timeline
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("A375_BL", "MM011_BL", "MM029_BL", "MM087_BL","MM099_BL","MM031_BL", "MM047_BL", "MM057_BL", "MM074_BL", "MM001_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/12.ten_lines_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM057 BL & NTC
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM057_NTC", "MM057_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/15.MM057_BL_NTC/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM074 BL & NTC
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_NTC", "MM074_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/16.MM074_BL_NTC/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 BL & NTC
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_NTC", "MM087_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/17.MM087_BL_NTC/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 TL_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/18.MM087_TL_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 w/o BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h", "MM087_NTC", "MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/19.MM087_SOX_wo_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 w/o BL & 72h
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_NTC", "MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/20.MM087_SOX_wo_BL_and_72h/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 w/o NTC & 72h
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_BL", "MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/21.MM087_SOX_wo_NTC_and_72h/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 NTC_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_NTC"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/22.MM087_NTC_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 24h_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_24h"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/23.MM087_SOX_24h_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 BL_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/24.MM087_BL_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 NTC and TL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_NTC","MM087_TL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/25.MM087_NTC_and_TL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

<<<<<<< HEAD
# MM087 w/o BL and NTC
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_TL","MM087_SOX10_24h", "MM087_SOX10_48h", "MM087_SOX10_72h"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/26.MM087_wo_NTC_and_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 48h_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_48h"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/27.MM087_SOX_48h_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087 SOX10 72h_only
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_SOX10_72h"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/28.MM087_SOX_72h_only/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

=======
# A375_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("A375_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/30.A375_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM001_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM001_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/31.MM001_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM011_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM011_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/32.MM011_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM029_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM029_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/33.MM029_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM031_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM031_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/34.MM031_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM047_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM047_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/35.MM047_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM057_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM057_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/36.MM057_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM074_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM074_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/37.MM074_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# MM087_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM087_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/38.MM087_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)

# 39.MM099_BL
seurat_obj<-SubsetData(MMlines_scRNA_seq_40k_CPM_seurat_obj,ident.use=c("MM099_BL"))
seurat_obj@raw.data<-seurat_obj@raw.data[,colnames(seurat_obj@data)]
save(file="/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/39.MM099_BL/00.RData/CPM_seurat_obj.RData",seurat_obj)
rm(seurat_obj)
>>>>>>> ce80842d38cde9233aad3075bcc8840072817574
# Continue from here on analysis on whole datasets

# Filter genes
# The gene should be expressed at least 1% of samples == 436 samples

dim(MMlines_scRNA_seq_40k_CPM_seurat_obj@raw.data)
# 32738 43582
filter_genes <- apply(raw_counts, 1, function(x) length(x[x > 1]) >= 435)
table(filter_genes)
#filter_genes
#FALSE  TRUE
#22971  9767

#raw_counts_filt<-raw_counts[filter_genes,]

# Create a new seurat object from filtered genes
misclassified_samples<-subset(misclassified_samples, Cell_Line=="NA")
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj <- CreateSeuratObject(raw.data = filt_counts_cpm[filter_genes,! colnames(filt_counts_cpm) %in% rownames(misclassified_samples)],
  meta.data = full_meta_data[ ! rownames(full_meta_data) %in% rownames(misclassified_samples),])

MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$CE_ID<-factor(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$CE_ID)
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$Experiment<-factor(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$Experiment)
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$Cell_Line<-factor(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$Cell_Line)
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$CE_ID
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident<-factor(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident)

cpm_data<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@raw.data
save(file="00.RData/gene_filt_CPM_matrix.RData",cpm_data)
print("Exporting CPM matrix & meta data ")
write.table(file="10.TextData/CPM_matrix.tsv",cpm_data, sep="\t",quote=F)
t_cpm_data<-t(cpm_data)
write.table(file="10.TextData/CPM_transposed_matrix.tsv",t_cpm_data, sep="\t",quote=F)
meta_data<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data
save(file="00.RData/meta_data.RData",meta_data)
write.table(file="10.TextData/meta_data.tsv",meta_data, sep="\t",quote=F)
log_CPM<-log2(cpm_data+1)
write.table(file="10.TextData/log_CPM_matrix.tsv",log_CPM, sep="\t",quote=F)
save(file="00.RData/log_CPM_matrix.RData",log_CPM)
t_log_CPM<-t(log_CPM)
write.table(file="10.TextData/log_CPM_matrix_transposed.tsv",t_log_CPM, sep="\t",quote=F)

# SCENIC prep
library(RcisTarget.hg19.motifDatabases.20k)
org<-"hg19"
data(hg19_500bpUpstream_motifRanking)
data(hg19_direct_motifAnnotation)
allTFs <- hg19_direct_motifAnnotation$allTFs
inputTFs <- allTFs[allTFs %in% rownames(cpm_data)]
save(inputTFs, file="30.SCENIC/int/1.2_inputTFs.RData")

system("ln -s ../../00.RData/log_CPM_matrix.RData 30.SCENIC/int/1.1_exprMatrix_filtered.RData")


# Normalize & Scale
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj<-NormalizeData(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj<-ScaleData(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)
save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)


# Joy Plots
png(file="20.Plots/joy_plot_MITF_SOX10_SOX9.png", width=2400, height=1200,res=96)
features.plot <- c("MITF", "SOX10", "SOX9")
JoyPlot(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, features.plot = features.plot, group.by="CE_ID", nCol=3)
dev.off()

# Scale (Regress out UMI counts)
# MMlines_scRNA_seq_40k_CPM_seurat_obj <- ScaleData(object = MMlines_scRNA_seq_40k_CPM_seurat_obj, vars.to.regress = c("nUMI"))

pdf("20.Plots/variable_genes_with_CPM.pdf")
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj <- FindVariableGenes(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, do.plot = T)
dev.off()

### PCA
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj <- RunPCA(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, pc.genes = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@var.genes,pcs.compute = 20)

pca_cords<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,1:2]
write.table(file="10.TextData/seurat_pc1_pc2.tsv", pca_cords,sep="\t",quote=F)

png(file="20.Plots/PCA_elbow_plot.png",width=800, height=800,res=96)
PCElbowPlot(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,num.pc = 20)
dev.off()

png(file="20.Plots/PC1_heatmap.png",width=800, height=800,res=96)
PCHeatmap(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, col.use = colorRampPalette(c("red", "white", "blue")), cexRow=0.7, dendrogram="none", key=F)
dev.off()

png(file="20.Plots/PC2_heatmap.png",width=800, height=800,res=96)
PCHeatmap(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, col.use = colorRampPalette(c("red", "white", "blue")), cexRow=0.7, dendrogram="none", key=F)
dev.off()

png(file="20.Plots/PCA_plot_PC1_vs_PC2.png",width=800, height=800,res=96)
PCAPlot(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, dim.1 = 1, dim.2 = 2, group.by="CE_ID",no.legend=T)
dev.off()

save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)

# Color scale for nGenes
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
colorNgenes<-setNames(adjustcolor(colorPal(10),alpha=.8)[as.numeric(cut(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$nUMI, breaks=10,right=F,include.lowest=T))], rownames(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data))

# PCA with color nGenes
png(file="20.Plots/PCA_plot_PC1_vs_PC2_with_nGenes.png",width=800, height=800,res=96)
plot(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,1],
  MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,2],
  pch=20,
  col= colorNgenes[rownames(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings)],
  xlab="PC1",
  ylab="PC2",
  cex=1 )
dev.off()

png(file="20.Plots/PCA_plot_PC1_vs_PC2_with_CE_ID.png",width=800, height=800,res=96)
plot(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,1],
  MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,2],
  pch=20,
  col= as.factor(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@meta.data$CE_ID),
  xlab="PC1",
  ylab="PC2",
  cex=1 )
dev.off()

# save PCA coordinates
pca_cords<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$pca@cell.embeddings[,1:2]
write.table(file="10.TextData/seurat_pc1_pc2.tsv", pca_cords,sep="\t",quote=F)


# tSNE
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj <- RunTSNE(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj, do.fast = TRUE,check_duplicates = FALSE)

png(file="20.Plots/TSNE_plot.png",width=1600, height=800,res=96)
TSNEPlot(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,group.by="CE_ID")
dev.off()

MMlines_scRNA_seq_40k_CPM_seurat_obj<-RunDiffusion(MMlines_scRNA_seq_40k_CPM_seurat_obj,max.dim = 3)


tsne_palette<-c("#37195d",
"#cd3f96",
"#ce9a7f",
"#509987",
"#c5d4bc",
"#2e5157",
"#b55039",
"#8cac3c",
"#53713b",
"#6ce1c0",
"#9ed28f",
"#41499a",
"#c4a239",
"#d43869",
"#e04431",
"#67abd7",
"#bfc6dd",
"#583329",
"#713fde",
"#9c3fc5",
"#772327",
"#7669e0",
"#449342",
"#84652f",
"#67d4df",
"#d8e046",
"#78db48",
"#65527c",
"#d392d5",
"#df758f",
"#d173d6",
"#6380d8",
"#893d8e",
"#db46d5",
"#323c21",
"#7f807d",
"#9c5d70",
"#60df89",
"#2c1f36",
"#db8439",
"#d6a2b6",
"#42719b",
"#d9ce82",
"#6b204b",
"#a59bd4")

png(file="20.Plots/TSNE_plot_with_custom_colors.png",width=1024, height=768,res=96)
plot(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$tsne@cell.embeddings[,1],
  MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$tsne@cell.embeddings[,2],
  pch=20,
  col= tsne_palette[MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident],
  xlab=" ",
  ylab=" ",
  cex=1)
dev.off()

png(file="20.Plots/TSNE_plot_legend.png",width=400, height=1200,res=96)
plot.new()
legend("bottom",levels(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident), col=tsne_palette, pch=20,cex=1, pt.cex=3)
dev.off()

tsne_coords<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$tsne@cell.embeddings
write.table(file="10.TextData/seurat_tsne1_tsne2.tsv", tsne_coords,sep="\t",quote=F)



var_genes_norm_exp<-t(as.matrix(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@data)[MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@var.genes,])
write.table(file="10.TextData/var_genes_CPM_matrix_transposed.tsv",var_genes_norm_exp, sep="\t",quote=F)

save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)

# Seurat Diffusion Components
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj<-RunDiffusion(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,
  genes.use=MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@var.genes,
  max.dim = 3)

dmap_coords<-MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$dm@cell.embeddings[,1:2]
write.table(file="10.TextData/seurat_dc1_dc2.tsv", dmap_coords,sep="\t",quote=F)

library(scatterplot3d)
png(file="20.Plots/Diffusion_plot_with_custom_colors.png",width=1024, height=768,res=96)
scatterplot3d(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@dr$dm@cell.embeddings[,1:3],
pch=20,
box=F,
cex.symbols=1.5,
color=alpha(tsne_palette[MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident],0.75))
legend("right",levels(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj@ident),col=tsne_palette, pch=20,cex=0.7, pt.cex=1.5, bty="n")
dev.off()

save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)

# Clustering
png(file="20.Plots/Clustering.png",width=768, height=768,res=96)
MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj <- FindClusters(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,
  reduction.type = "pca",
  dims.use = 1:10,
  resolution = 0.1,
  print.output = 0,
  save.SNN = TRUE,
  plot.SNN=T)
dev.off()

save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)

png(file="20.Plots/TSNE_plot_with_clusters.png",width=900, height=800,res=96)
TSNEPlot(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,group.by="ident")
dev.off()

# Finding markers
marker_genes_for_all_clusters<-FindAllMarkers(MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,
only.pos = F,
min.pct = 0.25)

write.table(file="10.TextData/seurat_clustering_marker_genes.txt", marker_genes_for_all_clusters, sep="\t",quote=F)

library(dplyr)
top10 <- marker_genes_for_all_clusters %>% group_by(cluster) %>% top_n(10, avg_logFC)

png(file="20.Plots/seurat_cluster_heatmap_with_markers.png",width=768, height=1024,res=96)
DoHeatmap(object = MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj,
  genes.use = top10$gene,
  slim.col.label = TRUE,
  remove.key = TRUE,
  col.low="dodgerblue2",
  col.high="red3",
  col.mid="white")
dev.off()

# Saving RData
save(file="00.RData/CPM_gene_filt_seurat_obj.RData",MMlines_scRNA_seq_40k_CPM_gene_filt_seurat_obj)

####################

# Diffusion Maps with destiny
library(destiny)
DM_components<-DiffusionMap(var_genes_norm_exp,sigma="local", k=250, n_eigs= 100)

save(file="00.RData/DM_components_100_eigenvalues.RData",DM_components)

library("scatterplot3d")
attach(DM_components_eigenvectors); scatterplot3d(DM_components_eigenvectors[,1:3],pch=20,box=F,cex.symbols=1.5,color=alpha(c("burlywood1","darkorchid1","indianred1","limegreen","gold","dodgerblue2","violetred1","darkcyan","forestgreen", "darkorange", "magenta4", "hotpink", "red3", "skyblue", "darkblue","coral2","aquamarine3","greenyellow")[Experiment],0.5)); detach(DM_components_eigenvectors)
```
