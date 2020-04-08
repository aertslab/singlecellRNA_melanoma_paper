# module load R/3.3.2-foss-2014a-noX
.libPaths('/user/leuven/303/vsc30324/lcb/zkalender/R/x86_64-pc-linux-gnu-library/3.3')
library(Seurat)
# Mix_MM_57_74_87_SOX10_24h

load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Mix_MM_57_74_87_SOX10_24h_combined/Mix_MM_57_74_87_SOX10_24h/Mix_MM_57_74_87_SOX10_24h.RData")
Mix_MM_57_74_87_SOX10_24h_seurat_obj<-seurat_obj
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data/Mix_MM_57_74_87_SOX10_24h_demux_best_predictions.RData")
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("SOX10_24h")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
Mix_MM_57_74_87_SOX10_24h_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# Mix_MM_57_74_87_SOX10_48h
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Mix_MM_57_74_87_SOX10_48h_combined/Mix_MM_57_74_87_SOX10_48h/Mix_MM_57_74_87_SOX10_48h.RData")
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data/Mix_MM_57_74_87_SOX10_48h_demux_best_predictions.RData")
demux_best_predictions$SNG.1ST<-gsub("A375","NA",demux_best_predictions$SNG.1ST)
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("SOX10_48h")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
Mix_MM_57_74_87_SOX10_48h_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# Mix_MM_lines_TGF_only
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/HiSeq4000_20171123/cellranger_output/Mix_MM_lines_TGF_only/Mix_MM_lines_TGF_only.RData")
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data/Mix_MM_lines_TGF_only_demux_best_predictions.RData")
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("TGF_only")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
Mix_MM_lines_TGF_only_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# Mix_MM_lines_TGF_TNF
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/HiSeq4000_20171123/cellranger_output/Mix_MM_lines_TGF_TNF/Mix_MM_lines_TGF_TNF.RData")
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data/Mix_MM_lines_TGF_TNF_demux_best_predictions.RData")
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("TGF_TNF")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
Mix_MM_lines_TGF_TNF_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# MM057_SOX10
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/10x_MM057_SOX10.RData")
seurat_obj@meta.data$Cell_Line<-"MM057"
seurat_obj@meta.data$Experiment<-c("SOX10_72h")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM057_SOX10_seurat_obj<-seurat_obj
rm(seurat_obj)

# MM074_SOX10
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/10x_MM074_SOX10.RData")
seurat_obj@meta.data$Cell_Line<-"MM074"
seurat_obj@meta.data$Experiment<-c("SOX10_72h")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM074_SOX10_seurat_obj<-seurat_obj
rm(seurat_obj)

# MM087_NTC
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/10x_MM087_NTC.RData")
seurat_obj@meta.data$Cell_Line<-"MM087"
seurat_obj@meta.data$Experiment<-c("NTC")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM087_NTC_seurat_obj<-seurat_obj
rm(seurat_obj)

# MM087_SOX10
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/10x_MM087_SOX10.RData")
seurat_obj@meta.data$Cell_Line<-"MM087"
seurat_obj@meta.data$Experiment<-c("SOX10_72h")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM087_SOX10_seurat_obj<-seurat_obj
rm(seurat_obj)

# MM057_TL
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/MM057_timeline.RData")
seurat_obj@meta.data$Cell_Line<-"MM057"
seurat_obj@meta.data$Experiment<-c("TL")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM057_TL_seurat_obj<-seurat_obj
rm(seurat_obj)

# MM087_TL
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/MM087_timeline.RData")
seurat_obj@meta.data$Cell_Line<-"MM087"
seurat_obj@meta.data$Experiment<-c("TL")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM087_TL_seurat_obj<-seurat_obj
rm(seurat_obj)

# Mix_MM_lines
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/Mix_MM_lines.RData")
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/00.combining_data/Mix_MM_lines_demux_best_prediction.RData")
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("BL")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
Mix_MM_lines_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# MM057_MM074_NTC_Run1
load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/MM057siNTC_MM074_siNTC.RData")
load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/MM057siNTC_MM074_siNTC_demux_best_prediction.RData")
demux_best_predictions$SNG.1ST<-gsub("A375|MM031|MM047|MM099","NA",demux_best_predictions$SNG.1ST)
seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
seurat_obj@meta.data$Experiment<-c("NTC")
seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
MM057_MM074_NTC_Run1_seurat_obj<-seurat_obj
rm(seurat_obj)
rm(demux_best_predictions)

# MM057_MM074_NTC_Run2
# load("/user/leuven/303/vsc30324/lcb/zkalender/Runs/Seurat_objects_from_previous_runs/MM057siNTC_MM074_siNTC_5.RData")
# load("/user/leuven/303/vsc30324/lcb/zkalender/scRNA_seq_melanoma/MM057siNTC_MM074_siNTC_5_demux_best_prediction.RData")
# demux_best_predictions$SNG.1ST<-gsub("MM001|MM029|MM031|MM047|MM099","NA",demux_best_predictions$SNG.1ST)
# seurat_obj@meta.data$Cell_Line<-demux_best_predictions[demux_best_predictions$BARCODE %in% rownames(seurat_obj@meta.data),2]
# seurat_obj@meta.data$Experiment<-c("NTC")
# seurat_obj@meta.data$CE_ID<-paste(seurat_obj@meta.data$Cell_Line,seurat_obj@meta.data$Experiment,sep="_")
# MM057_MM074_NTC_Run2_seurat_obj<-seurat_obj
# rm(seurat_obj)
# rm(demux_best_predictions)



# Mix_MM_57_74_87_SOX10_24h
# Mix_MM_57_74_87_SOX10_48h
a<-MergeSeurat(object1=Mix_MM_57_74_87_SOX10_24h_seurat_obj, object2=Mix_MM_57_74_87_SOX10_48h_seurat_obj, add.cell.id1="Mix_MM_57_74_87_SOX10_24h", add.cell.id2="Mix_MM_57_74_87_SOX10_48h", do.normalize=F)
# Mix_MM_lines_TGF_only
a<-MergeSeurat(object1=a, object2=Mix_MM_lines_TGF_only_seurat_obj, add.cell.id2="Mix_MM_lines_TGF_only", do.normalize=F)
# Mix_MM_lines_TGF_TNF
a<-MergeSeurat(object1=a, object2=Mix_MM_lines_TGF_TNF_seurat_obj, add.cell.id2="Mix_MM_lines_TGF_TNF",do.normalize=F)
# MM057_SOX10
a<-MergeSeurat(object1=a, object2=MM057_SOX10_seurat_obj, add.cell.id2="MM057_SOX10_72h", do.normalize=F)
# MM074_SOX10
a<-MergeSeurat(object1=a, object2=MM074_SOX10_seurat_obj, add.cell.id2="MM074_SOX10_72h", do.normalize=F)
# MM087_NTC
a<-MergeSeurat(object1=a, object2=MM087_NTC_seurat_obj, add.cell.id2="MM087_NTC", do.normalize=F)
# MM087_SOX10
a<-MergeSeurat(object1=a, object2=MM087_SOX10_seurat_obj, add.cell.id2="MM087_SOX10_72h", do.normalize=F)
# MM057_TL
a<-MergeSeurat(object1=a, object2=MM057_TL_seurat_obj, add.cell.id2="MM057_TL",do.normalize=F)
# MM087_TL
a<-MergeSeurat(object1=a, object2=MM087_TL_seurat_obj, add.cell.id2="MM087_TL", do.normalize=F)
# Mix_MM_lines
a<-MergeSeurat(object1=a, object2=Mix_MM_lines_seurat_obj, add.cell.id2="Mix_MM_lines", do.normalize=F)
# MM057_MM074_NTC_Run1
a<-MergeSeurat(object1=a, object2=MM057_MM074_NTC_Run1_seurat_obj, add.cell.id2="MM057_MM074_NTC_Run1", do.normalize=F)

# MM057_MM074_NTC_Run2
MMlines_scRNA_seq_40k_seurat_obj<-a
save(file="all_data_combined.RData",MMlines_scRNA_seq_40k_seurat_obj)
