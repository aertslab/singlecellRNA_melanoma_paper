## Scripts related to the manuscript

### single cell migration
#### 
	single_cell_migration.Rmd

### scRNA-seq Data 10x
	#### CellRanger 
	10xMix_MM_lines_cellranger.sh

	#### Demultiplexing of pooled samples with demuxlet
	Mix_MM_57_74_87_SOX10_24h_combined_demuxlet.txt
	HiSeq4000_20171123_parse_demuxlet_output.R

	#### Analysis with Seurat
	all_analysis_combination_of_matrices.R
	seurat_processing_40K.R
	
	#### aggregate single-cell counts per cell line / treatment / melanoma cell state
	aggregate_and_normalize_counts.R

	#### DC trajectory alignment
	functions.R
	Melanoma_Trajectory_Alignment_Scripts.Rmd
	Melanoma_Trajectory_Alignment_Scripts.Rmd
	Melanoma_Heatmaps.Rmd
	
	#### running pySCENIC
	Melanoma_pySCENIC_100x_Scripts_new_trk.Rmd
	
	#### get the regulons for each of the 100 runs (motif and track regulons separately)
	get_regulons_from_pyscenic_table.R
		
	#### make data base with AUCell rankings for each 10x subset
	make_AUCell_ranking_db.R

	#### filter regulons from 100 scenic runs (per 10x subset)
	filter_regulons.R
	
 	#### boxplots AUCell scores on motif and track regulons 100 pySCENIC runs subset 10 baselines
	calc_AUCell_100runs.R

	#### run functional enrichment on regulons (10 baselines subset)
	regulon_functional_enrichment.R

	#### filter regulons to networks
	filter_regulons_to_network.R
	
	#### special selection 10 baselines regulons
	ten_baselines_selected_regulons_aucell_heatmap_and_violinplots.R
	
	#### UMAP on AUCell values filtered regulons
	calculate_aucell_umaps.R
	
	#### calculate AUCell scores on public gene signatures subset 10 baselines
	calculate_aucell_signatures_ten_baselines.R
    
	#### calculate AUCell scores on public gene signatures knockdowns
	calculate_aucell_signatures_knockdowns.R
	
	#### statistical test of differential expression of genes / gene sets / regulons
	wilcoxon_regulons_ten_baselines.R
	wilcoxon_signatures_ten_baselines.R
	wilcoxon_genes.R
 
### scRNA-seq Data Drop-seq
	#### Nemesh script
	DropseqMM057_NemeshPipeline_SOX10.sh
	
	#### Analysis with Seurat
	DropSeqMM057_Seurat.R

	#### pySCENIC
	DropSeqMM057_runGRNboost2.py
	DropSeqMM057_pySCENIC.txt

### public scRNA-seq data
	#### score tumor biopsy data from Jerby-Arnon et al. 2018 (GSE115978)
	SKCM_GSE115978_score_with_MMregulons.R

	#### run pySCENIC on tumor biopsy data from Tirosh et al. 2016 (GSE72056)
	MM-GSE72056_pySCENIC.txt

### bulk RNA-seq
  ####
  	bulk_RNAseq_map_count.sh
	Bulk_RNA_new_MMlines.Rmd
  	Bulk_RNA_TGFbTNFa.Rmd
  	Bulk_RNA_THZ2.Rmd
  	Bulk_RNA_KDs.Rmd

### OmniATAC-seq
  ####
  	omni_atac_peak_calling.sh
  	link_atac_peaks_to_genes.R
  	ATAC_with_AC_hg19_500bp.sh

