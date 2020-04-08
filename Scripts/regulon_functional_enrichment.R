# task: run clusterProfiler on filtered set of regulons 10 baselines subset

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
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(parallel)

sample <- "12.ten_lines_BL"
wd <- file.path(staging, 
	"kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets", sample,
	"RegFunction")
dir.create(wd, recursive = TRUE)
setwd(wd)

source("dotplot2.R")	#modified dotplot function from enricher package to reorder regulons

regdir <- file.path(staging,
	"kspan/scRNA_seq_melanoma/analysis_40k/MM40k100xSCENIC/MM40k100xSCENIC_woDoublets", sample)

read.gmt <- function (fileName) {
	tmp <- unname(sapply(readLines(fileName), function(x) strsplit(as.character(x), "\t")))
	tmp <- tmp[which(lengths(tmp) > 0)]
	names(tmp) <- sapply(tmp, function(x) x[1])
	lapply(tmp, function(x) x[3:length(x)])
}

# prepare dbs {
## get MSigDB hallmark gene sets
#H: hallmark gene sets
#C1: positional gene sets
#C2: curated gene sets
#C3: motif gene sets
#C4: computational gene sets
#C5: GO gene sets
#C6: oncogenic signatures
#C7: immunologic signatures
m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, human_gene_symbol)
m_df$gs_name <- gsub("HALLMARK_", "", m_df$gs_name)

## get all expressed genes
universe <- read.table(file.path(staging, "kspan/10xCombined/Mix_MM_lines/outs/filtered_gene_bc_matrices/hg19/genes.tsv"))
entrez_universe <- mapIds(org.Hs.eg.db, keys = universe$V1, keytype = "ENSEMBL", column = "ENTREZID")
entrez_universe <- unique(entrez_universe[!is.na(entrez_universe)])
symbol_universe <- unique(universe$V2)
#}

# [1] run on regulon set all_pos_regulons_recurrent_100x_jw_new_trk {
# get regulons & convert ids {
regulons <- readRDS(file.path(regdir, "all_pos_regulons_recurrent_100x_jw_new_trk_regulons.rds.gz"))
#0498718172fc020d4625e3d2ac9df067645a5d2a all_pos_regulons_recurrent_100x_jw_new_trk_regulons.rds.gz

entrez_regulons <- lapply(regulons, function(reg) {
	x <- as.vector(mapIds(org.Hs.eg.db, keys = reg, keytype = "SYMBOL", column = "ENTREZID"))
	x[!is.na(x)]
	})
uniprot_regulons <- lapply(regulons, function(reg) {
	x <- as.vector(mapIds(org.Hs.eg.db, keys = reg, keytype = "SYMBOL", column = "UNIPROT"))
	x[!is.na(x)]
	})
#}
save.image("input_and_dbs.RData.gz", compress = TRUE)

enrichments <- list()
## KEGG all regulons
#################### {
#https://yulab-smu.github.io/clusterProfiler-book/chapter11.html
#ck[["kegg"]] <- compareCluster(geneCluster = entrez_regulons, 
#	fun = "enrichKEGG", universe = entrez_universe,
#	qvalueCutoff = 0.01)
#head(as.data.frame(ck))
#x <- as.data.table(ck)
#length(unique(x$Cluster))
#	#129
#	#out of 364 have some enrichment
#dotplot(ck, showCategory = 5)

### process in parallel
enr <- mclapply(names(entrez_regulons), function(regname) {
	tryCatch({
		thisenr <- compareCluster(geneCluster = entrez_regulons[regname], fun = "enrichKEGG", 
				universe = entrez_universe, qvalueCutoff = 0.01)
     	}, error = function(e) { print(e)
        }, warning = function(w) { print("Warning")
        }, finally = { print("Done")
        })
	}, mc.cores = 30)
#rm errors
errors <- sapply(enr, function(e) {
	class(e)[[1]]!="compareClusterResult"
})
enr <- enr[!errors]

enrdf <- enr[[1]]@compareClusterResult
for (i in 2:length(enr)) {
	enrdf <- rbind(enrdf, enr[[i]]@compareClusterResult)
}
res <- new("compareClusterResult",
	compareClusterResult = enrdf,
	geneClusters = entrez_regulons,
	.call = match.call(expand.dots=TRUE)
)

pdf("all_pos_regulons_recurrent_100x_jw_new_trk_enrichKEGG_dotplot_top3.pdf", width = 20, height = 12)
	p <- dotplot(res, showCategory = 3, font.size = 10)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle("KEGG top 3")
dev.off()

enrichments[["kegg"]] <- res
save(enrichments, file = "all_pos_regulons_recurrent_100x_jw_new_trk_enrichments.rds.gz", compress = "gzip")
#}

## GO all regulons
################## {
#ck[["go"]] <- compareCluster(geneCluster = regulons, 
#	fun = "enrichGO", OrgDb = "org.Hs.eg.db", universe = symbol_universe,
#	ont = "BP", keyType = "SYMBOL")
#head(as.data.frame(ck[["go"]]))

#pdf("all_mtf_regulons_enrichGO_BP_dotplot_top3.pdf", width = 20, height = 12)
#	p <- dotplot(ck[["go"]], showCategory = 3, font.size = 9)
#	p$data[, "Cluster"] <- vapply(strsplit(
#			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
#	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
#	ggplot2::ggtitle("GO biologcal process top 3")

### process in parallel
enr <- mclapply(names(regulons), function(regname) {
	tryCatch({
		thisenr <- compareCluster(geneCluster = regulons[regname], fun = "enrichGO", 
				OrgDb = "org.Hs.eg.db", universe = symbol_universe, on = "BP", 
				keyType = "SYMBOL")
     	}, error = function(e) { print(e)
        }, warning = function(w) { print("Warning")
        }, finally = { print("Done")
        })
	}, mc.cores = 30)
#rm errors
errors <- sapply(enr, function(e) {
	class(e)[[1]]!="compareClusterResult"
})
table(errors)
	#FALSE  TRUE 
  	#223   141 
enr <- enr[!errors]

enrdf <- enr[[1]]@compareClusterResult
for (i in 2:length(enr)) {
	enrdf <- rbind(enrdf, enr[[i]]@compareClusterResult)
}

res <- new("compareClusterResult",
	compareClusterResult = enrdf,
	geneClusters = entrez_regulons,
	.call = match.call(expand.dots=TRUE)
)

pdf("all_pos_regulons_recurrent_100x_jw_new_trk_enrichGO_dotplot_top2.pdf", width = 20, height = 12)
	p <- dotplot(res, showCategory = 2, font.size = 9)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle("GO biological process top 3")
dev.off()

enrichments[["go"]] <- res
save(enrichments, file = "all_pos_regulons_recurrent_100x_jw_new_trk_enrichments.rds.gz", compress = "gzip")


## plot in two batches
# replace very long name(s)
enrdf[which(enrdf$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"), "Description"] <- "adaptive immune response, immunoglobulin"

# split data frame for plotting purpose
enrdf1 <- enrdf[enrdf$Cluster %in% unique(enrdf$Cluster)[1:111],]
enrdf2 <- enrdf[enrdf$Cluster %in% unique(enrdf$Cluster)[112:223],]

res1 <- new("compareClusterResult",
	compareClusterResult = enrdf1,
	geneClusters = entrez_regulons,
	.call = match.call(expand.dots=TRUE)
)
pdf("all_pos_regulons_recurrent_100x_jw_new_trk_enrichGO_dotplot_top2_1.pdf", width = 20, height = 16)
	p <- dotplot(res1, showCategory = 2, font.size = 9)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle("GO biological process top 3")
dev.off()

res2 <- new("compareClusterResult",
	compareClusterResult = enrdf2,
	geneClusters = entrez_regulons,
	.call = match.call(expand.dots=TRUE)
)

pdf("all_pos_regulons_recurrent_100x_jw_new_trk_enrichGO_dotplot_top2_2.pdf", width = 20, height = 16)
	p <- dotplot(res2, showCategory = 2, font.size = 9)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle("GO biological process top 3")
dev.off()
#}

## MSigDB
######### {
#ck[["msigdb"]] <- compareCluster(geneCluster = regulons, 
#	fun = "enricher", TERM2GENE = m_df)

#pdf("all_mtf_regulons_enrich_Hallmark_dotplot_top2.pdf", width = 20, height = 12)
#	p <- dotplot(cmh, showCategory = 2, font.size = 9)
#	p$data[, "Cluster"] <- vapply(strsplit(
#			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
#	p + ggplot2::theme(axis.text.x=element_text(angle=45, hjust=1, size = 7)) +
#		ggplot2::ggtitle("hallmark signatures top 2")
#dev.off()

### process in parallel
enr <- mclapply(names(regulons), function(regname) {
	tryCatch({
		thisenr <- compareCluster(geneCluster = regulons[regname], fun = "enricher", 
				TERM2GENE = m_df)
     	}, error = function(e) { print(e)
        }, warning = function(w) { print("Warning")
        }, finally = { print("Done")
        })
	}, mc.cores = 30)
#rm errors
errors <- sapply(enr, function(e) {
	class(e)[[1]]!="compareClusterResult"
})
table(errors)
	#FALSE  TRUE 
	#  180   184 
enr <- enr[!errors]

enrdf <- enr[[1]]@compareClusterResult
for (i in 2:length(enr)) {
	enrdf <- rbind(enrdf, enr[[i]]@compareClusterResult)
}

res <- new("compareClusterResult",
	compareClusterResult = enrdf,
	geneClusters = entrez_regulons,
	.call = match.call(expand.dots=TRUE)
)

pdf("all_pos_regulons_recurrent_100x_jw_new_trk_enrichHallmark_dotplot_top3.pdf", width = 24, height = 8)
	p <- dotplot(res, showCategory = 3, font.size = 9)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle("KEGG top 3")
dev.off()

enrichments[["hallmark"]] <- res
save(enrichments, file = "all_pos_regulons_recurrent_100x_jw_new_trk_enrichments.rds.gz", compress = "gzip")
#}
saveRDS(enrichments, file = "all_pos_regulons_recurrent_100x_jw_new_trk_enrichments.rds.gz", 
	compress = "gzip")

# plot top20 aucell regulons from calc_AUCell_100runs.R {
auc <- readRDS(file.path(regdir, "MM40kRSS_tenBL/medianAUC_100runs_three_states.rds.gz"))

nlevels <- length(unique(auc$rn))
cutoff <- 80
recregs <- as.vector(auc[, .N, by = .(variable)][N >= cutoff*nlevels, variable])
print(paste0("regulons found >= ", cutoff, "x: ", length(recregs)))

auc_100 <- auc[variable %in% recregs]

#calculate zscore per regulon
auc_100[, zscore := scale(value), by = variable]

lines <- unique(auc_100$rn)
top20 <- sapply(lines, function(line) {
	regorder <- as.vector(auc_100[rn == line, median(zscore, na.rm = T), by = variable][order(V1, decreasing = T), variable])
	gsub("_pos", "", regorder[1:20])
}, simplify = FALSE, USE.NAMES = TRUE)

#subset enrichment to top regulons
subenr <- sapply(enrichments, function(db) {
	topenr <- sapply(top20, function(topregs) {
		subdb <- db
		ccr <- subdb@compareClusterResult
		ccr <- ccr[ccr$Cluster %in% topregs,]

		#order by topregs
		ccr$Cluster <- factor(ccr$Cluster, levels = topregs)
		levels(ccr$Cluster) <- topregs
		ccr <- ccr[order(ccr$Cluster),]
		
		#replace very long GO term with shorter one
		ccr[which(ccr$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"), "Description"] <- "adaptive immune response, immunoglobulin"

		terms <- ccr[which(ccr$Description %like% "antigen processing and presentation of peptide"), "Description"]
		terms <- gsub("antigen processing and presentation of peptide", "antigen pept.", terms)
		ccr[which(ccr$Description %like% "antigen processing and presentation of peptide"), "Description"] <- terms

		subdb@compareClusterResult <- ccr
		subdb
	}, simplify = FALSE, USE.NAMES = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE)

pdf("all_pos_regulons_recurrent_100x_jw_new_trk_top20_per_group_enrichKEGG_GO_Hallmark_top3.pdf", width = 10)
for (db in names(subenr)) {
	for (line in names(subenr[[db]])) {

	p <- dotplot(subenr[[db]][[line]], showCategory = 3, font.size = 9)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p$data[, "Cluster"] <- factor(p$data[, "Cluster"],
		levels = top20[[line]][top20[[line]] %in% p$data[, "Cluster"]])
	p <- p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle(paste0(db, " enr top20 regs ", line, " top 3 terms"))
	print(p)
	}
}
dev.off()
#}

for (db in names(enrichments)) {
	print(db)
	tab <- enrichments[[db]]@compareClusterResult
	tab <- tab[order(tab$Cluster),]
	write.table(tab, file = paste0("functional_enrichment_regsI_", db, ".tsv"),
		quote = F, sep = "\t", row.names = F)
}
#}

# [2] rerun complete analysis selected regulons {
regulons <- readRDS(file.path(regdir, "all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk_regulons.rds.gz"))

entrez_regulons <- lapply(regulons, function(reg) {
	x <- as.vector(mapIds(org.Hs.eg.db, keys = reg, keytype = "SYMBOL", column = "ENTREZID"))
	x[!is.na(x)]
	})

enrichments <- list()
enrichments[["kegg"]] <- compareCluster(geneCluster = entrez_regulons, 
		fun = "enrichKEGG", universe = entrez_universe,
		qvalueCutoff = 0.01)
enrichments[["go"]] <- compareCluster(geneCluster = regulons, 
		fun = "enrichGO", OrgDb = "org.Hs.eg.db", universe = symbol_universe,
		ont = "BP", keyType = "SYMBOL")
enrichments[["msigdb"]] <- compareCluster(geneCluster = regulons, 
		fun = "enricher", TERM2GENE = m_df)
saveRDS(enrichments, file = "functional_enrichment_all_pos_regulons_recurrent_100x_jw_new_trk_plusSelRegs_andUSF2trk.rds.gz", compress = "gzip")

# write to tables
for (db in names(enrichments)) {
	print(db)
	tab <- enrichments[[db]]@compareClusterResult
	tab <- tab[order(tab$Cluster),]
	write.table(tab, file = paste0("functional_enrichment_all_regulons_", db, ".tsv"),
		quote = F, sep = "\t", row.names = F)
}
#}

# [3] run on signatures {
#signatures from jasper slack 200124
#${staging}kspan/scRNA_seq_melanoma/analysis_40k/scMMpaper/scRNA_paper_20200123.gmt

sigs <- read.gmt(file.path(staging, "kspan/scRNA_seq_melanoma/analysis_40k/scMMpaper/scRNA_paper_20200123.gmt"))
#take non-tf signatures
sigs <- sigs[70:81]
names(sigs)
# [1] "INT_cell_state_old_bulk_RNA"         "MEL_cell_state_old_bulk_RNA"        
# [3] "Curated_INT_cell_state_old_bulk_RNA" "Curated_MEL_cell_state_old_bulk_RNA"
# [5] "INT_vs_other_old_bulk_RNA"           "other_vs_INT_old_bulk_RNA"          
# [7] "New_MM_lines_PTM_pos_corr_p0.05"     "New_MM_lines_PTM_neg_corr_p0.05"    
# [9] "New_MM_lines_Cohort_ B"              "Curated_new_MM_lines_Cohort_ B"     
#[11] "New_MM_lines_MES"                    "New_MM_lines_MEL"               

entrez_sigs <- lapply(sigs, function(reg) {
	x <- as.vector(mapIds(org.Hs.eg.db, keys = reg, keytype = "SYMBOL", column = "ENTREZID"))
	x[!is.na(x)]
	})

enrichments <- list()
enrichments[["kegg"]] <- compareCluster(geneCluster = entrez_sigs, 
		fun = "enrichKEGG", universe = entrez_universe,
		qvalueCutoff = 0.01)
enrichments[["go"]] <- compareCluster(geneCluster = sigs, 
		fun = "enrichGO", OrgDb = "org.Hs.eg.db", universe = symbol_universe,
		ont = "BP", keyType = "SYMBOL")
enrichments[["msigdb"]] <- compareCluster(geneCluster = sigs, 
		fun = "enricher", TERM2GENE = m_df)
saveRDS(enrichments, file = "signatures_enrichments.rds.gz", compress = "gzip")

pdf("signatures_enrichKEGG_GO_Hallmark_top3.pdf", width = 10)
for (db in names(enrichments)) {
	p <- dotplot(enrichments[[db]], showCategory = 3, font.size = 10)

	terms <- p$data[, 'Description']
	terms <- gsub("antigen processing and presentation", "antig. proc. and pres.", terms)
	p$data[, 'Description'] <- terms

	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p <- p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle(paste0(db, " signatures top 3 terms"))
	print(p)
}
dev.off()

pdf("signatures_enrichKEGG_GO_Hallmark_top5.pdf", width = 12)
for (db in names(enrichments)) {
	p <- dotplot(enrichments[[db]], showCategory = 5, font.size = 10)

	terms <- p$data[, 'Description']
	terms <- gsub("antigen processing and presentation", "antig. proc. and pres.", terms)
	p$data[, 'Description'] <- terms


	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p <- p + ggplot2::theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) +
		ggplot2::ggtitle(paste0(db, " signatures top 5 terms"))
	print(p)
}
dev.off()


### plot only selected signatures
#INT_cell_state_old_bulk_RNA
#New_MM_lines_MEL
#New_MM_lines_MES
selsigs <- c("INT_cell_state_old_bulk_RNA", "New_MM_lines_MEL", "New_MM_lines_MES")
enrichments <- readRDS("signatures_enrichments.rds.gz")

enr <- lapply(enrichments, function(enr) {
	ccr <- enr@compareClusterResult
	ccr <- ccr[ccr$Cluster %in% selsigs,]
	ccr <- droplevels(ccr)

	#take top 10 terms for each cluster
	topterms <- lapply(levels(ccr$Cluster), function(l) {
		x <- ccr[ccr$Cluster == l, ]
		x <- x[order(x$qvalue),]
		x[1:10,]
	})
	topterms <- do.call('rbind', topterms)
})
enr <- do.call('rbind', enr)
enr <- enr[!is.na(enr$Cluster),]

dbs <- vapply(strsplit(rownames(enr), "\\."), '[', 1, FUN.VALUE=character(1))
dbs <- paste0(toupper(dbs), ": ")
enr$Description <- paste0(dbs, enr$Description)

#issue: dotplot appends #genes in database to gene set
#-> problem when concatenating results of different enrichments with different dbs
gr <- enr$GeneRatio
grr <- as.vector(sapply(gr, function(s) {eval(parse(text = s))}))
grr <- round(grr*1500)
grr <- paste0(grr, "/1500")
enr$GeneRatio <- grr

res <- new("compareClusterResult",
	compareClusterResult = enr,
	geneClusters = entrez_sigs[selsigs],
	.call = match.call(expand.dots=TRUE)
)
pdf("selected_signatures_enrichKEGG_GO_Hallmark.pdf", height = 10, width = 10)
	p <- dotplot(res, showCategory = 100, font.size = 10)
	p$data[, "Cluster"] <- vapply(strsplit(
			as.vector(p$data[, "Cluster"]), "\n"), '[', 1, FUN.VALUE=character(1))
	p <- p + ggplot2::theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 10)) +
		ggplot2::ggtitle("enrichment KEGG/GO/MSigDB hallmark")
	print(p)
dev.off()

# write to table
for (db in names(enrichments)) {
	print(db)
	tab <- enrichments[[db]]@compareClusterResult
	tab <- tab[order(tab$Cluster),]
	write.table(tab, file = paste0("functional_enrichment_signatures_", db, ".tsv"),
		quote = F, sep = "\t", row.names = F)
}
#}

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
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets
#[8] methods   base
#
#other attached packages:
# [1] ggplot2_3.2.1          dplyr_0.8.3            msigdbr_7.0.1
# [4] org.Hs.eg.db_3.8.2     AnnotationDbi_1.46.1   IRanges_2.18.3
# [7] S4Vectors_0.22.1       Biobase_2.44.0         BiocGenerics_0.30.0
#[10] clusterProfiler_3.12.0 data.table_1.12.4
#
#loaded via a namespace (and not attached):
# [1] viridis_0.5.1       httr_1.4.1          tidyr_0.8.3
# [4] tidygraph_1.1.2     bit64_0.9-7         jsonlite_1.6
# [7] viridisLite_0.3.0   splines_3.6.1       ggraph_2.0.0
#[10] assertthat_0.2.1    DO.db_2.9           BiocManager_1.30.4
#[13] rvcheck_0.1.7       triebeard_0.3.0     urltools_1.7.3
#[16] blob_1.2.0          progress_1.2.2      ggrepel_0.8.1
#[19] pillar_1.4.2        RSQLite_2.1.2       backports_1.1.5
#[22] lattice_0.20-38     glue_1.3.1          digest_0.6.21
#[25] RColorBrewer_1.1-2  polyclip_1.10-0     qvalue_2.16.0
#[28] colorspace_1.4-1    cowplot_1.0.0       Matrix_1.2-17
#[31] plyr_1.8.4          pkgconfig_2.0.3     purrr_0.3.2
#[34] GO.db_3.8.2         scales_1.0.0        ggplotify_0.0.4
#[37] europepmc_0.3       tweenr_1.0.1        enrichplot_1.4.0
#[40] BiocParallel_1.18.1 ggforce_0.3.1       tibble_2.1.3
#[43] farver_2.0.1        withr_2.1.2         UpSetR_1.4.0
#[46] lazyeval_0.2.2      magrittr_1.5        crayon_1.3.4
#[49] memoise_1.1.0       DOSE_3.10.2         MASS_7.3-51.4
#[52] xml2_1.2.2          prettyunits_1.0.2   tools_3.6.1
#[55] hms_0.5.1           stringr_1.4.0       munsell_0.5.0
#[58] compiler_3.6.1      gridGraphics_0.4-1  rlang_0.4.0
#[61] ggridges_0.5.1      grid_3.6.1          igraph_1.2.4.1
#[64] labeling_0.3        gtable_0.3.0        DBI_1.0.0
#[67] reshape2_1.4.3      graphlayouts_0.5.0  R6_2.4.0
#[70] gridExtra_2.3       bit_1.1-14          zeallot_0.1.0
#[73] fastmatch_1.1-0     fgsea_1.10.1        GOSemSim_2.10.0
#[76] stringi_1.4.3       Rcpp_1.0.2          vctrs_0.2.0
#[79] tidyselect_0.2.5
#}
