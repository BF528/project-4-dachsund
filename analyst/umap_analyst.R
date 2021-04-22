# Runs differential expression analysis across clusters
# Makes violin plots, heatmap, and UMAP representing each of the clusters
# Allison Nau

# Load packages:
library(Seurat)
# Say yes to installing miniconda
library(dplyr)

prog_data <- readRDS("pbmc_seurat.rda")

# Label cluster identities (will vary for each study):
# TODO this rename cluster isn't working
#new_ids <- c("Fake Alpha", "Fake Beta", "Fake Delta", "Fake Gamma",
#             "Fake Epsilon", "Fake Ductal", "Fake Acinar", "Fake Stellate",
#             "Fake Vascular", "Fake Macrophage", "Fake Cytotoxic T", 
#             "Fake Mast", "Fake 13")

#names(new_ids) <- levels(prog_data)
#prog_data2 <- RenameIdents(prog_data, new_ids)

# Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report 
# only the positive ones
# Methods available to use: https://satijalab.org/seurat/archive/v3.0/de_vignette.html
# Default: Wilcox (used)
# Also reasonable: "LR": logistic regression
# Minimum percent: detected in a minimum proportion of cells in either of the two pops
# Defaults: logfc.threshold=0.25, min.pct=0.1
# These parameters were chosen in an attempt to balance being able to find
# classically defined markers for a cluster that may be a mixed type
# and still finding genes that are descriptive of the cell type that is in that
# cluster
diff_markers <- FindAllMarkers(prog_data, only.pos=TRUE, min.pct=0.25, 
                               logfc.threshold=0.25)

# Filter on only significant markers:
diff_markers_sign <- diff_markers[diff_markers$p_val_adj<0.05,]

# Save to csv file:
write.csv(diff_markers, file="diff_genes.csv")
write.csv(diff_markers_sign, file="diff_genes_sign.csv")

# TODO: make list of genes using, then filter diff_markers to grab only 
# genes just as descriptive and save to csv


# Look at only the top ones:
top5 <- diff_markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top8 <- diff_markers %>% group_by(cluster) %>% top_n(n=8, wt=avg_log2FC)
top1 <- diff_markers %>% group_by(cluster) %>% top_n(n=1, wt=avg_log2FC)

# Save a list of genes:
top1genes <- top1$gene
top1genes
top8genes <- top8$gene
top8genes
top5genes <- top5$gene
top5genes

# Make violin plots for genes or interest:
# myvln <- VlnPlot(prog_data, features=c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
#                                       "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
#                                       "CD163", "CD68", "IgG", "CD3", "CD8",
#                                       "TPSAB1", "KIT", "CPA3"), pt.size=0, combine=TRUE)
# TODO: add back in "IgG", "CD3", "CD8",
#png(filename="myvln.png", width=1500, height=1000)
#myvln
#dev.off()



# Make violin plot of the top expressed gene of each cluster:
vplot1 <- VlnPlot(prog_data, features=top1genes, pt.size=0, combine=TRUE)

#Save violin plot of top genes:
png(filename="vplotTop1.png", width=1200, height=800)
vplot1
dev.off()

# Make giant violin plot of top5 genes:
vplot5 <- VlnPlot(prog_data, features=top5genes, pt.size=0, combine=TRUE)

png(filename="vplotTop5.png", width=6000, height=4000)
vplot5
dev.off()

# Make giant violin plot of top8 genes:
vplot8 <- VlnPlot(prog_data, features=top8genes, pt.size=0, combine=TRUE)

png(filename="vplotTop8.png", width=12000, height=8000)
vplot8
dev.off()

# Make heatmaps:
# TODO: add dendograms?
map <- DoHeatmap(prog_data, features=top5$gene)

# Save heatmap:
png(filename="heatmapTop5.png", width=1200, height=800)
map
dev.off()


#TODO: UMAP cluster names

# Run UMAP:
#TODO: why dims 1:10
my_umap1 <- RunUMAP(prog_data, dims=1:10)
my_umap1

# Plot UMAP:
umap_plot <- DimPlot(my_umap1, reduction="umap", label=TRUE)
umap_plot

# Save UMAP:
png(filename="umap.png")
umap_plot
dev.off()