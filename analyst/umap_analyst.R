# Load packages:
library(Seurat)
# Say yes to installing miniconda
library(dplyr)

prog_data <- readRDS("GSM2230760_seurat_anau.rda")

# Label cluster identities (will vary for each study):
# TODO this rename cluster isn't working
new_ids <- c("Fake Alpha", "Fake Beta", "Fake Delta", "Fake Gamma",
             "Fake Epsilon", "Fake Ductal", "Fake Acinar", "Fake Stellate",
             "Fake Vascular", "Fake Macrophage", "Fake Cytotoxic T", 
             "Fake Mast", "Fake 13")

names(new_ids) <- levels(prog_data)
prog_data2 <- RenameIdents(prog_data, new_ids)

# Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report 
# only the positive ones
# Methods available to use: https://satijalab.org/seurat/archive/v3.0/de_vignette.html
# Default: Wilcox
# "LR": logistic regression
# Minimum percent: detected in a minimum proportion of cells in either of the two pops
# logfc.threshold
diff_markers <- FindAllMarkers(prog_data2, only.pos=TRUE, min.pct=0.25, 
                               logfc.threshold=0.25)
# Minimum percent 0.8 came back with nothing, try 0.5
diff_markers_high_percent <- FindAllMarkers(prog_data2, only.pos=TRUE, min.pct=0.5, 
                               logfc.threshold=0.25)

diff_markers_sign <- diff_markers[diff_markers$p_val_adj<0.05,]

# Save to csv file:
write.csv(diff_markers, file="diff_genes.csv")
write.csv(diff_markers_sign, file="diff_genes_sign.csv")

# Look at only the top ones:
top5 <- diff_markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top3 <- diff_markers %>% group_by(cluster) %>% top_n(n=3, wt=avg_log2FC)
top10 <- diff_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top1 <- diff_markers %>% group_by(cluster) %>% top_n(n=1, wt=avg_log2FC)

# Save a list of genes:
top1genes <- top1$gene
top1genes

# Make violin plots for genes or interest:
# myvln <- VlnPlot(prog_data, features=c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
#                                       "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
#                                       "CD163", "CD68", "IgG", "CD3", "CD8",
#                                       "TPSAB1", "KIT", "CPA3"), pt.size=0.0001)

# TODO: add back in "IgG", "CD3", "CD8",

# Make violin plot of the top expressed gene of each cluster:
vplot1 <- VlnPlot(prog_data2, features=top1genes, pt.size=0, combine=TRUE)
vplot1

#Save violin plot of top genes:
png(filename="vplotTop1.png", width=1200, height=800)
vplot1
dev.off()

# Make heatmaps:
# TODO: add dendograms?
DoHeatmap(prog_data2, features=top3$gene)
DoHeatmap(prog_data2, features=top5$gene)
DoHeatmap(prog_data2, features=top10$gene)
map <- DoHeatmap(prog_data2, features=top5$gene)
map

# Save heatmap:
png(filename="heatmapTop5.png", width=1200, height=800)
map
dev.off()


#TODO: UMAP cluster names

# Run UMAP:
#TODO: why dims 1:10
my_umap1 <- RunUMAP(prog_data2, dims=1:10)
my_umap1

# Plot UMAP:
umap_plot <- DimPlot(my_umap1, reduction="umap", label=TRUE)
umap_plot

# Save UMAP:
png(filename="umap.png")
umap_plot
dev.off()
