# Runs differential expression analysis across clusters
# Makes violin plots, heatmap, featureplot, and UMAP representing each of the clusters
# Cluster names are specific for this dataset and this round of anlsysis, 
# must be modified for other datasets.
# Allison Nau

# Load packages:
library(Seurat)
library(dplyr)
library(janitor)

# Read in data:
prog_data <- readRDS("pbmc_seurat.rda")

# High interest markers (based on paper and our cluster identification)
# IgG, CD3, and CD8 not in differentially expressed markers across clusters
# ETV1 from other literature to define Gamma/PP cells
# PPP1R1A, REG3A, NKX6-1 used to try to distinguish subsets
high_interest <- c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
                    "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
                    "CD163", "CD68", "TPSAB1", "KIT", "CPA3", 
                    "ETV1",
                    "PPP1R1A", "NKX6-1", "RP5-964N17.1")


# Label cluster identities (will vary for each study and each analysis):
new_ids <- c("Alpha_1", "Beta_1", "Alpha_2", "Delta", "UnKn_1", 
             "Ductal", "Stellate", "Acinar", "Beta_2", "Gamma",
             "Macrophage", "Vascular")
names(new_ids) <- levels(prog_data)
prog_data <- RenameIdents(prog_data, new_ids)


# -----------------------------------------------------------------------------
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
# TODO diff_markers_s <- FindAllMarkers(prog_data_s, only.pos=TRUE, min.pct=0.25, 
#                               logfc.threshold=0.25)
diff_markers$gene

# Filter on only significant markers:
diff_markers_sign <- diff_markers[diff_markers$p_val_adj<0.05,]
diff_markers_sign$gene

# Save to csv file:
write.csv(diff_markers, file="diff_genes.csv")
write.csv(diff_markers_sign, file="diff_genes_sign.csv")

# -----------------------------------------------------------------------------
# For each genemake list of genes using the known markers used to define the
# cell, then filter diff_markers to grab only 
# genes just as descriptive and save to csv

# Function takes cluster name, marker gene used to define that cluster, and 
# then returns dataframe with data representing genes that are just as 
# definitive for that cluster.
novel_markers <- function(cluster, marker_gene){
  # Find the significance level for that gene:
  min_sign <- diff_markers_sign$p_val_adj[diff_markers_sign$cluster == cluster & 
                                            diff_markers_sign$gene == marker_gene]
  # Then filter on cluster and significance level:
  temp_df <- diff_markers_sign[(diff_markers_sign$cluster == cluster & 
                                  diff_markers_sign$p_val_adj<=min_sign),]
  return (temp_df)
}

# Cluster 0 (alpha_1)
c0 <- novel_markers("Alpha_1", "GCG")
# Cluster 1 (Beta_1)
c1 <- novel_markers("Beta_1", "INS")
# Cluster 2 (Alpha_2)
c2 <- novel_markers("Alpha_2", "GCG")
# Cluster 3 (Delta)
c3 <- novel_markers("Delta", "SST")
# Cluster 4 (Unknown)
# TODO
# Cluster 5 (Ductal)
c5 <- novel_markers("Ductal", "KRT19")
# Cluster 6 (Stellate)
c6 <- novel_markers("Stellate", "PDGFRB")
# Cluster 7 (Acinar)
c7 <- novel_markers("Acinar", "CPA1")
# Cluster 8 (Beta_2)
c8 <- novel_markers("Beta_2", "INS")
# Cluster 9 (Gamma)
c9 <- novel_markers("Gamma", "PPY")
# Cluster 10 ("Macrophage")
c10 <- novel_markers("Macrophage", "CD68")
# Cluster 11 (Vascular)
# Which gene is of lower significance?
VWF <- diff_markers_sign$p_val_adj[diff_markers_sign$cluster == "Vascular" & 
                                       diff_markers_sign$gene == "VWF"]
VWF  # 0 , so use the other as the threshold
PECAM1 <- diff_markers_sign$p_val_adj[diff_markers_sign$cluster == "Vascular" & 
                                     diff_markers_sign$gene == "PECAM1"]
PECAM1  # 0 too
c11 <- novel_markers("Vascular", "PECAM1")

# Combine novel markers of all  # TODO: skipping cluster 4
novel_markers_all <- rbind(c0, c1, c2, c3, c5, c6, c7, c8, c9, c10, c11)

# Save to csv:
write.csv(novel_markers_all, file="novel_markers_all.csv")

# Since that is way too many for gamma, come up with a count so table can be
# abbreviated:
novel_counts <- tabyl(novel_markers_all, cluster)
write.csv(novel_counts, file="novel_counts.csv")

# TODO: Cluster 4 is NOT delta

# -----------------------------------------------------------------------------
# Look at only the top geness:
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


# -----------------------------------------------------------------------------
# Make violin plots for genes or interest:
myvln <- VlnPlot(prog_data, features=high_interest, pt.size=0, combine=TRUE)
myvln

# Save plot
png(filename="myvln.png", width=1500, height=1500)
myvln
dev.off()

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


# -----------------------------------------------------------------------------
# Make heatmap of top 5 genes:
map <- DoHeatmap(prog_data, features=top5$gene, size=5, angle=75)

# Save heatmap:
png(filename="heatmapTop5.png", width=1200, height=800)
map
dev.off()


# -----------------------------------------------------------------------------
# Run UMAP:
my_umap1 <- RunUMAP(prog_data, dims=1:10)
my_umap1

# Plot UMAP:
umap_plot <- DimPlot(my_umap1, reduction="umap", label=TRUE, pt.size=1)
umap_plot

# Save UMAP:
png(filename="umap.png")
umap_plot
dev.off()


# -----------------------------------------------------------------------------
# Make feature plot of high interest markers: 
feat_plot = FeaturePlot(prog_data, features=high_interest)
feat_plot

# Save feature plot:
png(filename="feature_plot.png", width=1200, height=1200)
feat_plot
dev.off()

# -----------------------------------------------------------------------------
# Work on distinguishing closely related clusters:
# Cluster 0 & 2 (Alpha)
clust2_to0 <- FindMarkers(prog_data, ident.1 = "Alpha_1", ident.2 = "Alpha_2", 
                          min.pct=0.5, logfc.threshold=0.5, min.diff.pct=0.5)
# Significant genes:
clust2_to0_sign <- clust2_to0[clust2_to0$p_val_adj<0.05,]
# Save csv:
write.csv(clust2_to0_sign, file="clust2_to0_sign.csv")
# Grab the top 8 genes that are different:
clust2_to0_top8 <- clust2_to0 %>% top_n(n=8, wt=avg_log2FC)
clust2_to0_top8genes <- row.names(clust2_to0_top8)
clust2_to0_top8genes
# Make violin plot of top8 genes:
clust2_to0_vplot8 <- VlnPlot(prog_data, features=clust2_to0_top8genes, pt.size=0, combine=TRUE)
png(filename="clust2_to0_vplotTop8.png", width=1000, height=1000)
clust2_to0_vplot8
dev.off()
# Make violin plot of all genes:
clust2_to0_vplot <- VlnPlot(prog_data, features=row.names(clust2_to0_sign), 
                            pt.size=0, combine=TRUE)
png(filename="clust2_to0_vplot_alpha.png", width=1000, height=1000)
clust2_to0_vplot
dev.off()


# ------------------------------------------------
# Cluster 1 & 8 (Beta)
clust8_to1 <- FindMarkers(prog_data, ident.1 = "Beta_1", ident.2 = "Beta_2", 
                          min.pct=0.5, logfc.threshold=0.5, min.diff.pct=0.5)
# Significant genes:
clust8_to1_sign <- clust8_to1[clust8_to1$p_val_adj<0.05,]
# Save csv:
write.csv(clust8_to1_sign, file="clust8_to1_sign.csv")
# Grab the top 8 genes that are different:
clust8_to1_top8 <- clust8_to1 %>% top_n(n=8, wt=avg_log2FC)
clust8_to1_top8genes <- row.names(clust8_to1_top8)
clust8_to1_top8genes
# Make violin plot of top8 genes:
clust8_to1_vplot8 <- VlnPlot(prog_data, features=clust8_to1_top8genes, 
                             pt.size=0, combine=TRUE)
png(filename="clust8_to1_vplotTop8.png", width=1000, height=1000)
clust8_to1_vplot8
dev.off()
# Make violin plot of all genes:
clust8_to1_vplot <- VlnPlot(prog_data, features=row.names(clust8_to1_sign), pt.size=0, combine=TRUE)
png(filename="clust8_to1_vplot_beta.png", width=1000, height=1000)
clust8_to1_vplot
dev.off()


# ------------------------------------------------
# Unknown (Note lower logfc.threshold, lower min_diff.pct threshold and only positive)
# Only 1 gene that is positive, the rest are LOWER
clust4 <- FindMarkers(prog_data, ident.1 = "UnKn_1", 
                          min.pct=0.5, logfc.threshold=0.5, min.diff.pct=0.2)
# Significant genes:
clust4_sign <- clust4[clust4$p_val_adj<0.05,]
# Save csv:
write.csv(clust4_sign, file="clust4_sign.csv")
# Grab the top genes that are different (only 1 in this case):
clust4_top9 <- clust4 %>% top_n(n=9, wt=avg_log2FC)
clust4_top9genes <- row.names(clust4_top9)
clust4_top9genes
# Make violin plot of top8 genes:
clust4_vplot9 <- VlnPlot(prog_data, features=clust4_top9genes, pt.size=0, combine=TRUE)
png(filename="clust4_vplotTop8.png", width=1000, height=1000)
clust4_vplot9
dev.off()
# Make violin plot of all genes:
clust4_vplot <- VlnPlot(prog_data, features=row.names(clust4_sign), pt.size=0, combine=TRUE)
png(filename="clust4_vplot.png", width=1000, height=1000)
clust4_vplot
dev.off()
