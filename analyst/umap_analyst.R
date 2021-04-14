library(Seurat)
library(dplyr)

prog_data <- readRDS("GSM2230760_seurat.rda")

# Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report 
# only the positive ones
# TODO: minimum percent 0.25?? logfc.threshold
diff_markers <- FindAllMarkers(prog_data, only.pos=TRUE, min.pct=0.25, 
                               logfc.threshold=0.25)

# Look at only the top ones:
grptop <- diff_markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

#TODO: Look at VLN plots??
# TODO: 
for (feature in c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
                  "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
                  "CD163", "CD68", 
                  "TPSAB1", "KIT", "CPA3")){
temp_vln <- VlnPlot(prog_data, features=c(feature))
show(temp_vln)
# TODO fix the paste
save.image(paste(feature, ".png", sep=""))
}


# TODO: try to make plots smaller, including maybe removing plotting the individual points

myvln <- VlnPlot(prog_data, features=c("GCG", "INS",  "SST", "PPY", "GHRL", "KRT19",
                                       "CPA1", "PDGFRB", "VWF", "PECAM1", "CD34",
                                       "CD163", "CD68", "IgG", "CD3", "CD8",
                                       "TPSAB1", "KIT", "CPA3"))
myvln

#TODO: fix saving 
save.image("myvln.png")
save.image("myvln.jpg")

top3 <- diff_markers %>% group_by(cluster) %>% top_n(n=3, wt=avg_log2FC)
# TODO: try original object

# Try heatmap:
map <- DoHeatmap(prog_data, features=top10$gene)
map
save.image("map.png")
#TODO: delete? DoHeatmap(diff_markers, features=top10$gene) + NoLegend()


#TODO: fix all images

# Run UMAP:
#TODO: why dims 1:10
my_umap1 <- RunUMAP(prog_data, dims=1:10)
#TODO my_umap2 <- RunUMAP(prog_data, dims=1:20)

# Plot UMAP:
DimPlot(my_umap1, reduction="umap", label=TRUE)
#TODO DimPlot(my_umap, reduction="umap2")
