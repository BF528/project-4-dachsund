library(Seurat)
library(dplyr)

prog_data <- readRDS("GSM2230760_seurat.rda")

# Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report 
# only the positive ones
# TODO: minimum percent 0.25??
diff_markers <- FindAllMarkers(prog_data, only.pos=TRUE, min.pct=0.25, )

# Run UMAP:
#TODO: why dims 1:10
my_umap1 <- RunUMAP(prog_data, dims=1:10)
#TODO my_umap2 <- RunUMAP(prog_data, dims=1:20)

# Plot UMAP:
DimPlot(my_umap1, reduction="umap")
#TODO DimPlot(my_umap, reduction="umap2")
