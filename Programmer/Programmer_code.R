#set the working environment
setwd('/projectnb2/bf528/users/dachshund/project_4/project-4-dachsund/programmer')

#loads the package needed to import alevin data
#install.packages(Seurat)
#BiocManager::install("tximport")
library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(tximport)
library(EnsDb.Hsapiens.v79)
library(Matrix)
library(tidyverse)
library(patchwork)
#all packages and scripts were obtained from the publishing website provided by BU clustering computing network

#give location of alevin output file to variable
files <- '/projectnb2/bf528/users/dachshund/project_4/project-4-dachsund/data_curator/alevin_output/alevin/quants_mat.gz'
file.exists(files)

#uses txiport to import alevin file and txi variable
txi <- tximport(files, type="alevin")
tx <- txi$counts
# get ensemble ids, and remove unnecessary extensions
ensembleid <-  rownames(tx) 
ensembleid <- sub("[.][0-9]*","",ensembleid) 
rownames(tx) <- ensembleid

# convert gene names of ensembl ids 
symbols <- select(EnsDb.Hsapiens.v79, keys= ensembleid, keytype = "GENEID", columns = c("SYMBOL","GENEID")) 

# subset to rows that have the matches, and rename matrix with genenames
tx <- tx[rownames(tx) %in% symbols$GENEID,]
rownames(tx) <- symbols$SYMBOL

#view count and the rownames
tx[1:10, 1:2]
View(symbols)

#create Seurat object using txi$counts which is the cell count and assinging it to variable pbmc

pbmc <- CreateSeuratObject(counts = tx , min.cells = 3, min.features = 200, project = "10X_PBMC")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") #calculates percentages of counts and does quality check
head(pbmc@meta.data, 5)
#plot initial numbers of cell count and gene expression count per cell pre filtering and normalization
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
      
#filtering genes higher than 200 to pbmc
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20) 
#VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3) 
#normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


#account for variation
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#set a standard for gene variation
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
        
pbmc <- RunPCA(pbmc, features=VariableFeatures(object=pbmc))

#cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <-RunUMAP(pbmc, dims=1:10)
DimPlot(pbmc, reduction="umap")

saveRDS(pbmc, file="/projectnb2/bf528/users/dachshund/project_4/project-4-dachsund/programmer/pbmc_seurat.rda")

#create barplot for per cluster count
makegraph<- (table(Idents(pbmc)))
barplot(makegraph, main="Cluster sizes by cell counts", xlab="Cluster number", ylab="Cell Count")
        

