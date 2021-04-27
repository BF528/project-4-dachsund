# Project Description

A brief description of what this repository is for and what it contains

# Contributors

List contributor names and github user names, or email addresses if desired

# Repository Contents

### Programs from Analyst, Allison Nau:

#### analyst_install.R #### 
Install packages required to run umap_analyst.R

Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript analyst_install.R
```

#### umap_analyst.R ####
Performs differential expression analysis on "pbmc_seurat.rda" outputed from a previous step. 
Generates violin plots, a UMAP, feature plots, and heatmap used to identify cell types for each cluster. 
Clusters of the same major cell type are also compared to each other. The one unknown cluster is compared to 
all other clusters. This script is specific for this dataset, and parameters will have to be changed based on the 
input file.

Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript umap_analyst.R
```

Outputs: plots and tables with differential expression results.