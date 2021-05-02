# Project Description

Single-cell RNA sequencing can be harnassed to elucidate the complexities of individual cells and to discover novel and rare cell types. It is useful for understanding transcriptomics at a single-cell resolution. Baron et al. (2017) implemented inDrop single-cell RNA-seq methodology to analyze over 12,000 individual pancreatic cells from four human donors and two mouse strains. This droplet-based technique allows for isolation of individual cells and incorporation of unique molecular identifiers (UMI) in preparation for RNA-sequencing. The authors sought to use this technology to better understand the heterogeneity of cell types within the pancreas in hopes of identifying previously known as well as rare and novel cell types. Characterization of this cellular diversity can be clinically applied to better understand the roles of these cell types in dysfuntion of disease.  

In this project, we sought to reproduce the original results of Baron et al. and determine 15 pancreatic cell types. To do so, we processed the raw sequencing data to generate a UMI count matrix using salmon, performed quality control and clustering of cells using Seurat, identified marker genes for distinct cell populations, and characterized the biological function of these cell clusters. 

This repository contains the scripts used to process and analyze the data.

Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell systems, 3(4), 346–360.e4. https://doi.org/10.1016/j.cels.2016.08.011

# Contributors

Data Curator: Mae Rose Gott  
Programmer: Abhishek Thakar  
Analyst: Allison Nau  
Biologist: Sheila Yee  

# Repository Contents

### Programs from Programmer, Abhishek Thakar:

#### Programmer_code.R ####
Using Seurat and tximport packages to get alevin output files and create seurat object for further downstream analysis
Created cluster sizes by cell counts and filtered data to low quality cells and low variance counts. Includes normalization prior to low variance filter step.

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
