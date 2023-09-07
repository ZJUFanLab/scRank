# scRank
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/)

<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/workflow.png'>


Cells respond divergently to drugs due to the heterogeneity among cell populations,thus it is crucial to identify the drug-responsive cell population for accurately
elucidating the mechanism of drug action, which is a great challenge yet. Here, we
address it with scRank using a target-perturbed gene regulatory network (tpGRN) to rank and infer drug-responsive cell population towards in-silico drug perturbation for single-cell transcriptomic data under disease condition. scRank enables the inference of drug-responsive cell types for single-cell data under disease condition, providing new insights into the mechanism of drug action. 

## Installation
- install dependent packages `devtools` and [`rTensor`](https://github.com/rikenbit/rTensor)
)
```{r}
#install.packages("devtools")
#devtools::install_github("rikenbit/rTensor")
devtools::install_github("ZJUFanLab/scRank")
```

## Overview
scRank method consists of two components, wherein the first is to reconstruct the gene regulatory network from expression ptrofiles using `Constr_net` function and the second step is to estimate the effect of the in silico drug perturbation for GRNs in each cell type using `rank_celltype` function. 

scRank start with create a S4 object by `CreateScRank` function:
- the `input` is the gene expression profil eand `meta` is the cell type information. 
- `cell_type` is the column name of the cell type information in `meta` 
- `species` is the species of the data. ("mouse" or "human")
- `drug` is the drug name and `target` is the target gene of the drug. `drug` could be any inhibitor in our database `utile_database`. if you know the specific target gene of the drug, you can input the target gene into `target` without inputing `drug`.

```{r}
CreateScRank <- function(input,
                         meta,
                         cell_type,
                         species,
                         drug,
                         target)
```

The format of the `input` is as follows:
1. gene expression profile formatted by matrix or data frame, where the column is gene and the row is cell.
2. Seurat object with metadata containing cell type information

The `meta` is required if `input` is not a Seurat objectas, where its format as follows:
1. a dataframe with row names as cell names matched with column names of `input` and column names as cell type information cooresponding to the `cell_type` argument.

## Tutorial
In this tutorial, we will demonstrate how to  infer the drug-responsive cell type by scRank based on a demo dataset (GSE110894) containing BET inhibitor resistant and sensitive leukaemic cells.

<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/original_data.png'>

### 1. Load the data and create a scRank object
we load the demo dataset from Seurat object, the drug target is known as Brd4.



```{r}
seuratObj <- system.file("extdata", "AML_object,rda", package="scRank")
obj <- CreateScRank(input = seuratObj,
                    species = 'mouse',
                    cell_type = 'label',
                    target = 'Brd4')
```

### 2. Construct the gene regulatory network
```{r}
obj <- scRank::Constr_net(obj)
```

### 3. Rank the cell types
```{r}
obj <- scRank::rank_celltype(obj)
```

the final infered rank of cell types that determine the drug response is stored in `obj@cell_type_rank`

### 4. Visualize the result
For visulizing the rank of cell types in dimension reduction space, we can use the `plot_dim` function.

```{r}
plot_dim(obj)
```
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/scRank_data.png'>

For visulizing the modularized drug-target-gene related subnetwork in specific cell type, we can use the `plot_net` function, where the parameter `mode` can be "heatmap" or "network" for different visualization.

```{r}
plot_net(obj, mode = "heatmap", cell_type = "sensitive")
plot_net(obj, mode = "heatmap", cell_type = "resistant")
```
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/sensitive_net.png'>
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/resistant_net.png'>

