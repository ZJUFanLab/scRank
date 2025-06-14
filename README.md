# scRank
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/)

<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/workflow.png'>


Cells respond divergently to drugs due to the heterogeneity among cell populations,thus it is crucial to identify the drug-responsive cell population for accurately
elucidating the mechanism of drug action, which is a great challenge yet. Here, we
address it with scRank using a target-perturbed gene regulatory network (tpGRN) to rank and infer drug-responsive cell population towards in-silico drug perturbation for single-cell transcriptomic data under disease condition. scRank enables the inference of drug-responsive cell types for single-cell data under disease condition, providing new insights into the mechanism of drug action. 

## Installation
- install dependent packages `devtools` and [`rTensor`](https://github.com/rikenbit/rTensor)

```{r}
#install.packages("devtools")
#devtools::install_github("rikenbit/rTensor")
devtools::install_github("ZJUFanLab/scRank")
```
## Note

### New
- Add `subsetScRank()` function. see [issue #14](https://github.com/ZJUFanLab/scRank/issues/14)

### Key Updates
- **Disease Relevance and Drug Effects Analysis:** Introducing the new `scRank_GSEA()` and `plot_drug_function()` functions for analyzing disease relevance and drug effects in the highest-ranking cell types.
- **Drug Type Specification:** Added a `type` parameter in `rank_celltype()` to specify modeling effects of either agonists or antagonists, enhancing the versatility of drug response modeling.
- **Efficient Large Matrix Manipulations:** Integration of the Python module "tensorly" in the `Constr_net()` function, with a new parameter `use_py`, to optimize large-scale data processing.
- **Enhanced Cell State Discernment:** Integration of the `scSHC` algorithm into the `CreateScRank()` function with an `if_cluster` parameter, improving the tool's ability to discern various cell states. [More about scSHC](https://github.com/igrabski/sc-SHC).
- **Incorporating Drug Resistance Mechanisms:** The `resistance_target` parameter in `rank_celltype()` allows for inputting targets of alternative pathways, aiding in the consideration of drug resistance mechanisms.
- **Flexible Edge Weight Adjustment:** Introduction of the `keep_ratio` parameter to adjust edge weights in the gene regulatory network, allowing for differential treatment of node types.

### To-Do
- **Packaging and Accessibility:** We are in the process of submitting scRank to Bioconductor or CRAN for enhanced accessibility.
## Overview
scRank method consists of two components, wherein the first is to reconstruct the gene regulatory network from expression ptrofiles using `Constr_net` function and the second step is to estimate the extent of the in silico drug perturbation for GRNs in each cell type using `rank_celltype` function. 

scRank start with create a S4 object by `CreateScRank` function:
- the `input` is the gene expression profil eand `meta` is the cell type information. 
- `cell_type` is the column name of the cell type information in `meta` 
- `species` is the species of the data. ("mouse" or "human")
- `drug` is the drug name and `target` is the target gene of the drug. `drug` could be found in our database `utile_database`. if you know the specific target gene of the drug, you can input the target gene into `target` without inputing `drug`.
- `type` characters meaning the MOAs of drug including antagonist or agonist. Default is antagonist.
- `if_cluster` A logical meaning whether clustering single-cell transcriptomic data. Default is `FALSE`.

```{r}
CreateScRank <- function(input,
                         meta,
                         cell_type,
                         species,
                         drug,
                         target,
                         type,
                         if_cluster,
                         var.genes)
```

The format of the `input` is as follows:
1. gene expression profile formatted by matrix or data frame, where the column is gene and the row is cell.
2. Seurat object with metadata containing cell type information

The `meta` is required if `input` is not a Seurat objectas, where its format as follows:
1. a dataframe with row names as cell names matched with column names of `input` and column names as cell type information cooresponding to the `cell_type` argument.

## Tutorial
In this tutorial, we will demonstrate how to  infer the drug-responsive cell type by scRank based on a demo dataset (GSE110894) containing BET inhibitor resistant and sensitive leukaemic cells.

<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/original_data.png' width='500'>

### 1. Load the data and create a scRank object
we load the demo dataset from Seurat object, the drug target is known as Brd4.



```{r}
seuratObj <- system.file("extdata", "AML_object.rda", package="scRank")
load(seuratObj)
obj <- CreateScRank(input = seuratObj,
                    species = 'mouse',
                    cell_type = 'labels',
                    target = 'Brd4')
```

### 2. Construct the gene regulatory network
```{r}
obj <- Constr_net(obj)
```

### 3. Rank the cell types
```{r}
obj <- rank_celltype(obj)
```

the final infered rank of cell types that determine the drug response is stored in `obj@cell_type_rank`

### 4. Visualize the result
To visualize the rank of cell types in dimension reduction space, we can use the `plot_dim` function.

```{r}
plot_dim(obj)
```
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/scRank_data.png' width='500'>

To visualize the modularized drug-target-gene subnetwork in a specific cell type, start by initializing the modularization using the `init_mod` function. Afterward, you can use the `plot_net` function to display the network. The `mode` parameter in `plot_net` allows you to choose between a "heatmap" or a "network" for different types of visualization.

```{r}
obj <- init_mod(obj)

plot_net(obj, mode = "heatmap", celltype = "sensitive")
plot_net(obj, mode = "heatmap", celltype = "resistant")
```
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/sensitive_net.png' width='500'>
<img src='https://github.com/ZJUFanLab/scRank/blob/main/img/resistant_net.png' width='500'>

## Reference
Chengyu Li, et al. scRank infers drug-responsive cell types from untreated scRNA-seq data using a target-perturbed gene regulatory network, Cell Reports Medicine, 2024, https://doi.org/10.1016/j.xcrm.2024.101568
