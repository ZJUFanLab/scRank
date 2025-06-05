#' @title  Definition of 'scRank' class
#' 
#' @description An S4 class containing the data, meta, and results of inferred GRN network, gene distance and cell type rank.
#' @slot data A list containing the raw and quality controled data.
#' @slot meta A list containing the raw and revised meta data.
#' @slot para A list containing the parameters.
#' @slot net A list containing the gene regulatory network of each cell type.
#' @slot cell_type_rank A data frame containing cell type rank result.
#' @import methods
#' @name scRank
#' @rdname scRank
#' @aliases scRank-class
#' @exportClass scRank

setClass('scRank', representation(
    data = 'list',
    meta = 'list',
    para = 'list',
    net = 'list',
    cell_type_rank = 'data.frame'
))

#' @title Show scRank object
#' @param object A \code{scRank} object.
#' @export
#' @docType methods
#' 
setMethod(
  f = 'show',
  signature = 'scRank',
  definition = function(object) {
    cat("An object of class", class(object),"created from a dataset",
        "with",nrow(object@data$data), "genes and",  ncol(object@data$data), "cells. \n")
    
    cat(paste0(length(object@para$gene4use)," genes are used in networks. \n"))
    
    cat(paste0("Cells are grouped by the ","'",object@para$cell_type,"'"," column for ranking!"))
  }
)

#' Subset ScRank Object by Expression
#'
#' @param object ScRank object
#' @param expression Logical expression to subset cells, e.g. condition == 'group1', 'condition' must exist in the column of metadata
#' @param ... Additional parameters to subset


subsetScRank <- function(object, expression, ...) {
  expr_quo <- enquo(expression)  # capture without evaluating
  
  new_seuratObj <- subset(
    x = object@data$seuratObj,
    subset = !!expr_quo,  # unquote into subset()
    ...
  )
  
  new_object <- object
  new_object@data$seuratObj <- new_seuratObj
  
  cells.use <- colnames(new_seuratObj)
  
  meta <- object@meta$rawmeta
  cells.use.index <- which(rownames(meta) %in% cells.use)
  
  new_object@meta$rawmeta <- meta[cells.use.index, , drop = FALSE]
  new_object@data$data <- object@data$data[, cells.use, drop = FALSE]
  
  if(!identical(colnames(new_object@data$data),rownames(new_object@meta$rawmeta))){
    stop("Unmatched cell names in data and meta")
  }
  
  # Handle optional slots - cell type ranking
  if (length(object@cell_type_rank) > 0) {
    cells.use.group <- unique(new_object@meta$rawmeta[,new_object@para$cell_type])
    
    # Check if cell_type column exists
    if (new_object@para$cell_type %in% colnames(new_object@meta$rawmeta)) {
      # Only keep cell types that exist in both the ranking and the subset
      available_types <- intersect(cells.use.group, rownames(new_object@cell_type_rank))
      
      if (length(available_types) > 0) {
        new_object@cell_type_rank <- object@cell_type_rank[available_types, , drop = FALSE]
        
        # Handle dpGRNmeta if it exists
        if (length(object@meta$dpGRNmeta) > 0) {
          available_dpGRN <- intersect(available_types, names(object@meta$dpGRNmeta))
          if (length(available_dpGRN) > 0) {
            new_object@meta$dpGRNmeta <- object@meta$dpGRNmeta[available_dpGRN]
          } else {
            new_object@meta$dpGRNmeta <- list()
          }
        }
      } else {
        new_object@cell_type_rank <- matrix(nrow = 0, ncol = ncol(object@cell_type_rank))
        new_object@meta$dpGRNmeta <- list()
        warnings("Zero matched result!")
      }
    }
  }
  
  # Handle network data
  if (length(object@net) > 0) {
    if (new_object@para$cell_type %in% colnames(new_object@meta$rawmeta)) {
      cells.use.group <- unique(new_object@meta$rawmeta[,new_object@para$cell_type])
      
      available_net_types <- intersect(cells.use.group, names(object@net))
      
      if (length(available_net_types) > 0) {
        new_object@net <- object@net[available_net_types]
      } else {
        new_object@net <- list()
      }
    }
  }
  
  
  return(new_object)
}
