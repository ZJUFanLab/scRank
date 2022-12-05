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
