#' @title Visualize the dimensional reduction plot with scRank prioritization result
#' 
#' @description Dimensional reduction plot superimposed with the drug response hierarchy across cell types in a single-cell transcriptomic dataset
#' @param object  scRank object generated from \code{\link{rank_celltype}}
#' @param coordinate a dataframe containing the low dimensional coordinates for each cells. The row name is the cell name identical with cell name in data of scRank object.
#' @param color a vector containing continuous color used to visualize cell type rank. Default is the viridis color.
#' @param point_size point size. Default is \code{0.5}
#' @param box_padding Amount of padding around bounding box, as unit or number. Default is \code{0.5}
#' @importFrom tibble rownames_to_column
#' @importFrom stats median
#' @importFrom dplyr mutate filter arrange desc pull left_join group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom magrittr %>%
#' @import Seurat ggplot2 viridis
#' @export

plot_dim <- function(object,
                     coordinate = NULL,
                     color = NULL,
                     point_size = 0.5,
                     box_padding = 0.5) {
  # check
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }

  if (is.null(object@cell_type_rank)) {
    stop("Error: Please run 'rank_celltype()' function before plotting result!")
  }

  # process cell type rank
  perb_score <- object@cell_type_rank

  # get cell type labels and tsne/umap coordinates from Seurat
  seuratObj <- object@data$seuratObj
  if (is.null(seuratObj@reductions$tsne) & is.null(seuratObj@reductions$umap)) {
    if (is.null(coordinate)) {
      stop("UMAP or TSNE coordinates is not provided in Seurat object. Please add it in 'coordinate'")
    } else {
      # get PC data
      reduc <- coordinate
      if (identical(rownames(reduc), colnames(object@data$data))) {
        stop("Please provide correct coordinates!")
      }
    }
  } else {
    if (is.null(seuratObj@reductions$umap)) {
      reduction <- "tsne"
    }
    reduction <- "umap"
    # get PC data
    reduc <- Seurat::Embeddings(seuratObj, reduction = reduction) %>% data.frame()
  }

  celltype <- unique(object@para$ct.keep)
  # add cluster
  cells <- object@meta$rawmeta %>%
    dplyr::filter(.data[[object@para$cell_type]] %in% celltype) %>%
    rownames(.)
  labels <- object@meta$rawmeta %>%
    dplyr::filter(.data[[object@para$cell_type]] %in% celltype) %>%
    dplyr::pull(.data[[object@para$cell_type]])
  reduc <- reduc[cells, ]
  reduc <- reduc %>% mutate(cluster = labels)

  plt_data <- reduc %>% dplyr::left_join(perb_score %>%
    tibble::rownames_to_column("cluster"),
  by = "cluster"
  )
  colnames(plt_data)[1:2] <- c("Dim1", "Dim2")

  # add highlight data
  highlight_label <- perb_score %>%
    arrange(desc(fill)) %>%
    rownames(.)

  highlight <- plt_data %>%
    filter(cluster %in% highlight_label) %>%
    group_by(cluster) %>%
    summarise(
      Dim1 = median(Dim1),
      Dim2 = median(Dim2)
    )

  p <- ggplot(data = plt_data, aes(x = Dim1, y = Dim2, fill = cluster)) +
    geom_point(size = point_size, aes(color = `top_rank%`)) +
    ggrepel::geom_text_repel(
      data = highlight,
      aes(Dim1, Dim2, label = cluster),
      color = "black",
      box.padding = box_padding,
      min.segment.length = 0
    ) +
    guides(
      fill = "none",
      color = guide_colorbar(
        nbin = 20, raster = FALSE, ticks = FALSE,
        title.position = "top", title.hjust = 0.5
      )
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      panel.border = element_rect(color = "black"),
      legend.justification = "top",
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 7.5),
      legend.key.size = unit(0.25, "lines"),
      legend.key.width = unit(0.4, "lines")
    )

  if (is.null(color)) {
    p <- p +
      viridis::scale_color_viridis(
        name = "Rank (%)",
        labels = c(0, 100),
        limits = c(0, 1),
        breaks = c(0, 1),
        na.value = "white"
      )
    p
  } else {
    p <- p +
      scale_color_gradientn(
        colours = color,
        name = "Rank (%)",
        labels = c(0, 100),
        limits = c(0, 1),
        breaks = c(0, 1),
        na.value = "white"
      )
  }

  return(p)
}

#' @title Visualize the modularized subnetwork with heatmap or network
#'
#' @description Modularization and visualization of drug-related subnetwork for cell types.
#' @param object  scRank object generated from \code{\link{Constr_net}}
#' @param celltype select a cell type to visualize its network.
#' @param mode select a mode to visualize network. Default is \code{heatmap}.
#' @param gene_set a vector containing interested gene set for visualizing network.
#' @param min_ModuleSize parameter in \code{cutreeDynamic} for the minimal the number of genes in module. Defalut is \code{10}.
#' @param highlight_gene highlighted gene in network. Default is \code{target gene}
#' @param vertex_label_cex node size in network plot. Default is \code{0.5}
#' @param charge charge in layout of node. Default is \code{0.01}
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @importFrom stats dist hclust cutree
#' @importFrom igraph graph_from_adjacency_matrix layout.graphopt degree V E
#' @import Seurat ggplot2 dynamicTreeCut ggpubr ComplexHeatmap
#' @export

plot_net <- function(object,
                     celltype = NULL,
                     mode = "heatmap",
                     gene_set = NULL,
                     min_ModuleSize = 10,
                     highlight_gene = NULL,
                     vertex_label_cex = 0.5,
                     charge = 0.01) {

  # check
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }

  if (length(object@net) < 1) {
    stop("Error: Please run 'Constr_net()' function before plotting network!")
  }

  if (is.null(celltype)) {
    stop("Error: Please assign value for 'celltype'")
  }

  # check gene
  if (is.null(highlight_gene)) {
    highlight_gene <- object@para$target
  }

  if (is.null(gene_set)) {
    top_celltype <- object@cell_type_rank %>%
      arrange(desc(rank)) %>%
      rownames(.) %>%
      .[1]
    target_link_gene <- object@net[[top_celltype]][highlight_gene, object@net[[top_celltype]][highlight_gene, ] > 0] %>% names()
    gene_set <- base::union(target_link_gene, highlight_gene)
  } else {
    gene_set <- intersect(rownames(object@net[[1]]), gene_set)
    gene_set <- base::union(gene_set, highlight_gene)
  }

  # module analysis
  new_net <- object@net[[top_celltype]][gene_set, gene_set]

  hclust_out <- stats::hclust(stats::dist(1 - as.matrix(new_net)))

  gene_groups <- dynamicTreeCut::cutreeDynamic(
    hclust_out,
    minClusterSize = min_ModuleSize,
    method = "tree",
    deepSplit = FALSE,
    useMedoids = FALSE
  )

  gene_clusters <- stats::cutree(hclust_out, length(unique(gene_groups)) - 1)

  # zoom in
  order_gene <- gene_clusters %>% sort()
  new_gene <- names(order_gene[order_gene == 7])
  new_net <- new_net[new_gene, new_gene]
  hclust_out <- stats::hclust(stats::dist(1 - as.matrix(new_net)))

  gene_groups <- dynamicTreeCut::cutreeDynamic(
    hclust_out,
    minClusterSize = min_ModuleSize,
    method = "tree",
    deepSplit = FALSE,
    useMedoids = FALSE
  )

  gene_clusters <- stats::cutree(hclust_out, length(unique(gene_groups)) - 1)

  order_gene <- gene_clusters %>%
    sort() %>%
    names()

  plot_data <- object@net[[celltype]][order_gene, order_gene]

  # plot
  if (mode == "heatmap") {
    # Generate annotations for rows
    nM <- do.call(c, sapply(1:length(as.numeric(table(gene_clusters))), function(i) {
      rep(paste0("Module", i), as.numeric(table(gene_clusters))[i])
    }))
    annotation_col <- data.frame(
      Module = factor(nM)
    )
    rownames(annotation_col) <- names(sort(gene_clusters))
    mat_colors <- list(Module = ggpubr::get_palette(palette = "lancet", k = length(unique(nM))))
    names(mat_colors$Module) <- unique(nM)
    bk <- seq(from = 0, to = 1, by = 0.1)
    gaps <- sapply(unique(nM), function(x) {
      max(which(nM == x))
    }) %>%
      .[1:length(unique(nM)) - 1] %>%
      as.numeric()

    id_other <- which(rownames(plot_data) != highlight_gene)
    id_drug <- which(rownames(plot_data) == highlight_gene)
    rownames(plot_data)[id_other] <- ""
    rownames(plot_data)[id_drug] <- paste0("<", highlight_gene)

    suppressWarnings(ComplexHeatmap::pheatmap(abs(plot_data),
      cluster_cols = F, cluster_rows = F,
      show_rownames = T, show_colnames = F,
      breaks = bk,
      border_color = NA,
      name = "abs(cor)",
      annotation_col = annotation_col,
      annotation_colors = mat_colors,
      annotation_names_col = F,
      annotation_legend = T,
      gaps_row = gaps,
      gaps_col = gaps,
      number_color = "black",
      legend = F
    ))
  } else if (mode == "network") {
    Y <- plot_data
    Y <- as.matrix(Y)
    module_number <- as.numeric(gene_clusters[which(names(gene_clusters %>% sort()) == highlight_gene)])
    module_show <- gene_clusters[gene_clusters == module_number] %>% names()
    Y <- igraph::graph_from_adjacency_matrix(Y[module_show, module_show], weighted = T, diag = FALSE)

    gY <- rownames(Y[])
    set.seed(12)

    lNet <- igraph::layout.graphopt(Y, charge = charge)

    color_use <- ggpubr::get_palette(palette = "lancet", k = length(unique(gene_clusters)))[module_number]

    fcolor <- ifelse(gY %in% highlight_gene, "forestgreen", color_use)

    Y_deg <- igraph::degree(Y, mode = c("ALL"))
    igraph::V(Y)$degree <- Y_deg

    plot(Y,
      layout = lNet, edge.arrow.size = 0, # edge.size = abs(E(Y)$weight),
      edge.color = ifelse(igraph::E(Y)$weight > 0, "gray60", "#6565FF"), vertex.label = gY, vertex.size = 10, vertex.label.dist = 1,
      vertex.label.family = "Arial", vertex.label.color = "black", vertex.label.cex = vertex_label_cex, edge.width = abs(igraph::E(Y)$weight),
      vertex.frame.width = 1, vertex.color = fcolor, vertex.frame.color = fcolor
    )
  } else {
    stop("Please provide either 'heatmap' or 'network' in 'mode'")
  }
}
