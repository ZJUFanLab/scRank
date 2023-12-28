#' @title Visualize the dimensional reduction plot with scRank prioritization result
#' 
#' @description Dimensional reduction plot superimposed with the drug response hierarchy across cell types in a single-cell transcriptomic dataset
#' @param object  scRank object generated from \code{\link{rank_celltype}}
#' @param reductions Reduction name in the Seurat object. If it is NULL then present UMAP or TSNE, or input customized reduction name.
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
                     reductions = NULL,
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
    reduction_name <- NULL
    if ('umap' %in% names(seuratObj@reductions)) {
      reduction_name <- 'umap'
    } else if ('tsne' %in% names(seuratObj@reductions)) {
      reduction_name <- 'tsne'
    }
    # get PC data
    if(is.null(reductions)){
      if(is.null(reduction_name)){
        stop("Error! It seems neither `umap` or `tsne` in the reduction of seurat objct. Please specify the reduction name in `reductions`")
      } else{
        reduction = reduction_name
      }
    } else{
      reduction = reductions
    }
    reduc <- Seurat::Embeddings(seuratObj, reduction = reduction) %>% data.frame()
  }

  celltype_ <- unique(object@para$ct.keep)
  # add cluster
  cells <- object@meta$rawmeta %>%
    dplyr::filter(.data[[object@para$cell_type]] %in% celltype_) %>%
    rownames(.)
  labels <- object@meta$rawmeta %>%
    dplyr::filter(.data[[object@para$cell_type]] %in% celltype_) %>%
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

#' @title Initialization for modularizing subnetwork 
#'
#' @description Modularization of drug-related subnetwork for cell types.
#' @param object  scRank object generated from \code{\link{Constr_net}}
#' @param min_ModuleSize parameter in \code{cutreeDynamic} for the minimal the number of genes in module. Default is \code{10}.
#' @param customized_gene_set a vector contained gene name. Default is \code{Null}, but can also input customized gene set like disease-assoiciated genes.
<<<<<<< HEAD
#' @param zoom a logical value for zooming in target-related module. Default is \code{TRUE}
=======
#' @param perturbed_target target name of drug perturbation. Default is in scRank object.
#' @param zoom a logical value for zooming in target-related module. Default is \code{TRUE}
#' @param top_celltype name of cell type name. Default is in scRank object.
#' @return scRank object with moularized network
>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @importFrom stats dist hclust cutree
#' @import dynamicTreeCut ggpubr ComplexHeatmap
#' @export

init_mod <- function(object,
                     min_ModuleSize = 10, 
                     customized_gene_set = NULL,
<<<<<<< HEAD
                     zoom = T){
  
  
  # load gene set
  
  top_celltype <- object@cell_type_rank %>%
    arrange(desc(rank)) %>%
    rownames(.) %>%
    .[1]
  
  if(is.null(customized_gene_set)){
    target_link_gene <- object@net[[top_celltype]][object@para$target, object@net[[top_celltype]][object@para$target, ] > 0] %>% names()
    gene_set <- base::union(target_link_gene, object@para$target)
  } else{
    customized_gene_set <- base::intersect(customized_gene_set, object@para$gene4use)
    gene_set <- base::union(customized_gene_set, object@para$target)
=======
                     perturbed_target = NULL,
                     zoom = T,
                     top_celltype = NULL){
  
  # check
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }
  if (is.null(object@cell_type_rank)) {
    stop("Error: Please run 'rank_celltype()' function before plotting result!")
  }
  if (is.null(perturbed_target)) {
    perturbed_target = object@para$target
  }
  
  # load gene set
  if(is.null(top_celltype)){
    top_celltype <- object@cell_type_rank %>%
      arrange(desc(rank)) %>%
      rownames(.) %>%
      .[1]
  }
  
  if(is.null(customized_gene_set)){
    target_link_gene <- lapply(perturbed_target,function(x){
      tmp <- colnames(object@net[[top_celltype]])[abs(object@net[[top_celltype]][x,]) > 0]
      c(tmp,x)
    })
    
    if(length(perturbed_target) > 1){
      gene_set <- do.call("union",target_link_gene)
    } else{
      gene_set <- unlist(target_link_gene)
    }
    
  } else{
    customized_gene_set <- base::intersect(customized_gene_set, object@para$gene4use)
    gene_set <- base::union(customized_gene_set, perturbed_target)
>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
  }
  
  # initialize modularization
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
  gene_num = length(gene_clusters)
  if (zoom & gene_num > 100) {
    order_gene <- gene_clusters %>% sort()
<<<<<<< HEAD
    cent_loc <- order_gene[which(names(order_gene) == object@para$target)] %>% as.numeric()
    new_gene <- names(order_gene[order_gene == cent_loc])
=======
    cent_loc <- order_gene[which(names(order_gene) %in% perturbed_target)] %>% as.numeric()
    new_gene <- names(order_gene[order_gene %in% cent_loc])
>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
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
    order_gene <- gene_clusters %>% sort() %>% names()
  } else{
    order_gene <- gene_clusters %>% sort() %>% names()
  }
  
  object@para$gene_clusters = gene_clusters
  object@para$order_gene = order_gene
  object@para$subnet_node = gene_set
  
  return(object)
}


#' @title Visualize the modularized subnetwork with heatmap or network
#'
#' @description visualization of target-related subnetwork for cell types.
#' @param object  scRank object generated from \code{\link{Constr_net}}
#' @param celltype select a cell type to visualize its network.
#' @param mode select a mode to visualize network by \code{network} or \code{heatmap}. Default is \code{heatmap}.
#' @param gene_set a vector containing interested gene set for visualizing network.
#' @param highlight_gene highlighted gene in network. Default is \code{target gene}
#' @param vertex_label_cex node size in network plot. Default is \code{0.5}
#' @param charge charge in layout of node. Default is \code{0.01}
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @importFrom stats dist hclust cutree
#' @importFrom igraph graph_from_adjacency_matrix layout.graphopt degree V E
#' @import dynamicTreeCut ggpubr ComplexHeatmap
#' @export

plot_net <- function(object,
                     celltype = NULL,
                     mode = "heatmap",
                     vertex_label_cex = 0.5,
                     charge = 0.01,
                     highlight_gene = NULL,
                     gene_set = NULL) {

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
<<<<<<< HEAD

  highlight_gene <- object@para$target

  gene_set <- object@para$subnet_node
  
=======
  if (is.null(highlight_gene)) {
    highlight_gene <- object@para$target
  }
  
  if (is.null(gene_set)) {
    gene_set <- object@para$subnet_node
  }

>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
  order_gene = object@para$order_gene
  gene_clusters = object@para$gene_clusters

  plot_data <- object@net[[celltype]][order_gene, order_gene]

  # plot
  if (mode == "heatmap") {
    # Generate annotations for rows
    if(length(unique(gene_clusters)) < 2){
      nM <- rep(paste0("Module1"), length(gene_clusters))
    } else {
      nM <- do.call(c, sapply(1:length(as.numeric(table(gene_clusters))), function(i) {
        rep(paste0("Module", i), as.numeric(table(gene_clusters))[i])
      }))
    }
    
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

    id_other <- which(!rownames(plot_data) %in% highlight_gene)
    id_drug <- which(rownames(plot_data) %in% highlight_gene)
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
    plot_data <- object@net[[celltype]][order_gene, order_gene]
    Y <- plot_data
    Y <- as.matrix(Y)
    module_number <- as.numeric(gene_clusters[highlight_gene])
<<<<<<< HEAD
    module_show <- gene_clusters[gene_clusters == module_number] %>% names()
=======
    module_show <- gene_clusters[gene_clusters %in% module_number] %>% names()
>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
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

#' @title Gene set enrichement analysis for determing the disease relevance of cell type predicted by scRank.
#'
#' @description Use the ranked genes in the cell type and the significantly up-regulated genes in disease condition compared to the normal condition to perform the GSEA analysis. The final p value would reflect the significance of the disease relevance for the cell type.
#' @param object  scRank object generated from \code{\link{Constr_net}}
#' @param celltype select a cell type to measure its relevance to the disease.
#' @param disease_gene a vector containing the positive differential expressed gene names between health and disease state.
#' @param selectGenes a character vector containing the selected genes shown in the plot. default is the top 10 genes.
#' @return GSEA plotting with p value indicating the significance of the disease relevance
#' @importFrom fgsea fgsea
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @importFrom cowplot plot_grid
#' @importFrom ggrepel geom_text_repel
#' @export
scRank_GSEA <- function(object,
                        celltype = NULL,
                        disease_gene = NULL,
                        selectGenes = NULL
) {
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }
  
  if (is.null(celltype)) {
    stop("Error: Please assign value for 'celltype'")
  } 
  
  if (is.null(disease_gene)) {
    stop("Please provide the disease differential experssed genes in `disease_gene`")
  } else {
    # check disease gene
    if(is.vector(disease_gene)){
      disease_gene <- as.character(disease_gene)
    } else {
      stop("Please provide the disease differential experssed genes in `disease_gene`")
    }
  }
  
  if (is.null(object@meta$dpGRNmeta)){
    stop("Error: Please run 'rank_celltype()' function before!")
  }
  
  if (is.null(selectGenes)) {
    selectGenes <- disease_gene[1:10]
  } else {
    if (is.vector(selectGenes)) {
      selectGenes <- as.character(selectGenes)
      selectGenes <- selectGenes[selectGenes %in% disease_gene & selectGenes %in% rownames(object@net[[celltype]])]
    } else {
      stop("Please provide the selected genes in `selectGenes`")
    }
  }
  
  # get drug perturbed distance from celltype 
  df <- object@meta$dpGRNmeta[[celltype]]
  distance <- df %>% arrange(desc(distance)) %>% mutate(z = -1/log2(distance)) %>%
    pull(z, gene) %>% 
    scale() %>%
    .[,1]
  
  disease_gene <- intersect(disease_gene, df$gene)
  term2gene <- list("Disease DEG" = disease_gene)
  
  gsea_res <- fgsea(pathways = term2gene, 
                    stats    = distance,
                    eps      = 1e-10,
                    minSize  = 0,
                    maxSize  = 500)
  
  print(gsea_res)
  
  # get gsea result
  pval <- gsea_res$pval
  x <- gsea_res
  if(pval < 0.05) pcol <- "darkred" else pcol <- "darkblue"
  if(pval < 0.00001) pval_label <- '< 0.00001' else pval_label = paste0("= ",round(pval,5))
  
  
  geneList <- position <- NULL
  geneList <- distance
  geneSetID <- x$pathway
  geneSet <- disease_gene
  exponent <- 1
  df <- .gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$gene <- names(geneList)
  df$Description <- x$pathway
  
  # plot the ES 
  top_p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = runningScore, color = Description), size = 1) +
    geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
    scale_color_manual(values = pcol, labels = paste0("Disease DEG for ",celltype,"      ","(p.value ", pval_label,")")) +
    theme_bw() +
    ylab("Enrichment Score") +
    theme(
      panel.grid = element_blank(),
      legend.position = 'top',
      legend.title = element_blank(),
      axis.text.y = element_text(size = 12,family = 'Arial',color = 'black'),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
    )   
  
  # plot distribution of target gene set
  mid_p <- ggplot(df, aes(x = x)) +
    geom_linerange(aes(ymin = ymin, ymax = ymax, color = Description)) +
    scale_color_manual(values = 'black') +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = 'none',
      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = margin(t=-.1, b=0,unit="cm")
    )
  
  # define color cut range 
  v <- seq(1, sum(df$position),length.out=9)
  inv <- findInterval(rev(cumsum(df$position)),v)
  if( min(inv) == 0 ) inv <- inv + 1
  col <- c("#08519C", "#3182BD", "#6BAED6", "#BDD7E7", "#EFF3FF", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26")
  
  ymin <- min(mid_p$data$ymax)
  yy <- max(mid_p$data$ymax - mid_p$data$ymin) * .3
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  d <- data.frame(ymin = ymin, ymax = yy,
                  xmin = xmin,
                  xmax = xmax,
                  col = col[unique(inv)])
  
  col_p <- ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill = I(col) )) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t=0, b=0,unit="cm")
    ) 
  # plot selected gene name
  low_p <- ggplot(df %>% filter(!x %in% 1), aes(x, geneList, label = gene)) +
    geom_bar(position = "dodge", stat = "identity", color = "gray90", width=0.5) +
    geom_point(data = df %>% filter(gene %in% selectGenes), size = 2, color = 'white', fill = 'darkred', shape = 21) +
    ggrepel::geom_text_repel(data = df %>% filter(gene %in% selectGenes),
                             box.padding = unit(0.35, "lines"),
                             max.overlaps = Inf,
                             direction = "y",
                             angle = 90,
                             size = 2.5,
                             point.padding = unit(0.3, "lines")) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
    ylab("Scaled distance") +
    xlab("Ranked in ordered gene list") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm")
    )
  
  # combined plot
  rel_heights <- c(1.5, .5, .2, 1.5)
  plotlist <- list(top_p, mid_p, col_p, low_p)
  
  n <- 4
  p <- cowplot::plot_grid(plotlist = plotlist, ncol = 1, align="v", axis = 'tblr', rel_heights = rel_heights)
  return(p)
}

#' @title Plotting drug function on cell type
#' @description Plotting drug function on cell type based on the drug perturbed gene set.
#' @param object scRank object generated from \code{\link{Constr_net}}
#' @param celltype select a cell type to measure the drug curative efficacy.
#' @param category the category of gene set in MSigDB. Default is "H" representing the hallmark gene sets.
#' @param top_number the number of top pathways to be shown in the plot. Default is 5.
#' @param show_leading_edge whether to show the leading edge genes in the plot. Default is TRUE
#' @import msigdbr ggplot2 
#' @importFrom fgsea fgsea
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @importFrom crayon cyan
#' @importFrom scales rescale
#' @importFrom ggridges geom_density_ridges position_points_jitter theme_ridges
#' @importFrom tibble as_tibble

plot_drug_function <- function(object, 
                               celltype = NULL,
                               category = "H",
                               top_number = 5,
                               show_leading_edge = T){
  # check
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }
  if (is.null(object@meta$dpGRNmeta)) {
    stop("Error: Please run 'Constr_net()' function before plotting result!")
  }
  if(is.null(celltype) | !(celltype %in% object@para$ct.keep)) {
    stop("Error: Please input correct name of cell type with network!")
  } else{
    df <- object@meta$dpGRNmeta[[celltype]]
  }
  
  message(crayon::cyan("Initializing ..."))
  
  gseaParam <- 1
  # scale distance  
  distance <- df  %>% arrange(desc(distance)) %>% mutate(z = -1/log2(distance)) %>% 
    pull(z,gene) %>%
    scale() %>% .[,1]
  
  # prepare gene id 
  gene2id <- utile_database$Gene_Info[[object@para$species]] %>% filter(Symbol %in% names(distance)) %>% pull(GeneID,Symbol)
  gene2id <- gene2id[names(distance)] %>% na.omit()
  generanks <- distance[names(gene2id)]
  names(generanks) <- gene2id %>% as.vector()
  
  # prepare gene set
  if(object@para$species == 'human'){
    species = "Homo sapiens"
  } else if (object@para$species == 'mouse'){
    species = "Mus musculus"
  }
  m_df <- msigdbr(species = species, category = category)
  fgsea_sets <- m_df %>% split(x = .$entrez_gene, f = .$gs_name)
  
  message(crayon::cyan("Running ..."))
  # run GSEA
  fgseaRes <- fgsea(pathways = fgsea_sets,stats = generanks)
  topPathwaysUp <- fgseaRes %>% as_tibble() %>% filter(ES > 0) %>% top_n(top_number,-pval) %>% pull(pathway)
  
  message(crayon::cyan("Plotting ..."))
  # plotting
  pathways <- fgsea_sets[topPathwaysUp]
  
  rnk <- rank(-generanks)
  ord <- order(rnk)
  statsAdj <- generanks[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  
  pathways <- lapply(pathways, function(p) {
    p_new <- unname(as.vector(na.omit(match(p, names(statsAdj)))))
    names(p_new) <- p[!is.na(match(p, names(statsAdj)))]
    return(p_new)
  })
  pathways <- pathways[sapply(pathways, length) > 0]
  
  
  pathway.dataframes <- do.call('rbind',lapply(names(pathways), function(name){
    fgseaRes_tibble <- fgseaRes %>% as_tibble() %>% filter(pathway %in% name)
    df1 <- data.frame(Pathway = name, GeneID = names(pathways[[name]]), value = generanks[names(pathways[[name]])], statsAdj = unname(pathways[[name]]), rank = rnk[names(pathways[[name]])],
                      pval = fgseaRes_tibble$pval,
                      padj = fgseaRes_tibble$padj,
                      log2err = fgseaRes_tibble$log2err,
                      ES = fgseaRes_tibble$ES,
                      NES = fgseaRes_tibble$NES,
                      size = fgseaRes_tibble$size,
                      leadingEdge = fgseaRes_tibble$leadingEdge %>% unlist() %>% str_c(collapse = ","))
    df2 <- utile_database[["Gene_Info"]][[object@para$species]] %>% filter(GeneID %in% df1$GeneID) %>% distinct(GeneID,.keep_all = T) %>% select(GeneID, Symbol)
    df3 <- merge(df1,df2,by = 'GeneID') %>% select(GeneID,Symbol,Pathway,value,rank,statsAdj,pval,padj,NES,size,ES,log2err,leadingEdge)
  }))
  
  if(show_leading_edge){
    pathway.dataframes <- pathway.dataframes %>% group_by(Pathway) %>% filter(GeneID %in% unlist(strsplit(leadingEdge, split = ","))) %>% ungroup()
  }
  
  print(pathway.dataframes %>% distinct(Pathway,.keep_all = T))
  
  p.ridges <- ggplot(pathway.dataframes, aes(x = value, y = Pathway,  fill = -log10(padj), color = -log10(padj))) +
    geom_density_ridges(
      jittered_points = TRUE, scale = .95, rel_min_height = .01,
      point_shape = "|", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0),
      alpha = 0.5
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), name = "Scaled Distance") +
    scale_fill_gradient(low = "grey50",high = "#E4AF03")	+
    scale_color_gradient(low = "grey50",high = "#E4AF03") +
    coord_cartesian(clip = "off") +
    theme_ridges(center_axis_labels = TRUE) 
  
  return(p.ridges)
}
