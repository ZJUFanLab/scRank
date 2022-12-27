#' @title scRank object
#'
#' @description create scRank object using single-cell transcriptomics data.
#' @param input gene expression profile formatted by matrix, data frame, or \code{Seurat}
#' @param meta  meta data describing the cells in input and its corresponding cell type; can be left as \code{NULL} if \code{input} is a \code{Seurat} object
#' @param cell_type the name of column which containing cell type labels for
#' each cell corresponding to gene-by-cell expression matrix
#' @param species characters meaning species of the single-cell transcriptomic data. \code{human} or \code{mouse}.
#' @param drug characters meaning the name of inhibitor drug.
#' @param target characters meaning the gene encoded target of interested drug.
#' @import Seurat
#' @importFrom tester is_numeric_matrix is_numeric_dataframe
#' @importFrom methods as
#' @importFrom dplyr select filter pull top_n
#' @importFrom magrittr %>%
#' @return scRank object
#' @export

CreateScRank <- function(input,
                         meta = NULL,
                         cell_type = "celltype",
                         species,
                         drug = NULL,
                         target = NULL) {
  # import data and meta from Seurat
  if (is(input, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Please confirm \"Seurat\" R package is installed correctly")
    }
    seuratObj <- input
    meta <- input@meta.data %>% droplevels()
    cell_types <- meta[[cell_type]]
    data <- as.matrix(Seurat::GetAssayData(
      object = input, slot = "counts",
      assay = "RNA"
    ))
  } else {
    # check meta
    if (is.null(meta)) {
      stop("Please provide a meta if not input a Seurat object!")
    }
    if (all(is.data.frame(meta), colnames(meta) %in% cell_type)) {
      stop("Please provide a correct meta formatted with data.frame!") # add demo_meta()
    }
    seuratObj <- NULL
    # cell_types <- meta %>%
    #   droplevels() %>%
    #   dplyr::select(cell_type) %>%
    #   as.character()

    # check input
    valid_format <- is(input, "dgCMatrix ") || tester::is_numeric_matrix(input) ||
      tester::is_numeric_dataframe(input)

    if (valid_format) {
      data <- methods::as(input, "dgCMatrix")
    }

    # check cells are matched in data and meta
    if (!identical(ncol(data), nrow(meta))) {
      stop(
        "The cells in meta and data are mismatched, where meta with ",
        nrow(meta), " cells and data with ", ncol(data), " cells."
      )
    }
    if (!identical(colnames(data), rownames(meta))) {
      stop("colnames(data) is not consistent with rownames(meta)")
    }
  }

  # Feature selecting

  # check HVG
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Please confirm \"Seurat\" R package is installed correctly",
      .call = FALSE
    )
  }
  if (is(input, "Seurat")) {
    if (nrow(input) > 2000) {
      hvg <- .hvg(input)
    } else {
      hvg <- .hvg(input, nfeatures = nrow(input))
    }
  } else {
    seuratObj <- CreateSeuratObject(data, meta.data = meta)
    if (nrow(seuratObj) > 2000) {
      hvg <- .hvg(seuratObj)
    } else {
      hvg <- .hvg(seuratObj, nfeatures = nrow(seuratObj))
    }
  }

  # check drug gene
  if (species == "human") {
    drug_gene <- utile_database$Drug_Target$human$Symbol
  }
  if (species == "mouse") {
    drug_gene <- utile_database$Drug_Target$mouse$mousegene
  }

  # check TF gene
  if (species == "human") {
    tf_gene <- utile_database$Gene_TF$human$Symbol
  }
  if (species == "mouse") {
    tf_gene <- utile_database$Gene_TF$mouse$Symbol
  }

  # check target gene
  if (is.null(drug)) {
    if (is.null(target)) {
      stop("Please input the target gene of inhibitor or inhibitor name")
    } else {
      if (target %in% rownames(data)) {
        target_gene <- target
      } else {
        stop("Please check if the target gene is in the gene expression profile. ")
      }
    }
  } else {
    drug_name <- toupper(drug)
    if (drug_name %in% utile_database$Drug_Info$drug_name) {
      target <- utile_database$Drug_Info %>%
        dplyr::filter(drug %in% drug_name) %>%
        dplyr::top_n(1, score) %>%
        dplyr::pull(gene_name)
    } else {
      stop("Drug name is not in our list, please assign the known inhibited
           target gene name in to CreateScRank(target=`known target`)")
    }
  }

  gene4use <- unique(c(target, hvg, tf_gene, drug_gene))

  # filtering mitochondrial and ribosomal gene
  gene4use <- gene4use[-grep(pattern = "^RP[[:digit:]]+|^RPL|^RPS|^MT-", toupper(gene4use))]

  # check if gene is detected in expression profile
  gene4use <- gene4use[gene4use %in% rownames(data)]

  # create scRank object
  object <- new("scRank",
    data = list(data = data, seuratObj = seuratObj),
    meta = list(rawmeta = meta),
    para = list(
      species = species,
      cell_type = cell_type,
      drug = drug,
      target = target,
      gene4use = gene4use
    )
  )

  return(object)
}

#' @title Constructing cell-type-specific GRN for single-cell transcriptomics data
#'
#' @description construct gene-gene co-expression adjacency matrix for every cell type with regression strategy.
#' @param object scRank object generated from \code{\link{CreateScRank}}.
#' @param min_cells the minimum number of cells in a particular cell type retained for constructing robust scGRN. Cell type with cells lower than \code{min_cells} will be droped. Default is 25.
#' @param select_ratio percent cells of total cells to be selected in random selection. Default is 0.5.
#' @param n_selection number of selection to create \code{n_selection} networks for a cell type. Default is 10.
#' @param cut_ratio threshold of top percent edge weight to be cut. Edge with weight lower than the threshold will be cut. Default is 0.95.
#' @param n.core Number of CPU cores to use. Default is a half of all cores.
#' @return scRank object containing the gene regulatory network for every cell type.
#' @import pbapply doParallel parallel foreach RSpectra methods
#' @importFrom crayon cyan
#' @importFrom Matrix t
#' @importFrom rTensor as.tensor cp
#' @importFrom stats quantile
#' @export

Constr_net <- function(object,
                       min_cells = 25,
                       select_ratio = 0.5,
                       n_selection = 10,
                       cut_ratio = 0.95,
                       n.core = NULL) {
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }

  # clear register
  if (is.null(n.core) || 1 < n.core) {
    # clear register
    doParallel::stopImplicitCluster()
  }
  # register parallel
  if (is.null(n.core)) {
    n.core <- round(parallel::detectCores() / 2)
  }
  n.core <- max(1, n.core)
  if (n.core != 1) {
    doParallel::registerDoParallel(cores = n.core)
  }

  cell_types <- object@meta$rawmeta[[object@para$cell_type]]

  ct.keep <- names(
    table(cell_types)[table(cell_types) > min_cells]
  )

  message(
    "--- Cell population ---\n", paste0(unique(ct.keep), "\n"),
    "--- --- --- --- --- ---\nwill be kept for constructing network!"
  )

  # Initialize Net for every cell type
  Net_set <- list()

  for (ct in ct.keep) {
    message(crayon::cyan(paste0("Constructing network for ", ct)))
    ct.data <- .get_ct_data(object, ct)
    ct.data <- as.matrix(ct.data)

    # number of random selected cells
    if (ncol(ct.data) > 100) {
      cells <- select_ratio * ncol(ct.data)
    } else {
      cells <- 25
    }

    # make networks from n selected cells
    set.seed(1)
    Net_set[[ct]] <- pbapply::pbsapply(seq_len(n_selection), function(n) {
      rdm_selected_cells <- sample(
        x = seq_len(ncol(ct.data)),
        size = cells,
        replace = TRUE
      )
      Mtx <- ct.data[, rdm_selected_cells]
      Mtx <- Mtx[rowSums(Mtx) > 0, ]
      ct.net <- .make_net(mat = Mtx, cut_ratio = cut_ratio, n.cores = n.core)
      genes_all <- rownames(ct.data)
      ct.net.final <- matrix(0, length(genes_all), length(genes_all))
      rownames(ct.net.final) <- colnames(ct.net.final) <- genes_all
      ct.net.final[rownames(ct.net), colnames(ct.net)] <- as.matrix(ct.net)
      ct.net.final <- as(ct.net.final, "dgCMatrix")
      return(ct.net.final)
    })
  }

  # integrate network using CANDECOMP/PARAFRAC Tensor Decomposition
  cat(crayon::cyan("Integrating sets of network ... It might take minutes to hours."))
  if (n.core > 1) {
    Net_integrat <- foreach::foreach(ct = ct.keep, .packages = c("rTensor", "methods")) %dopar% {
      Net_set_ct <- Net_set[[ct]]
      set.seed(1)
      O <- .integrat_net(Net_set_ct)
      return(O)
    }
  } else {
    Net_integrat <- pbapply::pblapply(ct.keep, function(ct) {
      set.seed(1)
      O <- .integrat_net(Net_set[[ct]])
    })
  }

  # process bi-direction
  Net_final <- lapply(Net_integrat, function(Net) {
    Net <- .process_biedge(Net)
  })

  names(Net_final) <- ct.keep

  object@net <- Net_final
  object@para$min_cells <- min_cells
  object@para$select_ratio <- select_ratio
  object@para$n_selection <- n_selection
  object@para$cut_ratio <- cut_ratio
  object@para$n.core <- n.core
  object@para$ct.keep <- ct.keep
  return(object)
}

#' @title Prioritizing cell type responsive to drug perturbation.
#'
#' @description Rank drug-responsive cell type based on perturbation score of drug using manifold alignment and network propagation.
#' @param object scRank object generated from \code{\link{CreateScRank}}.
#' @param n.core Number of CPU cores to use. Default is a half of all cores.
#' @param perturbed_target Target to be perturbed in dpGRN for ranking cell type. Default is the target in \code{object@para$target}
#' @param n_dim number of low dimensions in manifold space using manifold alignment. Default is 2.
#' @param n_hop the number of hop for network propagation. Default is 2.
#' @return scRank object containing the drug-responsive cell type rank.
#' @import dplyr doParallel parallel foreach methods RSpectra RhpcBLASctl
#' @importFrom magrittr %>%
#' @importFrom scales rescale
#' @importFrom stats dist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tidyr replace_na
#' @export

rank_celltype <- function(object,
                          n.core = NULL,
                          n_dim = 2,
                          n_hop = 2,
                          perturbed_target = NULL) {
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }

  if (!all(object@para$target %in% rownames(object@net[[1]])) || !all(perturbed_target %in% rownames(object@net[[1]]))) {
    stop("Drug target gene is not in the network, please add it into `CreateScRank(target = _target_)'")
  }

  # clear register
  if (is.null(n.core) || 1 < n.core) {
    # clear register
    doParallel::stopImplicitCluster()
  }
  # register parallel
  if (is.null(n.core)) {
    n.core <- round(parallel::detectCores() / 2)
  }
  n.core <- max(1, n.core)
  if (n.core != 1) {
    doParallel::registerDoParallel(cores = n.core)
  }

  if (is.null(perturbed_target)) {
    perturbed_target <- object@para$target
  }
  if (length(perturbed_target) > 2) {
    stop("The number of perturbed gene should not be more than two!")
  }

  # determing hop number
  if (n_hop == 2) {
    simple_sum <- F
    multi_layer <- F
  }
  if (n_hop == 1) {
    simple_sum <- T
    multi_layer <- F
  }
  if (n_hop == 3) {
    simple_sum <- F
    multi_layer <- T
  }
  if (n_hop > 3 | n_hop < 1) {
    stop("Error: The number of n_hop should not be larger than 3")
  }

  Net <- object@net
  ct.keep <- object@para$ct.keep
  # calculating global effect based on the distance between dpGRN and GRN
  if (n.core > 1) {
    lowDim <- foreach::foreach(ct = ct.keep, .packages = c("RhpcBLASctl", "RSpectra")) %dopar% {
      GRN <- Net[[ct]]
      dpGRN <- GRN
      dpGRN[perturbed_target, ] <- 0
      lowDim_net <- .align_net(GRN, dpGRN, ndim = n_dim, n.cores = n.core)
      return(lowDim_net)
    }
  } else {
    lowDim <- lapply(ct.keep, function(ct) {
      GRN <- Net[[ct]]
      dpGRN <- GRN
      dpGRN[perturbed_target, ] <- 0
      lowDim_net <- .align_net(GRN, dpGRN, ndim = n_dim, n.cores = n.core)
      return(lowDim_net)
    })
  }
  names(lowDim) <- ct.keep

  Dist <- lapply(ct.keep, function(ct) {
    lowDim_ct <- lowDim[[ct]]
    n <- nrow(lowDim_ct) / 2
    D <- sapply(seq_len(n), function(x) {
      dims_x <- lowDim_ct[x, ]
      dims_y <- lowDim_ct[(x + n), ]
      dist <- as.numeric(stats::dist(rbind(dims_x, dims_y)))
      return(dist)
    })
    O <- data.frame(
      gene = gsub("GRN_", "", rownames(lowDim_ct)[1:n]),
      distance = D
    )
    return(O)
  })

  names(Dist) <- ct.keep

  # adding network feature
  df <- .create_final_df(Dist, object@net, perturbed_target)
  
  # calculating perturbation score
  perb_score <- .calculate_score(df, ko = perturbed_target, Net = Net, simple_sum = simple_sum, multi_layer = multi_layer)

  perb_score <- as.data.frame(perb_score)
  perb_score <- perb_score %>% mutate(
    rank = rank(perb_score),
    `top_rank%` = rank / dplyr::n(),
    `top_rank%` = scales::rescale(`top_rank%`, c(0, 1)),
    fill = `top_rank%`
  )
  object@meta$dpGRNmeta <- df
  object@cell_type_rank <- perb_score
  return(object)
}
