#' @title scRank object
#'
#' @description create scRank object using single-cell transcriptomic data.
#' @param input gene expression profile formatted by matrix, data frame, or \code{Seurat}
#' @param meta  meta data describing the cells in input and its corresponding cell type; can be left as \code{NULL} if \code{input} is a \code{Seurat} object
#' @param cell_type the name of column which containing cell type labels for
#' each cell corresponding to gene-by-cell expression matrix
#' @param species characters meaning species of the single-cell transcriptomic data. \code{human} or \code{mouse}.
#' @param drug characters meaning the name of drug.
#' @param target characters meaning the gene encoded direct target of interested drug.
#' @param type characters meaning the MOAs of drug including antagonist or agonist. Default is antagonist.
#' @param if_cluster A logical meaning whether clustering single-cell transcriptomic data. Default is \code{FALSE}.
#' @param var.genes var.genes (optional) vector of gene names representing the subset of genes used for clustering.
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
                         target = NULL,
                         type = 'antagonist',
                         if_cluster = F,
                         var.gene = NULL) {
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
    if(all(is.null(meta),if_cluster)){
      message(crayon::cyan("Note: Metadata will be created after clustering"))
    } else if (is.null(meta)) {
      stop("Please provide a meta if not input a Seurat object!")
    } else if (all(is.data.frame(meta), colnames(meta) %in% cell_type)) {
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
    } else{
      data <- input
    }
    
    # create meta using scSHC to create metadata
    if(all(is.null(meta),if_cluster)){
      if (requireNamespace("scSHC", quietly = TRUE)) {
        message(crayon::cyan("Running clustering by scSHC..."))
        if (is.null(var.genes)){
          clusters <- scSHC::scSHC(data)
        } else {
          commone_features <- intersect(rownames(data),var.genes)
          clusters <- scSHC::scSHC(data[commone_features,],num_features = length(commone_features))
        }
        cells <- names(clusters[[1]])
        labels <- unname(clusters[[1]])
        meta <- data.frame(celltype = labels)
        rownames(meta) <- cells
      } else {
        stop("scSHC is required for this functionality. Please install it using devtools::install_github(\"igrabski/sc-SHC\").")
      }
    }
    
    # check cells are matched in data and meta
    if (!identical(ncol(data), nrow(meta))) {
      stop(
        paste0("The cells in meta and data are mismatched, where meta with ",
        nrow(meta), " cells and data with ", ncol(data), " cells."
      ))
    }
    if (!identical(colnames(data), rownames(meta))) {
      stop("colnames(data) is not consistent with rownames(meta)")
    }
  }
  # check drug MOA
  if(is.null(type)){
    stop("Error! Please provide correct type of drug (antagonist or agnoist).")
  } else if (! type %in% c('antagonist','agonist')){
    stop("Error! Please provide correct type of drug (antagonist or agnoist).")
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
      if (all(target %in% rownames(data))) {
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
      gene4use = gene4use,
      type = type
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
#' @param keep_ratio percent of the weaker edge weight should be kept. Default is 0.25. 
#' @param n.core Number of CPU cores to use for constructing network. Default is all cores - 1.
#' @param n.core_cp Number of CPU cores to use for integrating network. Default is 4. Recommendation is the number lower than the number of cell type.
#' @param use_py Whether to use python to conduct CP decomposition
#' @param env When using python in \code{use_py}, please define the python environment of installed tensorly and hdf5. Default is the 'base' environment. Anaconda is recommended.
#' @return scRank object containing the gene regulatory network for every cell type.
#' @import pbapply doParallel parallel foreach RSpectra methods readr
#' @importFrom crayon cyan
#' @importFrom Matrix t
#' @importFrom rTensor as.tensor cp
#' @importFrom stats quantile
#' @export

Constr_net <- function(object,
                       select_ratio = 0.5,
                       n_selection = 10,
                       cut_ratio = 0.95,
                       keep_ratio = 0.25,
                       min_cells = 25,
                       n.core = NULL,
                       n.core_cp = 4,
                       use_py = F,
                       env = 'base') {
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
    n.core <- parallel::detectCores() - 1
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
      cells <- min_cells
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
  cat(crayon::cyan("Integrating sets of network ... It might take minutes to hours.\n"))
  if(use_py){
    cat(crayon::cyan("Using tensorly, please install the tensorly and h5py (python package) in environment of",env ,"!", "\n"))
    # clear register
    if (is.null(n.core_cp) || 1 < n.core_cp) {
      # clear register
      doParallel::stopImplicitCluster()
    }
    # python
    require("rhdf5")
    require("reticulate")
    Net_integrat <- .integrat_net_tensorly(Net_set,env = env)
  } else{
    # clear register
    if (is.null(n.core_cp) || 1 < n.core_cp) {
      # clear register
      doParallel::stopImplicitCluster()
    }
    # register parallel
    if (n.core_cp != 1) {
      doParallel::registerDoParallel(cores = n.core_cp)
    }    
    
    if (n.core_cp > 1) {
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
  }
  # process bi-direction
  Net_final <- lapply(Net_integrat, function(Net) {
    Net <- .process_biedge(Net,ratio = keep_ratio)
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
#' @param n.core Number of CPU cores to use. Default is all cores - 1.
#' @param perturbed_target Target to be perturbed in dpGRN for ranking cell type. Default is the target in \code{object@para$target}
#' @param resistance_target Target in alternative activated pathway that counteract the effects of anti-cancer drugs. Defaul is \code{NULL}.
#' @param n_dim number of low dimensions in manifold space using manifold alignment. Default is 2.
#' @param n_hop the number of hop for network propagation. Default is 2.
#' @param type characters meaning the MOAs of drug including antagonist or agonist. Default is antagonist.
#' @param method 1 means using the scRank method, 2 means using direct target gene expression, 3 means using downstream target expression. Default is 1.
#' @return scRank object containing the drug-responsive cell type rank.
#' @import dplyr doParallel parallel foreach methods RSpectra RhpcBLASctl Seurat dynamicTreeCut
#' @importFrom magrittr %>%
#' @importFrom scales rescale
#' @importFrom stats dist hclust cutree
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tidyr replace_na
#' @export

rank_celltype <- function(object,
                          n.core = NULL,
                          n_dim = 2,
                          n_hop = 2,
                          perturbed_target = NULL,
                          resistance_target = NULL,
                          type = NULL,
                          method = 1) {
  if (!is(object, "scRank")) {
    stop("Invalid class for object: must be 'scRank'!")
  }
  
  if (is.null(perturbed_target)) {
    perturbed_target <- object@para$target
  }
  if (length(perturbed_target) > 2) {
    stop("The number of perturbed gene should not be more than two!")
  }
  
  if(method == 2){
    # use direct expression to rank cell type
    if(length(perturbed_target) == 1){
      avg_exp <- Seurat::AverageExpression(object@data$seuratObj, group.by = object@para$cell_type, features = perturbed_target) %>% .$RNA %>% t() %>% 
        as.data.frame() %>% rename(Avg_exp = V1) %>% 
        mutate(
          rank = rank(Avg_exp),
          `top_rank%` = rank / dplyr::n(),
          `top_rank%` = scales::rescale(`top_rank%`, c(0,1)),
          fill = `top_rank%`,
          method = 'direct expression'
        )
    }else{
      avg_exp <- Seurat::AverageExpression(object@data$seuratObj, group.by = object@para$cell_type, features = perturbed_target) %>% .$RNA %>% t() %>% 
        as.data.frame() %>% mutate(Avg_exp = rowMeans(select(., 1:2))) %>% select(Avg_exp) %>%
        mutate(
          rank = rank(Avg_exp),
          `top_rank%` = rank / dplyr::n(),
          `top_rank%` = scales::rescale(`top_rank%`, c(0,1)),
          fill = `top_rank%`,
          method = 'direct expression'
        )
    }
    object@cell_type_rank <- avg_exp
    return(object)
  }
  
  if(method == 3){
    Net <- object@net
    ct.keep <- object@para$ct.keep
    tmp <- list()
    tmp.features <- list()
    # use downstream expression to rank cell type
    for(ct in ct.keep){
      if(is.null(sapply(Net,.get_downstream_genes,target = perturbed_target) %>% .[[ct]])){
        tmp[[ct]] <- data.frame(celltype = ct, mean = 0, scale = 0, rank = 0, `top_rank%` = 0)
      } else{
        G.tmp <- .cluster_gene(Net[[ct]],cent_gene = perturbed_target)
        G.tmp <- G.tmp[G.tmp %in% G.tmp[perturbed_target]]
        if(length(G.tmp) > 100){
          G.tmp <- .cluster_gene(Net[[ct]][names(G.tmp),names(G.tmp)],cent_gene = perturbed_target)
          G.tmp <- G.tmp[G.tmp %in% G.tmp[perturbed_target]]
        }
        tmp.features[[ct]] <- names(G.tmp)[!names(G.tmp) %in% perturbed_target]
        tmp[[ct]] <- .cal_module_score(object@data$seuratObj, group_by = object@para$cell_type,features = list(c(tmp.features[[ct]]))) %>% filter(celltype %in% ct) %>% as.data.frame()
      }
    }
    avg_exp <- bind_rows(tmp) %>% mutate(rank = rank(mean)) %>% mutate(`top_rank%` = rescale(rank,c(0,1))) %>%
      mutate(fill = `top_rank%`, method = 'downstream expression') %>% column_to_rownames('celltype') %>% select(mean,rank,`top_rank%`,fill,method)
    object@cell_type_rank <- avg_exp
    return(object)
  }
  
  if(method != 1){
    stop("The method parameter is wrong")
  }
  
  if (!all(object@para$target %in% rownames(object@net[[1]])) || !all(perturbed_target %in% rownames(object@net[[1]]))) {
    stop("Drug target gene is not in the network, please add it into `CreateScRank(target = _target_)'")
  }
  
  # define MOA of drug
  if(is.null(type)){
    MOA = object@para$type
  } else if (type %in% c('antagonist','agonist')){
    MOA = type
  } else {
    stop("Error! Please provide correct type of drug (antagonist or agnoist).")
  }
  
  # clear register
  if (is.null(n.core) || 1 < n.core) {
    # clear register
    doParallel::stopImplicitCluster()
  }
  # register parallel
  if (is.null(n.core)) {
    n.core <- parallel::detectCores() - 2
  }
  n.core <- max(1, n.core)
  if (n.core != 1) {
    doParallel::registerDoParallel(cores = n.core)
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
      if(MOA == 'antagonist'){
        dpGRN[perturbed_target, ] <- 0
      } else if (MOA == 'agonist') {
        if(length(perturbed_target) > 1){
          for(i in perturbed_target){
            dpGRN[i,dpGRN[i,] > 0] <- 1 * sign(dpGRN[i,dpGRN[i,] > 0])
            dpGRN[i,dpGRN[i,] > 0] <- 1 * sign(dpGRN[i,dpGRN[i,] > 0])
          }
        } else {
          dpGRN[perturbed_target,dpGRN[perturbed_target,] > 0] <- 1 * sign(dpGRN[perturbed_target,dpGRN[perturbed_target,] > 0])
        }
      }
      lowDim_net <- .align_net(GRN, dpGRN, ndim = n_dim, n.cores = n.core)
      return(lowDim_net)
    }
  } else {
    lowDim <- lapply(ct.keep, function(ct) {
      GRN <- Net[[ct]]
      dpGRN <- GRN
      if(MOA == 'antagonist'){
        dpGRN[perturbed_target, ] <- 0
      } else if (MOA == 'agonist') {
        if(length(perturbed_target) > 1){
          for(i in perturbed_target){
            dpGRN[i,dpGRN[i,] > 0] <- 1 * sign(dpGRN[i,dpGRN[i,] > 0])
            dpGRN[i,dpGRN[i,] > 0] <- 1 * sign(dpGRN[i,dpGRN[i,] > 0])
          }
        } else {
          dpGRN[perturbed_target,dpGRN[perturbed_target,] > 0] <- 1 * sign(dpGRN[perturbed_target,dpGRN[perturbed_target,] > 0])
        }
      }
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
  if(!is.null(resistance_target)){
    if_in_used_gene <- resistance_target %in% object@para$gene4use
    if(!if_in_used_gene){
      stop("Error: The resistance target is not created in GRN. Please input the gene for `target` in `CreateScRank()`")
    }
    resis_score <- .calculate_score(df, ko = resistance_target, Net = Net, simple_sum = simple_sum, multi_layer = multi_layer)
    resis_score <- as.data.frame(resis_score)
    colnames(resis_score) <- "resis_score"
    integrate_score <- merge(perb_score,resis_score) %>% 
      mutate(integrate_score = perb_score/resis_score) %>%
      mutate(
        rank = rank(integrate_score),
        `top_rank%` = rank / dplyr::n(),
        `top_rank%` = scales::rescale(`top_rank%`, c(0, 1)),
        fill = `top_rank%`,
        method = 'scRank'
      )
    object@meta$dpGRNmeta <- df
    object@cell_type_rank <- integrate_score
    return(object)
  } else{
    perb_score <- perb_score %>% mutate(
      rank = rank(perb_score),
      `top_rank%` = rank / dplyr::n(),
      `top_rank%` = scales::rescale(`top_rank%`, c(0, 1)),
      fill = `top_rank%`,
      method = 'scRank'
    )
    object@meta$dpGRNmeta <- df
    object@cell_type_rank <- perb_score
    return(object)
  }
}


