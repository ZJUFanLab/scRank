.hvg <- function(input_data,nfeatures = 2000){
    # seurat object
    if(length(VariableFeatures(input_data)) > 0){
        return(Seurat::VariableFeatures(input_data) %>% .[1:2000])
    } else{
        input_data <- Seurat::FindVariableFeatures(input_data,
                                                   selection.method = "vst", 
                                                   nfeatures = nfeatures,
                                                   verbose = FALSE)
        return(Seurat::VariableFeatures(input_data) %>% .[1:2000])
    }
}

.get_ct_data <- function(object,celltype){
    cells <- rownames(object@meta$rawmeta)[object@meta$rawmeta[[object@para$cell_type]] == celltype]
    data <- object@data$data[object@para$gene4use,cells]
}

.make_net <- function(mat,cut_ratio,n.cores){
    # prepare matrix for PCA
    X <- scale(Matrix::t(mat))
    
    # PCR
    if(n.cores > 1){
        Beta <- .PCregression_doParallel(X,nGene = ncol(X),nPC = 3)
        Beta <- do.call(cbind,Beta)
    } else{
        Beta <- .PCregression(X,nGene = ncol(X),nPC = 3)
    }

    Beta <- t(Beta)
    
    A <- 1 - diag(ncol(X))
    for(n in seq_len(ncol(X))){
        A[n, A[n,] == 1] <- Beta[n,]
    }
    
    # scale adjacency matrix
    A_ <- abs(A)
    A <- A/max(A_)
    
    # cut edge
    A[A_ < stats::quantile(A_, cut_ratio)] <- 0
    diag(A) <- 0
    
    colnames(A) <- rownames(A) <- rownames(mat)
    
    A <- as(A,'dgCMatrix')
    
    return(A)
}

.PCregression_doParallel <- function(X, nGene, nPC = 3){
    Beta <- foreach::foreach(n = seq_len(nGene)) %dopar% {
        y <- X[,n]
        Xn <- X[,-n]
        
        V <- RSpectra::svds(Xn, nPC)$v
        X_ <- Xn %*% V
        
        # OLS
        B_ <- Matrix::t(
            Matrix::t(X_)/ (apply(X_,2,function(X) {sqrt(sum(X^2))}))^2
                )
        B_ <- colSums(y * B_)
        B <- V %*% (B_)
        return(B)
    }
}

.PCregression <- function(X, nGene, nPC = 3){
    Beta <- sapply(seq_len(nGene),function(n){
        y <- X[,n]
        Xn <- X[,-n]
        
        V <- RSpectra::svds(Xn, nPC)$v
        X_ <- Xn %*% V
        
        # OLS
        B_ <- t(
            t(X_)/ (apply(X_,2,function(X) {sqrt(sum(X^2))}))^2
        )
        B_ <- colSums(y * B_)
        B <- V %*% B_
        return(B)
    })
}

.integrat_net <- function(net_list,nComp = 5,resolution = 1){
    # initialize tensor chi
    gene_num <- sqrt(unique(lengths(net_list)))
    gene_name <- rownames(net_list[[1]])
    Chi <- array(data = 0, dim = c(gene_num,
                                   gene_num,
                                   1,
                                   length(net_list))
    )
    
    for(i in 1:length(net_list)){
        Chi[,,,i] <- as.matrix(net_list[[i]])
    }
    
    set.seed(1)
    Chi <- rTensor::as.tensor(Chi) %>% 
        rTensor::cp(num_components = 5,max_iter = 1e3,tol = 1e-05)
    Chi_ <- rowSums(Chi$est@data,dims = 2)
    Chi_ <- Chi_/length(net_list)
    Chi_ <- Chi_/max(abs(Chi_))
    Chi_ <- round(Chi_,digits = resolution)
    Chi_ <- methods::as(Chi_, 'dgCMatrix')
    rownames(Chi_) <- colnames(Chi_) <- gene_name
    
    return(Chi_)
}

.integrat_net_tensorly <- function(net_list,nComp = 5,resolution = 1,env){
  gene_name = rownames(net_list[[1]][[1]])
  # write input data
  tensor_list <- lapply(net_list, function(net) {
    gene_num <- sqrt(unique(lengths(net)))
    Chi <- array(data = 0, dim = c(gene_num, gene_num, 1, length(net)))
    for (i in 1:length(net)) {
      Chi[,,,i] <- as.matrix(net[[i]])
    }
    return(Chi)
  })
  
  if (!dir.exists("tensorly_data")) {
    dir.create("tensorly_data")
  }
  path <- getwd()
  for (i in names(tensor_list)) {
    h5write(tensor_list[[i]], file = paste0(path,"/tensorly_data/tensors.h5"), name = paste0("tensor_", i))
  }
  
  res_py <- c("import tensorly as tl", "import numpy as np", "import h5py")
  res_py <- c(res_py,paste0("with h5py.File('",path,"/tensorly_data/tensors.h5', 'r') as file:"))
  res_py <- c(res_py,"    tensor_data = [np.array(file[name]) for name in file.keys()]")
  res_py <- c(res_py,"    names_data = [name for name in file.keys()]")
  res_py <- c(res_py,"results = []")
  res_py <- c(res_py,"for tensor in tensor_data:")
  res_py <- c(res_py,"    tensor_squeezed = np.squeeze(tensor, axis=1)")
  res_py <- c(res_py,"    factors = tl.decomposition.parafac(tensor_squeezed, rank=5, tol=1e-6, n_iter_max=1000, random_state=12)")
  res_py <- c(res_py,"    chi = tl.cp_to_tensor(factors)")
  res_py <- c(res_py,"    chi_ = np.sum(chi, axis = 0)/10")
  res_py <- c(res_py,"    chi_ = np.round(chi_/np.max(abs(chi_)),1)")
  res_py <- c(res_py,"    results.append(chi_)")
  res_py <- c(res_py,paste0("with h5py.File('",path,"/tensorly_data/results.h5', 'w') as file:")) 
  res_py <- c(res_py,"    for name, result in zip(names_data, results):")
  res_py <- c(res_py,"        dataset_name = name.replace('tensor_','result_')")
  res_py <- c(res_py,"        file.create_dataset(dataset_name, data=result)")
  readr::write_lines(x = res_py,file = paste0(path,"/tensorly_data/res_py.py"))
  reticulate::use_condaenv(env)
  reticulate::py_run_file(file = paste0(path,"/tensorly_data/res_py.py"))
  results_list <- list()
  for (i in names(tensor_list)) {
    result_name <- paste0("result_", i)
    results_list[[i]] <- h5read(paste0(path,"/tensorly_data/results.h5"), result_name)
    results_list[[i]] <- methods::as(results_list[[i]], 'dgCMatrix')
    rownames(results_list[[i]]) <- colnames(results_list[[i]]) <- gene_name
  }
  return(results_list)
}

.process_biedge <- function(Net,ratio = 0.25){
    B <- as.matrix(Net)
    B[abs(B) < abs(t(B))] <- 0
    B_ <- (ratio * as.matrix(Net)) + ((1-ratio) * B)
    B_ <- as.matrix(B_)
    diag(B_) <- 0
    B_ <- t(B_)
    return(B_)
}

.align_net <- function(Net_1,Net_2,ndim = 2,n.cores) {
    
    feature_names_1 <- list(rownames(Net_1),colnames(Net_1))
    feature_names_2 <- list(rownames(Net_2),colnames(Net_2))
    
    if(!identical(feature_names_1,feature_names_2) || 
       !identical(feature_names_1[[1]],feature_names_1[[2]]) ||
       !identical(feature_names_2[[1]],feature_names_2[[2]])){
        stop("Network is not shared with the same gene feature!")
    }
    
    Net_12 <- diag(nrow(Net_1))
    laplacian_W <- .manifold_setup(Wx = Net_1,Wy = Net_2,Wxy = Net_12)
    
    RhpcBLASctl::omp_set_num_threads(n.cores)
    RhpcBLASctl::blas_set_num_threads(n.cores)
    
    Eig <- RSpectra::eigs(laplacian_W,ndim*2,'SR')
    Eig$values <- as.numeric(Eig$values)
    Eig$vectors <- apply(Eig$vectors,2,as.numeric)
    # not guaranteed to be in sorted order
    idx <- order(Eig$values)
    Eig$values <- Eig$values[idx]
    Eig$vectors <- Eig$vectors[,idx]
    # discard any with really small eigenvalues
    Eig$vectors <- Eig$vectors[,Eig$values > 1e-8]
    low_embed <- Eig$vectors[,1:ndim]
    colnames(low_embed) <- c("Dim1","Dim2")
    rownames(low_embed) <- c(paste0('GRN_',feature_names_1[[1]]),
                             paste0('dpGRN_',feature_names_1[[2]]))
    return(low_embed)
}

.manifold_setup <- function(Wx,Wy,Wxy,mu = 0.9){
    Wx <- Wx + 1
    Wy <- Wy + 1
    Wxy <- mu * (sum(Wx) + sum(Wy)) / (2 * sum(Wxy)) * Wxy
    W <- rbind(cbind(Wx,Wxy),cbind(t(Wxy),Wy))
    W <- -W # minus sign leads to a copy
    diag(W) <- 0 # set diagonal to zero, incase it isn't already
    diag(W) <- -colSums(W) # re-negate to get positive degrees
    return(W)
}

.create_final_df <- function(df_res,net_list,ko){
    x <- 0
    total <- length(df_res)
    pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
    # Bins: Number of iterations to complete each 1% of the progress bar.
    bins <- ceiling(total / 100)
    df_extend <- list()
    for(j in names(df_res)) {
        x = x + 1
        d <- data.frame()
        
        mtx_grn <- net_list[[j]]
        gene <- df_res[[j]] %>% pull(gene)
        #--- add edge weight value ---
        R <- mtx_grn[ko,gene] 
        In <- colSums(mtx_grn)[gene]
        Out <- rowSums(mtx_grn)[gene]
        #--- add entropy & degree value ---
        mtx_weighted <- apply(mtx_grn, 1, function(X) abs(X)/sum(abs(X)))
        mtx_weighted <- apply(mtx_weighted, 1, tidyr::replace_na, 0)
        degree <- apply(mtx_grn, 1, function(X) sum(abs(X) > 0))
        S.tmp <- apply(mtx_weighted,1,log)
        S.tmp[is.infinite(S.tmp)] <- 0
        S.tmp <- S.tmp * mtx_weighted
        S.tmp <- apply(S.tmp,1,sum)
        S <- -S.tmp
        S <- S[gene]
        degree <- degree[gene]
        
        
        if(length(ko) == 2){
            d <- data.frame(gene = gene,entropy = S,degree = degree, tf_in_weight_1 = R[1,],tf_in_weight_2 = R[2,], out_weight= Out,in_weight=In)
        } else if (length(ko) == 1){
            d <- data.frame(gene = gene,entropy = S,degree = degree, tf_in_weight = R, out_weight= Out,in_weight=In)
        } else{
            message("ko gene must be 1 or 2")
        }
        
        df_extend[[j]] <- d
        ### Update the progress bar.
        if (x[1]%%bins == 0 | x[1] == total) {
            Sys.sleep(0.1)
            utils::setTxtProgressBar(pb, value = x[1])
        }
    }
    
    # 合并
    df_comb <- list()
    for(i in names(df_res)){
        df_comb[[i]] <- cbind(df_res[[i]],df_extend[[i]])
    }
    df_comb <- lapply(df_comb,function(x){
            res <- x[,-c(3,4)]
    })
    return(df_comb)
}

.calculate_score <- function(df_combined,ko,multi_ko = F,zscore = F,multi_layer = F,Net,simple_sum = F){
    if(length(ko) > 1){
        multi_ko = T
    } else{
        multi_ko = F
    }
    score <- c()
    if(simple_sum == FALSE){
        if(zscore == T){
            if(multi_ko == F){
                for (i in names(df_combined)) {
                    
                    WT <- Net[[i]]
                    
                    #--- compute the distance of all target gene  ---
                    tf_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko) %>%
                        mutate(score = exp(Z) * abs(out_weight) / degree) %>%
                        pull(score)
                    if(is.nan(tf_distance.tmp)){
                        tf_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko) %>% 
                            pull(Z) %>% exp(.)
                    }
                    marker_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight != 0) %>%
                        mutate(score = case_when(degree > 0     ~ exp(Z) * abs(tf_in_weight) + exp(Z) * abs(out_weight) / degree,
                                                 degree == 0     ~ exp(Z) * abs(tf_in_weight))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    
                    score.layer2 <- 0
                    if(multi_layer & marker_distance.tmp != 0){
                        node.layer2 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight != 0) %>% pull(gene)
                        score.layer2 <- sapply(node.layer2, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ exp(Z) * abs(node_in_weight) + exp(Z) * abs(out_weight) / degree,
                                                         degree == 0  ~ exp(Z) * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2 <- sum(score.layer2)
                    }
                    # global_tmp <- as.data.frame(df_c57) %>%
                    #   filter(! gene in i) %>%
                    #   pull(distance * out_weight * tf_) %>% sum()
                    
                    score[i] <- tf_distance.tmp + marker_distance.tmp + score.layer2
                }
            }
            if(multi_ko == T){
                for (i in names(df_combined)) {
                    
                    WT <- Net[[i]]
                    
                    #--- compute the distance of all target gene  ---
                    tf_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko[1]) %>%
                        mutate(score = exp(Z) * abs(out_weight) / degree) %>%
                        pull(score) 
                    if(is.nan(tf_distance_1.tmp)){
                        tf_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko[1]) %>% 
                            pull(Z) %>% exp(.)
                    }
                    marker_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight_1 != 0) %>%
                        mutate(score = case_when(degree > 0     ~ exp(Z) * abs(tf_in_weight_1) + exp(Z) * abs(out_weight) / degree,
                                                 degree == 0     ~ exp(Z) * abs(tf_in_weight_1))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    # ko 2
                    tf_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko[2]) %>%
                        mutate(score = exp(Z) * abs(out_weight) / degree) %>%
                        pull(score)
                    if(is.nan(tf_distance_2.tmp)){
                        tf_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko[2]) %>% 
                            pull(Z) %>% exp(.)
                    }
                    marker_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight_2 != 0) %>%
                        mutate(score = case_when(degree > 0     ~ exp(Z) * abs(tf_in_weight_2) + exp(Z) * abs(out_weight) / degree,
                                                 degree == 0     ~ exp(Z) * abs(tf_in_weight_2))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    
                    score.layer2.ko1 <- 0
                    if(multi_layer & marker_distance_1.tmp != 0){
                        node.layer2.ko1 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight_1 != 0) %>% pull(gene)
                        score.layer2.ko1 <- sapply(node.layer2.ko1, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ exp(Z) * abs(node_in_weight) + exp(Z) * abs(out_weight) / degree,
                                                         degree == 0  ~ exp(Z) * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2.ko1 <- sum(score.layer2.ko1)
                    }
                    score.layer2.ko2 <- 0
                    if(multi_layer & marker_distance_2.tmp != 0){
                        node.layer2.ko2 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight_2 != 0) %>% pull(gene)
                        score.layer2.ko2 <- sapply(node.layer2.ko2, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ exp(Z) * abs(node_in_weight) + exp(Z) * abs(out_weight) / degree,
                                                         degree == 0  ~ exp(Z) * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2.ko2 <- sum(score.layer2.ko2)
                    }
                    
                    score[i] <- tf_distance_1.tmp + marker_distance_1.tmp + tf_distance_2.tmp + marker_distance_2.tmp + score.layer2.ko1 + score.layer2.ko2
                }
            }
        }
        
        if(zscore == F){
            if(multi_ko == F){
                for (i in names(df_combined)) {
                    
                    WT <- Net[[i]]
                    
                    #--- compute the distance of all target gene  ---
                    tf_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko) %>%
                        mutate(score = distance * abs(out_weight) / degree) %>%
                        pull(score)
                    if(is.nan(tf_distance.tmp)){
                        tf_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko) %>% 
                            pull(distance) 
                    }
                    marker_distance.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight != 0) %>%
                        mutate(score = case_when(degree > 0     ~ distance * abs(tf_in_weight) + distance * abs(out_weight) / degree,
                                                 degree == 0     ~ distance * abs(tf_in_weight))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    
                    score.layer2 <- 0
                    if(multi_layer & marker_distance.tmp != 0){
                        node.layer2 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight != 0) %>% pull(gene)
                        score.layer2 <- sapply(node.layer2, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ distance * abs(node_in_weight) + distance * abs(out_weight) / degree,
                                                         degree == 0  ~ distance * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2 <- sum(score.layer2)
                    }
                    
                    score[i] <- tf_distance.tmp + marker_distance.tmp + score.layer2
                }
            }
            if(multi_ko == T){
                for (i in names(df_combined)) {
                    
                    WT <- Net[[i]]
                    
                    #--- compute the distance of all target gene  ---
                    tf_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko[1]) %>%
                        mutate(score = distance * abs(out_weight) / degree) %>%
                        pull(score) 
                    if(is.nan(tf_distance_1.tmp)){
                        tf_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko[1]) %>% 
                            pull(distance) 
                    }
                    marker_distance_1.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight_1 != 0) %>%
                        mutate(score = case_when(degree > 0     ~ distance * abs(tf_in_weight_1) + distance * abs(out_weight) / degree,
                                                 degree == 0     ~ distance * abs(tf_in_weight_1))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    # ko 2
                    tf_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(gene %in% ko[2]) %>%
                        mutate(score = distance * abs(out_weight) / degree) %>%
                        pull(score)
                    if(is.nan(tf_distance_2.tmp)){
                        tf_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                            filter(gene %in% ko[2]) %>% 
                            pull(distance)
                    }
                    marker_distance_2.tmp <- as.data.frame(df_combined[[i]]) %>%
                        filter(tf_in_weight_2 != 0) %>%
                        mutate(score = case_when(degree > 0     ~ distance * abs(tf_in_weight_2) + distance * abs(out_weight) / degree,
                                                 degree == 0     ~ distance * abs(tf_in_weight_2))) %>%
                        filter(!is.nan(score)) %>% 
                        pull(score) %>% sum()
                    
                    score.layer2.ko1 <- 0
                    if(multi_layer & marker_distance_1.tmp != 0){
                        node.layer2.ko1 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight_1 != 0) %>% pull(gene)
                        score.layer2.ko1 <- sapply(node.layer2.ko1, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ distance * abs(node_in_weight) + distance * abs(out_weight) / degree,
                                                         degree == 0  ~ distance * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2.ko1 <- sum(score.layer2.ko1)
                    }
                    score.layer2.ko2 <- 0
                    if(multi_layer & marker_distance_2.tmp != 0){
                        node.layer2.ko2 <- as.data.frame(df_combined[[i]]) %>%
                            filter(tf_in_weight_2 != 0) %>% pull(gene)
                        score.layer2.ko2 <- sapply(node.layer2.ko2, function(node){
                            node_in_weight <- WT[node,]
                            layer3.node <- names(node_in_weight)
                            as.data.frame(df_combined[[i]]) %>% mutate(node_in_weight = node_in_weight[df_combined[[i]]$gene]) %>%
                                filter(node_in_weight != 0) %>%
                                mutate(score = case_when(degree > 0  ~ distance * abs(node_in_weight) + distance * abs(out_weight) / degree,
                                                         degree == 0  ~ distance * abs(node_in_weight))) %>%
                                filter(!is.nan(score)) %>%
                                pull(score) %>% sum()
                        })
                        score.layer2.ko2 <- sum(score.layer2.ko2)
                    }
                    
                    
                    score[i] <- tf_distance_1.tmp + marker_distance_1.tmp + tf_distance_2.tmp + marker_distance_2.tmp + score.layer2.ko1 + score.layer2.ko2
                }
            }
        }
    }
    if(simple_sum){
        if(zscore){
            for (i in names(df_combined)) {
                #--- compute the distance of all target gene  ---
                all_distance <- as.data.frame(df_combined[[i]]) %>%
                    pull(Z) %>% exp(.) %>% sum()
                score[i] <- all_distance
            }
        }
        if(zscore == FALSE){
            for (i in names(df_combined)){
                all_distance <- as.data.frame(df_combined[[i]]) %>%
                    pull(distance) %>% sum()
                score[i] <- all_distance
            }
        }
    }
    return(score)
<<<<<<< HEAD
=======
}

# import from enrichplot package
.gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}

.get_downstream_genes <- function(Net,target){
  if(length(target) == 1){
    downstream_genes <- names(Net[target,Net[target,] > 0])
  } else{
    downstream_genes1 <- names(Net[target[1],Net[target[1],] > 0])
    downstream_genes2 <- names(Net[target[2],Net[target[2],] > 0])
    downstream_genes <- union(downstream_genes1,downstream_genes2)
  }
  return(downstream_genes)
}

.cal_module_score <- function(seurat_object,features,group_by){
  seurat_object <- Seurat::AddModuleScore(
    object = seurat_object,
    features = features,
    name = 'downstream_Features',
    ctrl = 10,
    seed = 123,
    slot = 'count'
  )
  normalize_data <- function(data) {
    return ((data - min(data)) / (max(data) - min(data)))
  }
  
  df <- seurat_object@meta.data %>% group_by(.data[[group_by]]) %>% summarise(mean = mean(.data$downstream_Features1)) %>%
    mutate(scale = normalize_data(.data$mean)) %>%
    mutate(
      rank = rank(.data$scale),
      `top_rank%` = rank / dplyr::n(),
      `top_rank%` = scales::rescale(.data$`top_rank%`, c(0, 1))
    ) %>% 
    rename_with(~ "celltype", .cols = names(.)[1])
  return(df)
}

.cluster_gene <- function(Net,cent_gene){
  # Net : adjacency matrix
  # cent_gene : 
  gene_used <- lapply(cent_gene,function(x){
    tmp <- colnames(Net)[abs(Net[x,]) > 0]
    c(tmp,x)
  })
  if(length(gene_used) > 1){
    gene_used <- do.call("union",gene_used)
  } else{
    gene_used <- unlist(gene_used)
  }
  hclust_out <- hclust(stats::dist(1-as.matrix(Net[gene_used,gene_used])))
  gene_groups = dynamicTreeCut::cutreeDynamic(
    hclust_out, 
    minClusterSize = 10, 
    method = "tree",
    deepSplit = FALSE,
    useMedoids = FALSE
  )
  gene_clusters = stats::cutree(hclust_out, length(unique(gene_groups))-1)
  
  # order the gene by gene groups
  order_gene <- gene_clusters %>% sort()
  
  return(order_gene)
>>>>>>> 0357c2efb26c5ed859a881d8f252f6f7b3a06921
}