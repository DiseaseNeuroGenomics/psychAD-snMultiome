
.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ggplot2)
library(dplyr)
library(MOFA2)
library(variancePartition)
library(clue)
read_big_csv <- function(fh,row1asrowname=T){
  library(data.table)
  dat=fread(fh)
  dat=as.data.frame(dat)
  if(row1asrowname){
    rownames(dat)=dat[,1]
    dat[,1]=NULL
  }
  return(dat)
}
aggregate_blocked_factors <- function(
    runs,
    cor_method = "pearson",
    max_iter = 50L,
    tol = 1e-7,
    min_shared_genes = 100L
) {
  
  
  if (is.null(names(runs))) {
    names(runs) <- paste0("run", seq_along(runs))
  }
  
  run_names <- names(runs)
  n_runs <- length(runs)
  
  n_factors <- ncol(runs[[1]]$scores)
  
  if (!all(vapply(
    runs,
    function(x) ncol(x$scores) == n_factors &&
    ncol(x$weights) == n_factors,
    logical(1)
  ))) {
    stop("All runs must have the same number of factors.")
  }
  
  ## Use donors present in every run
  common_donors <- Reduce(
    intersect,lapply(runs, function(x) rownames(x$scores)))
  
  common_donors <- rownames(runs[[1]]$scores)[rownames(runs[[1]]$scores) %in% common_donors]
  
  
  ## ------------------------------------------------------------
  ## Helper functions
  ## ------------------------------------------------------------
  
  zscore_columns <- function(x) {
    x <- as.matrix(x)
    
    center <- colMeans(x, na.rm = TRUE)
    scale_value <- apply(x, 2, stats::sd, na.rm = TRUE)
    
    if (any(!is.finite(scale_value) | scale_value <= 0)) {
      stop("At least one factor has zero or invalid variance.")
    }
    
    z <- sweep(x, 2, center, FUN = "-")
    z <- sweep(z, 2, scale_value, FUN = "/")
    
    list(
      x = z,
      center = center,
      scale = scale_value
    )
  }
  
  fisher_mean_cor <- function(x) {
    x <- x[is.finite(x)]
    
    if (length(x) == 0L) {
      return(NA_real_)
    }
    
    eps <- 1e-12
    x <- pmax(pmin(x, 1 - eps), -1 + eps)
    
    tanh(mean(atanh(x)))
  }
  
  match_factors <- function(reference_scores, target_scores) {
    correlation_matrix <- suppressWarnings(
      stats::cor(
        reference_scores,
        target_scores,
        use = "pairwise.complete.obs",
        method = cor_method
      )
    )
    
    correlation_matrix[!is.finite(correlation_matrix)] <- 0
    
    permutation <- as.integer(
      clue::solve_LSAP(
        abs(correlation_matrix),
        maximum = TRUE
      )
    )
    
    matched_cor <- correlation_matrix[
      cbind(seq_len(nrow(correlation_matrix)), permutation)
    ]
    
    sign_multiplier <- ifelse(matched_cor < 0, -1, 1)
    
    list(
      permutation = permutation,
      sign = sign_multiplier,
      raw_cor = matched_cor,
      mean_abs_cor = mean(abs(matched_cor))
    )
  }
  
  pairwise_cor_matrix <- function(x, minimum_n = 3L) {
    x <- as.matrix(x)
    nr <- ncol(x)
    
    result <- matrix(
      NA_real_,
      nrow = nr,
      ncol = nr,
      dimnames = list(colnames(x), colnames(x))
    )
    
    diag(result) <- 1
    
    if (nr < 2L) {
      return(result)
    }
    
    for (i in seq_len(nr - 1L)) {
      for (j in (i + 1L):nr) {
        keep <- is.finite(x[, i]) & is.finite(x[, j])
        
        if (sum(keep) >= minimum_n) {
          value <- suppressWarnings(
            stats::cor(
              x[keep, i],
              x[keep, j],
              method = cor_method
            )
          )
        } else {
          value <- NA_real_
        }
        
        result[i, j] <- value
        result[j, i] <- value
      }
    }
    
    result
  }
  
  ## ------------------------------------------------------------
  ## 2. Normalize score scale
  ##
  ## If:
  ## X ~= scores %*% t(weights)
  ##
  ## scores are converted to SD = 1 and weights are multiplied by
  ## the original score SD. This preserves the factor contribution.
  ## ------------------------------------------------------------
  
  normalized_runs <- lapply(runs, function(x) {
    original_scores <- x$scores[common_donors, , drop = FALSE]
    
    if (anyNA(original_scores)) {
      stop("Missing donor factor scores are not currently supported.")
    }
    
    normalized <- zscore_columns(original_scores)
    
    normalized_weights <- sweep(
      x$weights,
      2,
      normalized$scale,
      FUN = "*"
    )
    
    list(
      scores = normalized$x,
      weights = normalized_weights,
      score_center = normalized$center,
      score_scale = normalized$scale,
      original_factor_names = colnames(original_scores)
    )
  })
  
  ## ------------------------------------------------------------
  ## 3. Choose the medoid run
  ##
  ## This avoids arbitrarily using run 1 as the reference.
  ## ------------------------------------------------------------
  
  run_similarity <- matrix(
    1,
    nrow = n_runs,
    ncol = n_runs,
    dimnames = list(run_names, run_names)
  )
  
  if (n_runs > 1L) {
    for (i in seq_len(n_runs - 1L)) {
      for (j in (i + 1L):n_runs) {
        matched <- match_factors(
          normalized_runs[[i]]$scores,
          normalized_runs[[j]]$scores
        )
        
        run_similarity[i, j] <- matched$mean_abs_cor
        run_similarity[j, i] <- matched$mean_abs_cor
      }
    }
  }
  
  medoid_similarity <- if (n_runs > 1L) {
    (rowSums(run_similarity) - 1) / (n_runs - 1)
  } else {
    1
  }
  
  reference_index <- which.max(medoid_similarity)
  reference_run <- run_names[reference_index]
  
  factor_names <- normalized_runs[[reference_index]]$original_factor_names
  
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(n_factors))
  }
  
  template <- normalized_runs[[reference_index]]$scores
  colnames(template) <- factor_names
  
  align_one_run <- function(x, current_template) {
    matched <- match_factors(
      current_template,
      x$scores
    )
    
    aligned_scores <- x$scores[
      , matched$permutation,
      drop = FALSE
    ]
    
    aligned_weights <- x$weights[
      , matched$permutation,
      drop = FALSE
    ]
    
    aligned_scores <- sweep(
      aligned_scores,
      2,
      matched$sign,
      FUN = "*"
    )
    
    aligned_weights <- sweep(
      aligned_weights,
      2,
      matched$sign,
      FUN = "*"
    )
    
    colnames(aligned_scores) <- factor_names
    colnames(aligned_weights) <- factor_names
    
    list(
      scores = aligned_scores,
      weights = aligned_weights,
      permutation = matched$permutation,
      sign = matched$sign,
      raw_cor = matched$raw_cor,
      original_factor = x$original_factor_names[
        matched$permutation
      ]
    )
  }
  
  ## ------------------------------------------------------------
  ## 4. Iterative matching to a consensus template
  ## ------------------------------------------------------------
  
  converged <- FALSE
  iteration <- 0L
  
  for (iteration in seq_len(max_iter)) {
    aligned_runs <- lapply(
      normalized_runs,
      align_one_run,
      current_template = template
    )
    
    new_template <- Reduce(
      "+",
      lapply(aligned_runs, function(x) x$scores)
    ) / n_runs
    
    new_template <- zscore_columns(new_template)$x
    colnames(new_template) <- factor_names
    
    difference <- max(
      abs(new_template - template),
      na.rm = TRUE
    )
    
    template <- new_template
    
    if (difference < tol) {
      converged <- TRUE
      break
    }
  }
  
  ## Re-align once to the final template
  aligned_runs <- lapply(
    normalized_runs,
    align_one_run,
    current_template = template
  )
  
  names(aligned_runs) <- run_names
  
  ## ------------------------------------------------------------
  ## 5. Alignment table
  ## ------------------------------------------------------------
  
  factor_alignment <- do.call(
    rbind,
    lapply(seq_along(aligned_runs), function(r) {
      data.frame(
        run = run_names[r],
        consensus_factor = factor_names,
        original_factor =
          aligned_runs[[r]]$original_factor,
        sign_multiplier =
          aligned_runs[[r]]$sign,
        correlation_before_sign =
          aligned_runs[[r]]$raw_cor,
        correlation_after_sign =
          abs(aligned_runs[[r]]$raw_cor),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(factor_alignment) <- NULL
  
  ## ------------------------------------------------------------
  ## 6. Donor score consensus
  ## ------------------------------------------------------------
  
  score_array <- array(
    NA_real_,
    dim = c(
      length(common_donors),
      n_factors,
      n_runs
    ),
    dimnames = list(
      donor = common_donors,
      factor = factor_names,
      run = run_names
    )
  )
  
  for (r in seq_along(aligned_runs)) {
    score_array[, , r] <- aligned_runs[[r]]$scores
  }
  
  ## Mean of the five aligned standardized scores
  consensus_scores <- apply(
    score_array,
    c(1, 2),
    mean,
    na.rm = TRUE
  )
  
  consensus_score_mean <- colMeans(consensus_scores)
  consensus_score_sd <- apply(
    consensus_scores,
    2,
    stats::sd
  )
  
  ## Optional SD=1 version for downstream regression
  consensus_scores_z <- zscore_columns(consensus_scores)$x
  colnames(consensus_scores_z) <- factor_names
  rownames(consensus_scores_z) <- common_donors
  
  ## ------------------------------------------------------------
  ## 7. Donor-score reproducibility
  ##
  ## Per-round score:
  ## correlation of one run with the consensus of the other runs.
  ##
  ## Final score:
  ## Fisher-z mean of all pairwise run correlations.
  ## ------------------------------------------------------------
  
  donor_loo_cor <- matrix(
    NA_real_,
    nrow = n_runs,
    ncol = n_factors,
    dimnames = list(run_names, factor_names)
  )
  
  donor_full_consensus_cor <- donor_loo_cor
  
  donor_pairwise_cor <- vector("list", n_factors)
  names(donor_pairwise_cor) <- factor_names
  
  donor_final_reproducibility <- numeric(n_factors)
  donor_pairwise_median <- numeric(n_factors)
  donor_pairwise_min <- numeric(n_factors)
  
  for (k in seq_len(n_factors)) {
    score_matrix <- score_array[, k, , drop = FALSE]
    score_matrix <- matrix(
      score_matrix,
      nrow = length(common_donors),
      ncol = n_runs,
      dimnames = list(common_donors, run_names)
    )
    
    for (r in seq_len(n_runs)) {
      other_runs <- setdiff(seq_len(n_runs), r)
      
      loo_consensus <- rowMeans(
        score_matrix[, other_runs, drop = FALSE]
      )
      
      donor_loo_cor[r, k] <- suppressWarnings(
        stats::cor(
          score_matrix[, r],
          loo_consensus,
          method = cor_method
        )
      )
      
      donor_full_consensus_cor[r, k] <- suppressWarnings(
        stats::cor(
          score_matrix[, r],
          consensus_scores[, k],
          method = cor_method
        )
      )
    }
    
    pair_cor <- pairwise_cor_matrix(
      score_matrix,
      minimum_n = 3L
    )
    
    donor_pairwise_cor[[k]] <- pair_cor
    
    pair_values <- pair_cor[lower.tri(pair_cor)]
    pair_values <- pair_values[is.finite(pair_values)]
    
    donor_final_reproducibility[k] <-
      fisher_mean_cor(pair_values)
    
    donor_pairwise_median[k] <- if (length(pair_values)) {
      median(pair_values)
    } else {
      NA_real_
    }
    
    donor_pairwise_min[k] <- if (length(pair_values)) {
      min(pair_values)
    } else {
      NA_real_
    }
  }
  
  ## ------------------------------------------------------------
  ## 8. Gene weight consensus
  ##
  ## Genes excluded in a particular chromosome block remain NA.
  ## The consensus is calculated over runs in which the gene was
  ## included.
  ## ------------------------------------------------------------
  
  all_genes <- Reduce(
    union,
    lapply(aligned_runs, function(x) rownames(x$weights))
  )
  
  weight_array <- array(
    NA_real_,
    dim = c(
      length(all_genes),
      n_factors,
      n_runs
    ),
    dimnames = list(
      gene = all_genes,
      factor = factor_names,
      run = run_names
    )
  )
  
  for (r in seq_along(aligned_runs)) {
    genes_in_run <- rownames(aligned_runs[[r]]$weights)
    
    weight_array[
      genes_in_run,
      ,
      r
    ] <- aligned_runs[[r]]$weights
  }
  
  consensus_weights <- apply(
    weight_array,
    c(1, 2),
    function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    }
  )
  
  weight_n_runs <- apply(
    is.finite(weight_array),
    c(1, 2),
    sum
  )
  
  weight_sd <- apply(
    weight_array,
    c(1, 2),
    function(x) {
      x <- x[is.finite(x)]
      
      if (length(x) >= 2L) {
        stats::sd(x)
      } else {
        NA_real_
      }
    }
  )
  
  weight_se <- weight_sd / sqrt(weight_n_runs)
  
  weight_sign_agreement <- apply(
    weight_array,
    c(1, 2),
    function(x) {
      x <- x[is.finite(x)]
      
      if (length(x) == 0L) {
        return(NA_real_)
      }
      
      positive <- sum(x > 0)
      negative <- sum(x < 0)
      nonzero <- positive + negative
      
      if (nonzero == 0L) {
        return(1)
      }
      
      max(positive, negative) / nonzero
    }
  )
  
  ## These weights correspond to consensus_scores_z.
  ##
  ## consensus_scores =
  ##     consensus_scores_z * consensus_score_sd
  ##
  ## Therefore:
  ## consensus_scores %*% weights
  ## =
  ## consensus_scores_z %*%
  ##     (weights * consensus_score_sd)
  consensus_weights_for_z_scores <- sweep(
    consensus_weights,
    2,
    consensus_score_sd,
    FUN = "*"
  )
  
  ## ------------------------------------------------------------
  ## 9. Gene-weight reproducibility by factor
  ## ------------------------------------------------------------
  
  weight_loo_cor <- matrix(
    NA_real_,
    nrow = n_runs,
    ncol = n_factors,
    dimnames = list(run_names, factor_names)
  )
  
  weight_pairwise_cor <- vector("list", n_factors)
  names(weight_pairwise_cor) <- factor_names
  
  weight_final_reproducibility <- numeric(n_factors)
  weight_pairwise_median <- numeric(n_factors)
  weight_pairwise_min <- numeric(n_factors)
  
  for (k in seq_len(n_factors)) {
    weight_matrix <- weight_array[, k, , drop = FALSE]
    weight_matrix <- matrix(
      weight_matrix,
      nrow = length(all_genes),
      ncol = n_runs,
      dimnames = list(all_genes, run_names)
    )
    
    for (r in seq_len(n_runs)) {
      other_runs <- setdiff(seq_len(n_runs), r)
      
      other_matrix <- weight_matrix[
        ,
        other_runs,
        drop = FALSE
      ]
      
      n_other <- rowSums(is.finite(other_matrix))
      
      loo_weight <- rowMeans(
        other_matrix,
        na.rm = TRUE
      )
      
      loo_weight[n_other < 2L] <- NA_real_
      
      keep <- is.finite(weight_matrix[, r]) &
        is.finite(loo_weight)
      
      if (sum(keep) >= min_shared_genes) {
        weight_loo_cor[r, k] <- suppressWarnings(
          stats::cor(
            weight_matrix[keep, r],
            loo_weight[keep],
            method = cor_method
          )
        )
      }
    }
    
    pair_cor <- pairwise_cor_matrix(
      weight_matrix,
      minimum_n = min_shared_genes
    )
    
    weight_pairwise_cor[[k]] <- pair_cor
    
    pair_values <- pair_cor[lower.tri(pair_cor)]
    pair_values <- pair_values[is.finite(pair_values)]
    
    weight_final_reproducibility[k] <-
      fisher_mean_cor(pair_values)
    
    weight_pairwise_median[k] <- if (length(pair_values)) {
      median(pair_values)
    } else {
      NA_real_
    }
    
    weight_pairwise_min[k] <- if (length(pair_values)) {
      min(pair_values)
    } else {
      NA_real_
    }
  }
  
  ## ------------------------------------------------------------
  ## 10. Summary reproducibility table
  ## ------------------------------------------------------------
  
  factor_reproducibility <- data.frame(
    factor = factor_names,
    
    ## Primary factor reproducibility score
    reproducibility_score =
      donor_final_reproducibility,
    
    donor_pairwise_median =
      donor_pairwise_median,
    
    donor_pairwise_min =
      donor_pairwise_min,
    
    weight_reproducibility_score =
      weight_final_reproducibility,
    
    weight_pairwise_median =
      weight_pairwise_median,
    
    weight_pairwise_min =
      weight_pairwise_min,
    
    consensus_score_sd =
      consensus_score_sd,
    
    stringsAsFactors = FALSE
  )
  
  safe_run_names <- make.names(run_names, unique = TRUE)
  
  for (r in seq_len(n_runs)) {
    factor_reproducibility[[paste0("donor_LOO_", safe_run_names[r])]] <- donor_loo_cor[r, ]
    
    factor_reproducibility[[paste0("donor_to_full_consensus_", safe_run_names[r])]] <- donor_full_consensus_cor[r, ]
    
    factor_reproducibility[[paste0("weight_LOO_", safe_run_names[r])]] <- weight_loo_cor[r, ]
  }
  
  ## ------------------------------------------------------------
  ## 11. Wide donor score table
  ## ------------------------------------------------------------
  
  donor_scores_long <- do.call(
    rbind,
    lapply(seq_len(n_factors), function(k) {
      round_values <- score_array[, k, , drop = FALSE]
      
      round_values <- matrix(
        round_values,
        nrow = length(common_donors),
        ncol = n_runs,
        dimnames = list(
          common_donors,
          paste0("score_", safe_run_names)
        )
      )
      
      result <- data.frame(
        donor = common_donors,
        factor = factor_names[k],
        round_values,
        consensus = consensus_scores[, k],
        consensus_z = consensus_scores_z[, k],
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      result
    })
  )
  
  rownames(donor_scores_long) <- NULL
  
  ## ------------------------------------------------------------
  ## 12. Wide gene-weight table
  ## ------------------------------------------------------------
  
  gene_weights_long <- do.call(
    rbind,
    lapply(seq_len(n_factors), function(k) {
      round_values <- weight_array[, k, , drop = FALSE]
      
      round_values <- matrix(
        round_values,
        nrow = length(all_genes),
        ncol = n_runs,
        dimnames = list(
          all_genes,
          paste0("weight_", safe_run_names)
        )
      )
      
      result <- data.frame(
        gene = all_genes,
        factor = factor_names[k],
        round_values,
        
        ## Corresponds to consensus_scores
        consensus = consensus_weights[, k],
        
        ## Corresponds to consensus_scores_z
        consensus_for_z_score =
          consensus_weights_for_z_scores[, k],
        
        n_runs = weight_n_runs[, k],
        weight_sd = weight_sd[, k],
        weight_se = weight_se[, k],
        sign_agreement =
          weight_sign_agreement[, k],
        
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      result
    })
  )
  
  rownames(gene_weights_long) <- NULL
  
  ## ------------------------------------------------------------
  ## Return everything
  ## ------------------------------------------------------------
  
  list(
    reference_run = reference_run,
    converged = converged,
    iterations = iteration,
    
    run_similarity = run_similarity,
    factor_alignment = factor_alignment,
    factor_reproducibility =
      factor_reproducibility,
    
    ## Convenient final tables
    donor_scores_long = donor_scores_long,
    gene_weights_long = gene_weights_long,
    
    ## Consensus matrices
    consensus_scores = consensus_scores,
    consensus_scores_z = consensus_scores_z,
    consensus_weights = consensus_weights,
    consensus_weights_for_z_scores =
      consensus_weights_for_z_scores,
    
    ## Each aligned round
    aligned_scores = setNames(
      lapply(aligned_runs, function(x) x$scores),
      run_names
    ),
    
    aligned_weights = setNames(
      lapply(aligned_runs, function(x) x$weights),
      run_names
    ),
    
    ## Full arrays
    score_array = score_array,
    weight_array = weight_array,
    
    ## Detailed correlation matrices
    donor_pairwise_cor = donor_pairwise_cor,
    weight_pairwise_cor = weight_pairwise_cor,
    
    ## Original score normalization parameters
    score_normalization = setNames(
      lapply(normalized_runs, function(x) {
        list(
          center = x$score_center,
          scale = x$score_scale
        )
      }),
      run_names
    )
  )
}
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')
ind=expand.grid(1:5,10,'top25')
runs=mapply(function(x,y,z){
  MOFAobject=load_model(sprintf('MOFA/MOFA_%s_%s_%s.hdf5',x,y,z))   
  list(scores=get_factors(MOFAobject)$single_group,
       weights=do.call(rbind,get_weights(MOFAobject)))
},as.character(ind[,3]),ind[,1],ind[,2],SIMPLIFY = F)
names(runs)=ind[,1]


agg=aggregate_blocked_factors(runs)
maps=lapply(split(data.frame(agg$factor_alignment),agg$factor_alignment$run),function(x){
  y=x$consensus_factor;names(y)=x$original_factor;y
})
varExps=mapply(function(x,y,z){
  MOFAobject=load_model(sprintf('MOFA/MOFA_%s_%s_%s.hdf5',x,y,z))   
  a=get_variance_explained(MOFAobject)
  b=a[[2]][[1]]
  rownames(b)=maps[[as.character(y)]][rownames(b)]
  list(total=a[[1]][[1]],perfactor=b)
},as.character(ind[,3]),ind[,1],ind[,2],SIMPLIFY = F)
names(varExps)=ind[,1]

tot_var=sapply(varExps,function(x)x[[1]])
colnames(tot_var)=paste0('X',colnames(tot_var))
factor_var=do.call(rbind,
                  lapply(maps[[1]],function(x){
                    z=sapply(varExps,function(y)y[[2]][x,])
                    colnames(z)=paste0('X',colnames(z))
                    z=data.frame(factor=x,assay=rownames(z),z,Mean=rowMeans(z))
                  }))

save(agg,maps,factor_var,tot_var,file='MOFA_top_25_10_agg.Rdata')
write.csv(agg$consensus_scores,file='MOFA_top_25_10_agg_factors.csv')
write.csv(agg$consensus_weights,file='MOFA_top_25_10_agg_weights.csv')


