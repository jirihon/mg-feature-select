##
# Feature selection for clustering of metagenomics data.
#
# Author: Jiri Hon <ihon@fit.vutbr.cz>
# Date: 2016/05/15
#

#' One point of joint distribution of the random variables associated with two 
#' clusterings.
#' @see doi:10.1016/j.jmva.2006.11.013, chapter 3.
#'
#' @param i Cluster index from first clustering.
#' @param j Cluster index from second clustering.
#' @param a First clustering.
#' @param b Second clustering.
#' @param n Number of clustered items.
#' @return P(A_i, B_j)
#' 
vi_P <- function(i, j, a, b, n) {
  return(sum(a == i & b == j) / n)
}

#' Mutual information between two clusterings.
#' @see doi:10.1016/j.jmva.2006.11.013, chapter 3.
#' 
#' @param a First clustering.
#' @param b Second clustering.
#' @param n Number of clustered items.
#' @return I(A, B)
#' 
vi_I <- function(a, b, n) {
  i_j <- expand.grid(i = levels(a), j = levels(b))
  pairs <- apply(i_j, 1, function(i_j, a, b, n) {
    i <- i_j[1]
    j <- i_j[2]
    
    P_1_2 = vi_P(i, j, a, b, n);
    P_1 = sum(a == i) / n;
    P_2 = sum(b == j) / n;
    
    if (P_1_2 > 0) {
      return(P_1_2 * log(P_1_2 / (P_1 * P_2)))
    } else {
      return(0)
    }
  }, a, b, n)
  return(sum(pairs))
}

#' Entropy associated with given clustering.
#' @see doi:10.1016/j.jmva.2006.11.013, chapter 3.
#' 
#' @param a Clustering.
#' @param n Number of clustered items.
#' @return Clustering entropy H(A)
#' 
vi_H <- function(a, n) {
  return(sum(sapply(levels(a), function(i, a, n) {
    P <- sum(a == i) / n
    return(-P * log(P))
  }, a, n)))
}

#' Variation of information.
#' @see doi:10.1016/j.jmva.2006.11.013, chapter 3.
#' 
#' @param a First clustering.
#' @param b Second clustering.
#' @param n Number of clustered items.
#' @return VI(A, B)
#' 
vi <- function(a, b, n) {
  return(vi_H(a, n) + vi_H(b, n) - 2 * vi_I(a, b, n))
}

#' Cosine dissimilarity.
#' 
#' @param x Matrix of row vectors.
#' @return Distance matrix.
#' 
cos_dist <- function(x) {
  library(inline)
  
  if (!is.matrix(x)) {
    stop("Matrix expected.")
  }
  cos_dist_c <- rcpp(signature(x = "numeric"), body = '
  Rcpp::NumericMatrix xm(x);
  int n = xm.nrow();
  int m = xm.ncol();
  Rcpp::NumericMatrix dm(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum_ab = 0;
      double sum_a2 = 0;
      double sum_b2 = 0;

      for (int k = 0; k < m; ++k) {
        double a = xm(i, k);
        double b = xm(j, k);
        sum_ab += a * b;
        sum_a2 += a * a;
        sum_b2 += b * b;
      }
      double sum_a2b2 = sum_a2 * sum_b2;
      if (sum_a2b2 > 0)
        dm(i, j) = 1 - (sum_ab / std::sqrt(sum_a2b2));
      else
        dm(i, j) = 1;
    }
  }
  return dm;')
  d <- as.dist(cos_dist_c(x))
  attr(d, 'method') <- 'cosine'
  return(d)
}

#' Find the best tree cutoff.
#' 
#' @param hc Clustering tree.
#' @param org Reference organism mapping.
#' @return Best cutoff (number of clusters) and corresponding VI.
#'
find_best_cutoff <- function(hc, org) {
  best_vi <- .Machine$double.xmax
  best_k <- 0
  
  for (k in 1:(2*length(levels(org)))) {
    cls <- factor(cutree(hc, k = k))
    sim <- vi(org, cls, length(org))
    if (sim < best_vi) {
      best_vi <- sim
      best_k <- k
    }
  }
  return(list(k = best_k, vi = best_vi))
}

#' Evaluate specific distance function on selected features.
#' 
#' @param sel_features Matrix of selected features.
#' @param dist_fn Distance function.
#' @param org Reference organism marking.
#' @return Clustering tree, cutoff and VI.
#' 
dist_eval <- function(sel_features, dist_fn, org) {
  library(fastcluster)
  dist <- do.call(dist_fn, list(sel_features))
  hc <- hclust(dist, method = "average")
  cutoff <- find_best_cutoff(hc, org)
  return(c(tree = list(hc), cutoff))
}

#' Evaluate specific feature combination.
#' 
#' @param names Combination of features.
#' @return Euc. VI and Cos. VI.
#' 
eval_feature_comb <- function(names, features, org) {
  sel_features <- as.matrix(features[ , names])
  euc_res <- dist_eval(sel_features, dist, org)
  cos_res <- dist_eval(sel_features, cos_dist, org)
  res <- list(features = names, euc = euc_res, cos = cos_res)
  print_feature_select_eval(res)
  return(res)
}

#' Print feature evaluation result.
#' 
#' @param res Result of feature evaluation.
#' 
print_feature_select_eval <- function(res) {
  cat(paste(res$features, collapse = " "), "\n")
  cat(sprintf("* EUC: k = %d, vi = %f", res$euc$k, res$euc$vi), "\n")
  cat(sprintf("* COS: k = %d, vi = %f", res$cos$k, res$cos$vi), "\n")
}

#' Load features from file
#' 
#' @param file Input CSV file.
#' @return Features.
#' 
load_features <- function(file) {
  features <- read.csv(file)
  col_names <- names(features)
  
  if (length(col_names) < 2) {
    stop("Input CSV file must contain at least two columns.")
  }
  if (!"Organism" %in% col_names) {
    stop("Missing Organism column in input CSV file.")
  }
  return(features)
}

#' Evaluate all feature combinations for hierarchical clustering.
#' 
#' @param file Input CSV file with features and reference marking.
#' @param min_features Minimal number of features.
#' @param max_features Maximal number of features.
#' @return Evaluations of all feature combinations.
#'
eval_all_feature_comb <- function(
  features, min_features = 3, max_features = .Machine$integer.max,
  parallel = FALSE)
{
  col_names <- names(features)
  # Pick reference organism marking from feature table
  org <- factor(features$Organism)
  org_ann <- data.frame(id = 1:length(org), group = org)
  
  # Consider all other columns to represent various features
  feature_names <- col_names[col_names != "Organism"]
  
  if (length(feature_names) < min_features) {
    stop("Not enough features.")
  }
  if (max_features < min_features) {
    stop("Invalid specification of maximum and minimum features.")
  }
  # Generate all feature combinations
  max_features <- min(max_features, length(feature_names))
  feature_comb <- unlist(
    lapply(min_features:max_features, combn, x = feature_names, simplify = FALSE),
    recursive = FALSE)
  
  if (parallel) {
    library(parallel)
    cl <- makeCluster(detectCores()-1, outfile = "");
    clusterExport(cl, c("cos_dist", "dist_eval", "print_feature_select_eval",
                        "find_best_cutoff", "vi", "vi_H", "vi_I", "vi_P"))
    res <- parLapply(cl, feature_comb, eval_feature_comb, features, org)
  } else {
    res <- lapply(feature_comb, eval_feature_comb, features, org)
  }
  return(res)
}

#' Find best features for specific metric.
#' 
#' @param feat_evals Feature evaluations.
#' @param metric
#' @return Best features, tree, cutoff and VI.
#' 
find_best_features <- function(feat_evals, metric) {
  min_vi <- .Machine$double.xmax
  best <- NULL
  for (res in feat_evals) {
    if (res[[metric]]$vi < min_vi) {
      best <- c(list(features = res$features), res[[metric]])
    }
  }
  return(best)
}

#' Plot clustering tree colored by reference clustering.
#' 
#' @param tree Clustering tree.
#' @param ref Reference clustering.
#' @param file Output file.
#' 
plot_tree <- function(tree, ref, file) {
  library(ape)
  library(ggtree)
  
  ref_ann <- data.frame(id = 1:length(ref), group = ref)
  
  pdf(file, 10, 20)
  print(ggtree(as.phylo(tree), aes(color = group), branch.length = 'none')
        %<+% org_ann + theme(legend.position = "right"))
  dev.off()
}

#' Plot 2D projection of clustering for comparison with reference.
#' 
#' @param x X values
#' @param y Y values
#' @param cls Clustering.
#' @param file Output file.
#' 
plot_clustering_2d <- function(x, y, cls, file) {
  library(ggplot2)
  
  data <- data.frame(x = x, y = y, group = factor(cls))
  pdf(file, 10, 8)
  print(ggplot(data, aes(x = x, y = y, color = group)) + geom_point())
  dev.off()
}
