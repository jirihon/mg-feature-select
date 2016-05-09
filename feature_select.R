
library(ggtree)
library(ape)
library(parallel)

setwd("/home/xhonji01/Projekty/MG/feature_select")

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
  return(as.dist(cos_dist_c(x)))
}

#' Select best features for hierarchical clustering.
#' 
#' @param file Input CSV file with features and reference marking.
#' @param min_features Minimal number of features.
#' @param max_features Maximal number of features.
#'
feature_select <- function(file, min_features = 3, max_features = NULL) {
  features <- read.csv(file)
  col_names <- names(features)
  
  if (length(col_names) < 2) {
    stop("Input CSV file must contain at least two columns.")
  }
  if (!"Organism" %in% col_names) {
    stop("Missing Organism column in input CSV file.")
  }
  # Pick reference organism marking from feature table
  org <- factor(features$Organism)
  org_ann <- data.frame(id = 1:length(org), group = org)
  
  # Consider all other columns to represent various features
  feature_names <- col_names[col_names != "Organism"]
  
  if (length(feature_names) < min_features) {
    stop("Not enough features.")
  }
  if (is.null(max_features)) {
    max_features <- length(feature_names)
  }
  if (max_features < min_features) {
    stop("Invalid specification of maximum and minimum features.")
  }
  
  # Generate all feature combinations
  max_features <- min(max_features, length(feature_names))
  feature_comb <- unlist(lapply(min_features:max_features, combn,
                                x = feature_names, simplify = FALSE),
                         recursive = FALSE)
  
  # n_cores <- detectCores()
  # cl <- makeCluster(n_cores, outfile = "");
  # clusterExport(cl, c("features", "org", "cos_dist", "vi", "vi_H", "vi_I", "vi_P"))
  
  res <- lapply(feature_comb, function(names) {
    library(fastcluster)
    cat(names, "\n")
    
    sel_features <- as.matrix(features[ , names])
    
    # Euclidean
    dist <- dist(sel_features, method = "euclidean")
    hc <- hclust(dist, method = "average")
    cls <- factor(cutree(hc, k = length(levels(org))))
    euc_sim <- vi(org, cls, length(org))
    
    # Cosine
    dist <- cos_dist(sel_features)
    hc <- hclust(dist, method = "average")
    cls <- factor(cutree(hc, k = length(levels(org))))
    cos_sim <- vi(org, cls, length(org))
    
    cat("Euclidean:", euc_sim, "Cosine:", cos_sim, "\n")
    
    return(c(euc_sim, cos_sim))
  })
  
  resm <- matrix(unlist(x), nrow = 2)
  
  # pdf("out/tree.pdf", 8, 15)
  # print(ggtree(as.phylo(hc), aes(color = group), branch.length = "none") %<+% org_ann + theme(legend.position = "right"))
  # dev.off()
}
