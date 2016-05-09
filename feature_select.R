
library(ggtree)
library(ape)
library(dplyr)

setwd("/home/xhonji01/Projekty/MG/feature_select")


features <- read.csv("features.csv")
col_names <- names(features)

if (length(col_names) < 2) {
  stop("Input CSV file must contain at least two columns.")
}
if (!"Organism" %in% col_names) {
  stop("Missing Organism column in input CSV file.")
}

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


# Pick reference organism marking from feature table
org <- factor(features$Organism)
org_ann <- data.frame(id = seq(1, length(org)), group = org)

# Consider all other columns to represent various features
feature_names <- col_names[col_names != "Organism"]

sel_features <- features[ , feature_names]
dist <- dist(sel_features)

hc <- hclust(dist, "average")
cls <- factor(cutree(hc, k = 10))

# Cartesian product of reference cluster ids and computed cluster ids
i_j <- expand.grid(levels(org), levels(cls))

# pdf("out/tree.pdf", 8, 15)
print(ggtree(as.phylo(hc), aes(color = group)) %<+% org_ann + theme(legend.position = "right"))
# dev.off()
