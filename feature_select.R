
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

vi_P <- function(i, j, a, b, n) {
  return(sum(a == i & b == j) / n)
}
vi_I_pair <- function(i, j, a, b, n) {
  P_1_2 = vi_P(i, j, a, b, n);
  
  P_1 = sum(a == i) / n;
  P_2 = sum(b == j) / n;
  
  if (P_1_2 > 0) {
    return(P_1_2 * log(P_1_2 / (P_1 * P_2)))
  } else {
    return(0)
  }
}
vi_I <- function(a, b, n) {
  i_j <- expand.grid(i = levels(a), j = levels(b))
  pairs <- plyr::mdply(i_j, 'vi_I_pair', a, b, n)
  return(sum(pairs$V1))
}
vi_H_single <- function(i, a, n) {
  P <- sum(a == i) / n
  return(-P * log(P))
}
vi_H <- function(a, n) {
  return(sum(sapply(levels(a), 'vi_H_single', a, n)))
}
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
