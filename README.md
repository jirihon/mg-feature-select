# Metagenomic Feature Select

R/Bioconductor script for hierachical clustering by various DNA sequence features and measuring its similarity with reference clustering.

## How to Run

1. Install the latest stable version of R (3.3.0) from [CRAN](http://mirrors.nic.cz/R/).

2. Run the latest R and install packages `Rcpp`, `inline`, `fastcluster` and `ggtree` using Bioconductor installer.

    ```R
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("Rcpp", "inline", "fastcluster", "ggtree"))
    ```
3. Download the script `feature_select.R` from GitHub into current working folder and load it into R.

    ```R
    source("feature_select.R")
    ```

4. Use functions `load_features`, `eval_all_feature_comb`, `find_best_features` and `plot_tree` to run the analysis and visualize the results. The most time consuming function `eval_all_feature_comb` can be run in parallel mode (set `parallel` option to `TRUE`). See the following example:

    ```R
    features <- load_features("features.csv")
    res <- eval_all_feature_comb(features, min_features = 2, max_features = 6, parallel = TRUE)
    best_euc <- find_best_features(res, 'euc') # for euclidean distance
    best_cos <- find_best_features(res, 'cos') # for cosine dissimilarity
    plot_tree(best_euc$tree, features$Organism, "tree_euc.pdf")
    plot_tree(best_cos$tree, features$Organism, "tree_cos.pdf")
    ```

5. Feel free to inspect the structure of `best_euc` and `best_cos` variables. It is just a list of three items: `features` contains a vector of used features, `k` is a cutoff value (number of clusters) of the hierarchical clustering tree, that gave the best Variantion of Information (VI) with reference clustering and `vi` contains the best VI value.

6. You can also plot a 2D projection of the clustering for comparison with reference:

    ```R
    ref_cls <- features$Organism
    plot_clustering_2d(features$HP1_aktivita, features$HP2_mobilita, ref_cls, "2d_ref.pdf")
    best_cls <- cutree(best_euc$tree, best_euc$k)
    plot_clustering_2d(features$HP1_aktivita, features$HP2_mobilita, best_cls, "2d_euc.pdf")
    ```
