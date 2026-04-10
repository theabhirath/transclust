# Perform clustering of isolates using a hard SNP distance cutoff.

This function uses a hard SNP distance cutoff to define clusters
naïvely, without using the phylogenetic tree.

## Usage

``` r
get_tn_clusters_snp_thresh(snp_dist, snp_thresh, hclust_method = "complete")
```

## Arguments

- snp_dist:

  A matrix of SNP distances between isolates constructed using a model
  of DNA evolution. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/transclust/reference/get_snp_dist_matrix.md)
  for a useful function to generate this.

- snp_thresh:

  A threshold for defining clusters.

- hclust_method:

  A string indicating the method to use for hierarchical clustering. See
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) for more
  details. Default is "complete".

## Value

A numeric vector indicating the cluster that each isolate belongs to.
