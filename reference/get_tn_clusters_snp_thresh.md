# Perform clustering of isolates using a hard SNP distance cutoff.

This function uses a hard SNP distance cutoff to define clusters
naïvely, without using the phylogenetic tree.

## Usage

``` r
get_tn_clusters_snp_thresh(
  snp_dist,
  snp_thresh,
  hclust_method = "complete",
  tree = NULL,
  monophyly_method = c("expand", "break_down")
)
```

## Arguments

- snp_dist:

  A matrix of SNP distances between isolates constructed using a model
  of DNA evolution. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/hospitraceR/reference/get_snp_dist_matrix.md)
  for a useful function to generate this.

- snp_thresh:

  A threshold for defining clusters.

- hclust_method:

  A string indicating the method to use for hierarchical clustering. See
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) for more
  details. Default is "complete".

- tree:

  An optional phylogenetic tree of class `phylo` covering the same
  isolates. When supplied, clusters are forced to be monophyletic with
  respect to this tree (e.g. an NJ or maximum-parsimony tree from
  [`get_phylo_tree()`](https://theabhirath.github.io/hospitraceR/reference/get_phylo_tree.md));
  otherwise monophyly is enforced on the dendrogram implied by the SNP
  distances themselves. Default is `NULL`.

- monophyly_method:

  How to make a non-monophyletic cluster monophyletic. `"expand"` (the
  default) grows the cluster to the smallest clade of the tree that
  contains all of its members, absorbing any intervening isolates.
  `"break_down"` instead splits the cluster into the largest
  monophyletic clades it already contains, without pulling in foreign
  isolates.

## Value

A numeric vector indicating the cluster that each isolate belongs to.
