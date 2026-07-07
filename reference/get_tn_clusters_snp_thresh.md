# Cluster isolates using a hard SNP distance cutoff

Defines clusters with a hard SNP distance cutoff, via hierarchical
clustering of the distance matrix cut at `snp_thresh`. This is the naive
baseline to the threshold-free
[`get_tn_clusters_sv_index()`](https://theabhirath.github.io/hospitraceR/reference/get_tn_clusters_sv_index.md).
The resulting clusters are then reconciled with a phylogenetic tree so
that each is monophyletic.

## Usage

``` r
get_tn_clusters_snp_thresh(
  snp_dist,
  snp_thresh,
  hclust_method = "single",
  tree = NULL,
  monophyly_method = c("break_down", "expand")
)
```

## Arguments

- snp_dist:

  A matrix of SNP distances between isolates. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/hospitraceR/reference/get_snp_dist_matrix.md).

- snp_thresh:

  SNP distance at which to cut the hierarchical clustering into
  clusters.

- hclust_method:

  Linkage method for hierarchical clustering, passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

- tree:

  An optional phylogenetic tree of class `phylo` over the same isolates
  (e.g. from
  [`get_phylo_tree()`](https://theabhirath.github.io/hospitraceR/reference/get_phylo_tree.md)).
  When supplied, clusters are forced to be monophyletic with respect to
  it; otherwise monophyly is enforced on the dendrogram implied by the
  SNP distances.

- monophyly_method:

  How a non-monophyletic cluster is reconciled with the tree. `"expand"`
  grows the cluster to the smallest clade containing all its members,
  absorbing any intervening isolates; `"break_down"` splits it into the
  largest monophyletic clades it already contains, pulling in no foreign
  isolates.

## Value

A numeric vector giving the cluster each isolate belongs to.
