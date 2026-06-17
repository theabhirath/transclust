# Compute per-cluster nearest-neighbour SNP distances

For every non-singleton cluster, the minimum SNP distance to an isolate
in another cluster (`min_inter_cluster`) and to any isolate outside the
cluster (`min_inter_isolate`). Unlike
[`cluster_pairwise_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_pairwise_distances.md),
this works over the whole clustering at once: it derives each cluster's
"other cluster" and "non cluster" comparison sets from the full
assignment vector. Singleton clusters are dropped via
[`remove_singleton_clusters()`](https://theabhirath.github.io/hospitraceR/reference/remove_singleton_clusters.md),
so their sequences never form a cluster of their own but still count as
unclustered isolates toward `min_inter_isolate`.

## Usage

``` r
cluster_inter_distances(clusters, snp_dist)
```

## Arguments

- clusters:

  A vector named by sequence IDs giving the cluster each sequence
  belongs to.

- snp_dist:

  A matrix of SNP distances between isolates. Its row/column names
  define the full universe of isolates, including those not assigned to
  any cluster.

## Value

A matrix with one row per cluster (named by cluster ID) and columns
`min_inter_cluster` and `min_inter_isolate`.
