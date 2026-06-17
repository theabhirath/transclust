# Compute intra-cluster SNP-distance summaries

The single source of truth for intra-cluster genetic-distance summaries.
Returns the mean, median and max pairwise SNP distance within a single
cluster (the upper triangle of the cluster's own submatrix, excluding
the diagonal). Inter-cluster / inter-isolate context is computed
separately, over all clusters at once, by
[`cluster_inter_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_inter_distances.md).

## Usage

``` r
cluster_pairwise_distances(cluster_seqs, snp_dist)
```

## Arguments

- cluster_seqs:

  Vector of sequence IDs in the cluster.

- snp_dist:

  Matrix of SNP distances between isolates.

## Value

A named numeric vector with `mean_genetic_distance`,
`median_genetic_distance` and `max_genetic_distance`.
