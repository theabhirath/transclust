# Compute intra-cluster SNP-distance summaries

Returns the mean, median and maximum pairwise SNP distance within a
single cluster, taken from the upper triangle of the cluster's own
submatrix (excluding the diagonal). Inter-cluster and inter-isolate
distances are computed separately, over all clusters at once, by
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
