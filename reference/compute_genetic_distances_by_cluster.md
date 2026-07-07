# Compute genetic distances by cluster

Summarize intra-cluster SNP distances (mean/median/max) for every
multi-patient cluster in an isolate lookup, via
[`cluster_pairwise_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_pairwise_distances.md).

## Usage

``` r
compute_genetic_distances_by_cluster(isolate_lookup, snp_dist)
```

## Arguments

- isolate_lookup:

  An isolate lookup table from
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

- snp_dist:

  Matrix of SNP distances between isolates.

## Value

Data frame with `mean_genetic_distance`, `median_genetic_distance`,
`max_genetic_distance` and `cluster`, one row per cluster.
