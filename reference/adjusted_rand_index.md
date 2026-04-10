# Adjusted Rand Index from Cluster Assignments

Computes the Adjusted Rand Index (ARI) between two cluster assignments.
This is a convenience wrapper around
[`ari_from_contingency()`](https://theabhirath.github.io/transclust/reference/ari_from_contingency.md).

## Usage

``` r
adjusted_rand_index(clusters1, clusters2)
```

## Arguments

- clusters1:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- clusters2:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

## Value

A numeric value representing the Adjusted Rand Index.

## See also

[`ari_from_contingency()`](https://theabhirath.github.io/transclust/reference/ari_from_contingency.md),
[`cluster_contingency_table()`](https://theabhirath.github.io/transclust/reference/cluster_contingency_table.md)
