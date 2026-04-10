# Adjusted Mutual Information from Cluster Assignments

Computes the Adjusted Mutual Information (AMI) between two cluster
assignments. This is a convenience wrapper around
[`ami_from_contingency()`](https://theabhirath.github.io/transclust/reference/ami_from_contingency.md).

## Usage

``` r
adjusted_mutual_information(clusters1, clusters2)
```

## Arguments

- clusters1:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- clusters2:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

## Value

A numeric value representing the Adjusted Mutual Information.

## See also

[`ami_from_contingency()`](https://theabhirath.github.io/transclust/reference/ami_from_contingency.md),
[`cluster_contingency_table()`](https://theabhirath.github.io/transclust/reference/cluster_contingency_table.md)
