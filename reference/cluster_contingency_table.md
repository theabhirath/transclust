# Create Contingency Table from Two Cluster Assignments

Generates a contingency table comparing two cluster assignments. The
table counts the number of isolates shared between each pair of clusters
from the two assignments. Only isolates present in both assignments are
considered.

## Usage

``` r
cluster_contingency_table(clusters1, clusters2)
```

## Arguments

- clusters1:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- clusters2:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

## Value

A matrix where rows correspond to clusters in `clusters1`, columns
correspond to clusters in `clusters2`, and cell values represent the
count of shared isolates.
