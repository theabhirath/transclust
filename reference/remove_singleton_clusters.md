# Remove singleton clusters from a vector of cluster assignments

Drops clusters made up of a single sequence, preserving the original
cluster values and isolate names of the rest. Note that a singleton is
counted by sequence here, unlike
[`get_non_single_patient_clusters()`](https://theabhirath.github.io/hospitraceR/reference/get_non_single_patient_clusters.md),
which counts distinct patients.

## Usage

``` r
remove_singleton_clusters(clusters)
```

## Arguments

- clusters:

  A numeric vector of cluster assignments.

## Value

The input vector with the singleton-cluster entries removed.
