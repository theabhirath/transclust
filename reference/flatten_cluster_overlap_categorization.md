# Flatten cluster-overlap categorization to a data frame

Converts the named vector output of
[`categorize_cluster_overlap()`](https://theabhirath.github.io/hospitraceR/reference/categorize_cluster_overlap.md)
into a flat data frame with one row per cluster. The patient-level
counterpart is
[`flatten_cluster_patient_categorization()`](https://theabhirath.github.io/hospitraceR/reference/flatten_cluster_patient_categorization.md).

## Usage

``` r
flatten_cluster_overlap_categorization(categorization)
```

## Arguments

- categorization:

  A named vector as returned by
  [`categorize_cluster_overlap()`](https://theabhirath.github.io/hospitraceR/reference/categorize_cluster_overlap.md),
  mapping each cluster ID to its overlap category.

## Value

A data frame with columns `cluster` and `category`.
