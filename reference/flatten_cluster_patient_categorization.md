# Flatten cluster-patient categorization to a data frame

Converts the named list output of
[`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md)
into a flat data frame with one row per cluster-patient pair.

## Usage

``` r
flatten_cluster_patient_categorization(categorization)
```

## Arguments

- categorization:

  A named list as returned by
  [`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md).

## Value

A data frame with columns `cluster`, `patient_id`, and `category`.
