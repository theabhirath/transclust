# Get non-single patient clusters from a vector of cluster assignments.

This function gets non-single patient clusters from a vector of cluster
assignments.

## Usage

``` r
get_non_single_patient_clusters(isolate_lookup)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md).

## Value

A numeric vector of cluster assignments with the non-single patient
clusters.
