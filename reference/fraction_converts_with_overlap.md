# Calculate the fraction of converts with overlap

This function calculates the fraction of converts with overlap for each
cluster.

## Usage

``` r
fraction_converts_with_overlap(cluster_overlap_df, isolate_lookup)
```

## Arguments

- cluster_overlap_df:

  A data frame with overlap information for isolate pairs. For more
  information, see
  [`cluster_isolate_overlap()`](https://theabhirath.github.io/transclust/reference/cluster_isolate_overlap.md).

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md).

## Value

A named vector with the fraction of converts with overlap for each
cluster.
