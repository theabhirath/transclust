# Calculate isolate-isolate overlap for clusters

This function considers every isolate in a cluster as a recipient, and
returns if there is any overlap when other isolates in the cluster are
the donors.

## Usage

``` r
cluster_isolate_overlap(isolate_lookup, iso_overlap_df)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md).

- iso_overlap_df:

  A data frame with overlap information for isolate pairs. For more
  information, see
  [`isolate_isolate_overlap()`](https://theabhirath.github.io/transclust/reference/isolate_isolate_overlap.md).

## Value

A data frame with columns: cluster, isolate_id, overlap.
