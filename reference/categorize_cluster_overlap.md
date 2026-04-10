# Categorize cluster overlap

This function categorizes the overlap of isolates in a cluster, based on
the earliest collected isolate from all patients deemed to be in the
cluster (hereafter referred to as the "index isolate"), and overlap
explanations for other isolates in the cluster. The categories are:

- "patient-to-patient": if the index isolate is admission-positive and
  there is overlap explanation for all other converts in the cluster.

- "weak-index": if the index isolate is not admission-positive but is
  the first surveillance for the patient after admission and there is
  overlap explanation for all other isolates in the cluster.

- "missing-intermediate": if the index isolate is admission-positive but
  at least one other convert in the cluster has no overlap explanation.

- "false-negative-index": if the index isolate is not admission-positive
  but there is overlap explanation for all other converts in the cluster
  barring one (deemed to be the false negative index).

- "missing-source": if the index isolate is not admission-positive and
  there is no overlap explanation for more than one convert in the
  cluster.

- "multiply-colonized-index": if the index isolate is not in the cluster
  but is admission-positive, this is a "multiply-colonized index" if
  there is overlap explanation for all other isolates in the cluster.

- "multiply-colonized-index-missing-intermediate": if the index isolate
  is not in the cluster but is admission-positive, and at least one
  other convert in the cluster has no overlap explanation.

- "inexplicable": catch-all category for cases that do not fit into the
  other categories. Currently, index and overlap explanations are not
  enough to explain these clusters.

## Usage

``` r
categorize_cluster_overlap(isolate_lookup, cluster_overlap_df, surv_df)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md).

- cluster_overlap_df:

  A data frame with overlap information for isolate pairs. For more
  information, see
  [`cluster_isolate_overlap()`](https://theabhirath.github.io/transclust/reference/cluster_isolate_overlap.md).

- surv_df:

  A data frame with surveillance information for isolates.

## Value

A vector with the categorization of the overlap for each cluster.
