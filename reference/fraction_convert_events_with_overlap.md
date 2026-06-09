# Calculate the fraction of convert events with overlap

This function calculates the fraction of convert events with overlap for
each cluster.

A convert event is an isolate where (1) the patient has no
admission-positive isolates, and (2) the patient had a prior negative
surveillance before this positive isolate. Criterion (2) is taken from
the `prev_surv_neg` flag in the lookup (see
[`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md)):
TRUE only when a surveillance exists before the isolate date and none of
those prior surveillances were positive. Checking the previous
surveillance date alone is not enough, since that prior surveillance may
itself have been positive (meaning the patient was already colonized
rather than converting). Each qualifying isolate counts as one convert
event, and a convert event with an overlap explanation is a convert
event with overlap. The per-cluster fraction is
`convert_events_with_overlap / convert_events`.

The per-cluster counts that make up each fraction are also returned, as
the attributes `n_overlap` (convert events with overlap, the numerator)
and `n_converts` (convert events, the denominator). These let callers
aggregate a pooled, convert-weighted fraction across clusters
(`sum(n_overlap) / sum(n_converts)`) instead of averaging the
per-cluster fractions, which would ignore how many convert events each
cluster contributes.

## Usage

``` r
fraction_convert_events_with_overlap(cluster_overlap_df, isolate_lookup)
```

## Arguments

- cluster_overlap_df:

  A data frame with overlap information for isolate pairs. For more
  information, see
  [`cluster_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/cluster_isolate_overlap.md).

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

## Value

A named numeric vector with the fraction of convert events with overlap
for each cluster (`NA` for clusters with no convert events), carrying
the per-cluster `n_overlap` and `n_converts` counts as attributes.
