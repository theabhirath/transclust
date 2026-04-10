# Fraction of Converts with Same Source from Isolate Lookups

Compares two isolate lookups by examining patients categorized as
"convert". This is the internal implementation that works directly with
isolate lookups.

## Usage

``` r
fraction_convert_same_source_from_lookups(
  isolate_lookup1,
  isolate_lookup2,
  surv_df_1,
  surv_df_2
)
```

## Arguments

- isolate_lookup1:

  A data frame from
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md)
  for the first cluster assignment.

- isolate_lookup2:

  A data frame from
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md)
  for the second cluster assignment.

- surv_df:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result.

## Value

A numeric value between 0 and 1 representing the fraction of convert
patients whose source (index or weak-index) is the same in both cluster
assignments. Returns NA if there are no common converts.

## See also

[`fraction_convert_same_source()`](https://theabhirath.github.io/transclust/reference/fraction_convert_same_source.md),
[`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md)
