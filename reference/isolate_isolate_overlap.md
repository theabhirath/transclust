# Calculate isolate-isolate overlap

This function calculates the overlap between isolates in a given trace
matrix.

## Usage

``` r
isolate_isolate_overlap(isolate_lookup, trace_mat)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their clusters assignments which has
  other relevant epidemiological information. For more information, see
  [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md).

- trace_mat:

  A matrix with rows representing days and columns representing
  patients.

## Value

A data frame with columns: iso_donor, iso_recipient, overlap_days.
