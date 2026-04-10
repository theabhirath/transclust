# Calculate isolate-isolate sequential overlap

This function calculates the sequential overlap between isolates in a
given trace matrix. Sequential overlap counts the number of days where
the recipient is at a location that the donor had previously been at (on
a strictly earlier day), capturing potential indirect transmission via
shared locations.

## Usage

``` r
isolate_isolate_sequential_overlap(isolate_lookup, trace_mat)
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
