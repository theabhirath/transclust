# Calculate isolate-isolate sequential overlap

For every ordered pair of isolates, counts the days on which the
recipient was at a location the donor had occupied on a strictly earlier
day, capturing potential indirect transmission via a shared location.
Pairs with any concurrent co-location are excluded, since that is direct
overlap and is handled by
[`isolate_isolate_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_overlap.md).

## Usage

``` r
isolate_isolate_sequential_overlap(isolate_lookup, trace_mat)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their cluster assignments with other
  relevant epidemiological information. See
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

- trace_mat:

  A patient-by-date matrix of location occupancy (patients in rows,
  dates in columns); a positive value marks a patient's location on a
  given day.

## Value

A data frame with columns `iso_donor`, `iso_recipient` and
`overlap_days`.
