# Calculate isolate-isolate co-location overlap

For every ordered pair of isolates, counts the days on which the two
patients were in the same place in `trace_mat`, within the window from
the donor's previous surveillance date up to the recipient's collection
date. This is the concurrent (spatiotemporal) overlap; for indirect
overlap via a shared location at different times, see
[`isolate_isolate_sequential_overlap()`](https://theabhirath.github.io/hospitraceR/reference/isolate_isolate_sequential_overlap.md).

## Usage

``` r
isolate_isolate_overlap(isolate_lookup, trace_mat)
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
