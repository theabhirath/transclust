# Get clusters containing more than one patient

Get clusters containing more than one patient

## Usage

``` r
get_non_single_patient_clusters(isolate_lookup)
```

## Arguments

- isolate_lookup:

  A lookup table for isolates and their cluster assignments with other
  relevant epidemiological information. See
  [`get_isolate_lookup()`](https://theabhirath.github.io/hospitraceR/reference/get_isolate_lookup.md).

## Value

A numeric vector of cluster IDs that contain isolates from more than one
patient.
