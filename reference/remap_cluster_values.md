# Remap cluster values: each unique value (including each special value) gets its own number in order of appearance.

This function remaps cluster values so that each unique value (including
each special value) gets its own number in order of appearance.

## Usage

``` r
remap_cluster_values(x, special_val = 0)
```

## Arguments

- x:

  A numeric vector of cluster assignments.

- special_val:

  A value to treat as special i.e. it will get its own number in order
  of appearance. This is useful for values that are not part of the
  cluster assignment, such as singleton clusters.

## Value

A numeric vector of remapped cluster assignments.
