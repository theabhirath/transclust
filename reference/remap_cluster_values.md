# Remap cluster values to sequential IDs in order of appearance

Renumbers cluster values so that each distinct value becomes a
sequential ID, assigned in order of first appearance. Each occurrence of
`special_val` is given its own separate ID rather than being grouped,
which is useful for values that stand for many independent isolates
(e.g. unclustered isolates) rather than a single cluster.

## Usage

``` r
remap_cluster_values(x, special_val = 0)
```

## Arguments

- x:

  A numeric vector of cluster assignments.

- special_val:

  A value whose every occurrence gets its own ID instead of being
  grouped.

## Value

A numeric vector of remapped cluster assignments.
