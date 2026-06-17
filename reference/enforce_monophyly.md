# Force cluster assignments to be monophyletic with respect to a tree.

Checks each multi-isolate cluster for monophyly on `tree` and reconciles
those that are not monophyletic using the chosen strategy. The result is
renumbered sequentially with
[`remap_cluster_values()`](https://theabhirath.github.io/hospitraceR/reference/remap_cluster_values.md),
and each reconciled cluster is reported via
[`message()`](https://rdrr.io/r/base/message.html) against its final
(remapped) cluster id, listing the isolates that make it up.

## Usage

``` r
enforce_monophyly(
  clusters,
  tree,
  monophyly_method = c("expand", "break_down"),
  special_val = NULL
)
```

## Arguments

- clusters:

  A named numeric vector of cluster assignments; the names are the
  isolate/tip labels used to locate members on the tree.

- tree:

  A phylogenetic tree of class `phylo`.

- monophyly_method:

  How to make a non-monophyletic cluster monophyletic. `"expand"` grows
  the cluster to the smallest clade of the tree that contains all of its
  members, absorbing any intervening isolates. `"break_down"` instead
  splits the cluster into the largest monophyletic clades it already
  contains, without pulling in foreign isolates.

- special_val:

  An optional cluster value to leave untouched, e.g. `0` for unclustered
  isolates that should not be treated as a single cluster. Default is
  `NULL`.

## Value

A numeric vector of cluster assignments that are monophyletic with
respect to the tree, renumbered sequentially.
