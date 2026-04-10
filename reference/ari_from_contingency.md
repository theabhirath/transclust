# Adjusted Rand Index from Contingency Table

Computes the Adjusted Rand Index (ARI) from a contingency table. The ARI
measures the similarity between two clusterings, adjusted for chance.
Values range from -1 to 1, where 1 indicates perfect agreement, 0
indicates random agreement, and negative values indicate less agreement
than expected by chance.

## Usage

``` r
ari_from_contingency(cont_table)
```

## Arguments

- cont_table:

  A contingency table (matrix) where rows correspond to clusters in the
  first assignment and columns to clusters in the second.

## Value

A numeric value representing the Adjusted Rand Index.

## Details

The Adjusted Rand Index is calculated as: \$\$ARI = \frac{Index -
Expected}{MaxIndex - Expected}\$\$

where Index is the sum of \\\binom{n\_{ij}}{2}\\ over all cells,
Expected is derived from the row and column marginals, and MaxIndex is
the average of row and column sum contributions. The values \\n\_{ij}\\
are from the contingency table, \\a_i\\ are row sums, \\b_j\\ are column
sums, and \\n\\ is the total number of isolates.

## References

Hubert, L. and Arabie, P. (1985). Comparing partitions. Journal of
Classification, 2(1), 193-218.
