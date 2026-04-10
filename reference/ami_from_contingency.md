# Adjusted Mutual Information from Contingency Table

Computes the Adjusted Mutual Information (AMI) from a contingency table.
The AMI measures the mutual information between two clusterings,
adjusted for chance. Values range from 0 to 1, where 1 indicates perfect
agreement and 0 indicates independence.

## Usage

``` r
ami_from_contingency(cont_table)
```

## Arguments

- cont_table:

  A contingency table (matrix) where rows correspond to clusters in the
  first assignment and columns to clusters in the second.

## Value

A numeric value representing the Adjusted Mutual Information.

## Details

The Adjusted Mutual Information is calculated as: \$\$AMI = \frac{MI -
E\[MI\]}{\max(H(U), H(V)) - E\[MI\]}\$\$

where \\MI\\ is the mutual information, \\E\[MI\]\\ is the expected
mutual information under random permutation, and \\H(U)\\ and \\H(V)\\
are the entropies of the two clusterings.

## References

Vinh, N. X., Epps, J., and Bailey, J. (2010). Information theoretic
measures for clusterings comparison: Variants, properties, normalization
and correction for chance. Journal of Machine Learning Research, 11,
2837-2854.
