# Get phylogenetic tree using neighbor-joining or maximum parsimony method.

This function constructs a neighbor-joining or maximum parsimony
phylogenetic tree from a DNA alignment object and a SNP distance matrix.

## Usage

``` r
get_phylo_tree(dna_aln, snp_dist, method = c("nj", "pars"))
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- snp_dist:

  A numeric matrix representing the SNP distance between sequences. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/transclust/reference/get_snp_dist_matrix.md)
  for a useful function to generate this.

- method:

  A string indicating the method to use for tree construction. Options
  are "nj" (neighbor-joining) or "pars" (maximum parsimony).

## Value

An object of class `phylo` representing the phylogenetic tree.
