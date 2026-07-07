# Build a phylogenetic tree by neighbor-joining or maximum parsimony

Builds a tree from a DNA alignment and its SNP distance matrix. A
neighbor-joining tree is computed first and rooted on the most divergent
isolate (the one with the largest mean SNP distance); with
`method = "pars"` this is then refined into a maximum-parsimony tree.

## Usage

``` r
get_phylo_tree(dna_aln, snp_dist, method = c("nj", "pars"))
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- snp_dist:

  A numeric matrix of SNP distances between sequences. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/hospitraceR/reference/get_snp_dist_matrix.md).

- method:

  Tree-construction method: `"nj"` (neighbor-joining) or `"pars"`
  (maximum parsimony).

## Value

An object of class `phylo` representing the phylogenetic tree.
