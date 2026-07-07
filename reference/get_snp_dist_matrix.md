# Compute a pairwise SNP distance matrix from a DNA alignment

Counts the number of nucleotide differences between every pair of
sequences, via
[`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html) with
model `"N"`.

## Usage

``` r
get_snp_dist_matrix(dna_aln, core = TRUE)
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- core:

  Logical; if `TRUE`, sites with missing data in any sequence are
  dropped before counting, giving a core-genome distance. If `FALSE`,
  missing data are handled per pair of sequences (pairwise deletion).

## Value

A numeric matrix of pairwise SNP distances between sequences.
