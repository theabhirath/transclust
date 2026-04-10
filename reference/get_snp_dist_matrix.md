# Get SNP distance matrix from DNA alignment constructed using a model of DNA evolution.

Get SNP distance matrix from DNA alignment constructed using a model of
DNA evolution.

## Usage

``` r
get_snp_dist_matrix(dna_aln, core = TRUE)
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- core:

  Logical: if `TRUE`, return the core SNP distance matrix, otherwise
  return the full SNP distance matrix.

## Value

A numeric matrix representing the SNP distance between sequences.
