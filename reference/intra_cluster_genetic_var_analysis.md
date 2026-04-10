# Analyze Intra-Cluster Genetic Variation

This function analyzes genetic diversity within clusters by identifying
variant sites and calculating mutational spectra (e.g. transition and
transversion rates). It uses a DNA alignment and a vector of cluster
assignments (named by sequence IDs) and returns a data frame summarizing
the diversity measures for each cluster as well as for the overall
population.

## Usage

``` r
intra_cluster_genetic_var_analysis(clusters, dna_aln)
```

## Arguments

- clusters:

  A named vector of cluster assignments (cluster IDs) with names
  corresponding to sequence IDs.

- dna_aln:

  A matrix representing the DNA sequence alignment, where rows are
  sequences (named by sequence IDs) and columns represent nucleotide
  positions.

## Value

A data frame with one row per cluster (plus a row for the overall
population) containing the number of variable sites and the rates for
six mutation types. Mutation rates are rounded to three decimals.
