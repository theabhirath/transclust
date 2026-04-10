# Identify transmission clusters based on the number of shared variants.

Clustering is performed to identify the maximal clusters containing a
single intake-positive patient that occurs before all cluster converts.
The clustering metric is the number of shared variants, and clusters can
have multiple intake-positive patients if they share an identical number
of variants with other cluster members or intake-positive patients occur
after converts. This clustering also requires that clusters be defined
by at least one shared variant that other isolates don't have.

## Usage

``` r
get_tn_clusters_sv_index(
  dna_aln,
  snp_dist,
  adm_seqs,
  adm_pos_pt_seqs,
  seq2pt,
  dates,
  tree
)
```

## Arguments

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- snp_dist:

  A matrix of SNP distances between isolates constructed using a model
  of DNA evolution. See
  [`get_snp_dist_matrix()`](https://theabhirath.github.io/transclust/reference/get_snp_dist_matrix.md)
  for a useful function to generate this.

- adm_seqs:

  A vector of sequence IDs which correspond to admission positive
  patient sequences.

- adm_pos_pt_seqs:

  A vector of all sequence IDs which correspond to admission-positive
  patients, either at intake or collected later. This will be a superset
  of `adm_seqs` by definition.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- dates:

  A vector of isolate dates named by sequence IDs.

- tree:

  A phylogenetic tree object of class `phylo` constructed from the DNA
  alignment. This can be constructed using the
  [`get_phylo_tree()`](https://theabhirath.github.io/transclust/reference/get_phylo_tree.md)
  or can be any other tree object constructed from the same isolates.

## Value

A numeric vector indicating the cluster that each isolate belongs to.

## References

Hawken, S. E., Yelin, R. D., Lolans, K., Pirani, A., Weinstein, R. A.,
Lin, M. Y., Hayden, M. K., & Snitkin, E. S. (2022). Threshold-free
genomic cluster detection to track transmission pathways in health-care
settings: A genomic epidemiology analysis. The Lancet Microbe, 3(9),
e652–e662.
[doi:10.1016/S2666-5247(22)00115-X](https://doi.org/10.1016/S2666-5247%2822%2900115-X)
