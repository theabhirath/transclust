# Fraction of Converts with Same Source Across Cluster Assignments

Compares two cluster assignments by examining patients categorized as
"convert". For each patient who is a convert in both assignments, checks
whether the source patient (categorized as "index" or "weak-index") in
their respective clusters is the same across both assignments.

## Usage

``` r
fraction_convert_same_source(
  clusters1,
  clusters2,
  dna_aln,
  seq2pt,
  adm_seqs,
  dates,
  surv_df
)
```

## Arguments

- clusters1:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- clusters2:

  A named vector of cluster assignments where names are isolate IDs and
  values are cluster numbers.

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- adm_seqs:

  A vector of sequence IDs corresponding to admission positive patients.

- dates:

  A named vector mapping sequence IDs to dates.

- surv_df:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result.

## Value

A numeric value between 0 and 1 representing the fraction of convert
patients whose source (index or weak-index) is the same in both cluster
assignments. Returns NA if there are no common converts.

## Details

The function:

1.  Creates isolate lookups for both cluster assignments

2.  Categorizes patients in both using
    [`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md)

3.  Identifies patients categorized as "convert" in both assignments

4.  For each such convert, finds the "index" or "weak-index" patient in
    their cluster for each assignment

5.  Returns the fraction where the source patient matches

## See also

[`fraction_convert_same_source_from_lookups()`](https://theabhirath.github.io/transclust/reference/fraction_convert_same_source_from_lookups.md),
[`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md),
[`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md)
