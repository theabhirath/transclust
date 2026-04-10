# Create a lookup table for isolates along with their cluster assignments.

This function creates a lookup table for isolates. The columns are:

- cluster: the cluster assignment for the isolate.

- isolate_id: the isolate ID.

- patient_id: the patient ID.

- date: the date of the isolate.

- adm_pos: a logical value indicating if the isolate is admission
  positive.

- prev_surv: the previous surveillance date for the patient.

## Usage

``` r
get_isolate_lookup(clusters, dna_aln, seq2pt, adm_seqs, dates, surv_df)
```

## Arguments

- clusters:

  A numeric vector of cluster assignments.

- dna_aln:

  A DNA alignment object of class `DNAbin`.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- adm_seqs:

  A vector of sequence IDs which correspond to admission positive
  patient sequences.

- dates:

  A named vector mapping sequence IDs to dates.

- surv_df:

  A data frame with surveillance information for isolates.

## Value

A data frame with columns: isolate_id, patient_id, date.
