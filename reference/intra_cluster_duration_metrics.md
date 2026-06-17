# Compute intra-cluster duration metrics

Calculate temporal metrics for spread of a single cluster.

## Usage

``` r
intra_cluster_duration_metrics(seqs, seq2pt, dates)
```

## Arguments

- seqs:

  Vector of sequence IDs in the cluster

- seq2pt:

  Named vector mapping sequence IDs to patient IDs

- dates:

  Vector of isolate dates named by sequence IDs

## Value

Named numeric vector with `cluster_start_date`,
`time_to_first_acquisition`, `time_to_last_acquisition` and
`median_time_to_acquisition`.
