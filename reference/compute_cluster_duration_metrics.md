# Compute Cluster Duration Metrics

Calculate temporal metrics for cluster spread

## Usage

``` r
compute_cluster_duration_metrics(seqs, seq2pt, dates)
```

## Arguments

- seqs:

  Vector of sequence IDs in the cluster

- seq2pt:

  Named vector mapping sequence IDs to patient IDs

- dates:

  Vector of isolate dates named by sequence IDs

## Value

List with cluster duration properties and earliest positive date for
each patient
