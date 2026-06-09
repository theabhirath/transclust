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

Vector with cluster duration properties
