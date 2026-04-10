# Package index

## Transmission cluster detection algorithms

A collection of algorithms for detecting transmission clusters in
epidemiological data.

- [`get_tn_clusters_snp_thresh()`](https://theabhirath.github.io/transclust/reference/get_tn_clusters_snp_thresh.md)
  : Perform clustering of isolates using a hard SNP distance cutoff.
- [`get_tn_clusters_sv_index()`](https://theabhirath.github.io/transclust/reference/get_tn_clusters_sv_index.md)
  : Identify transmission clusters based on the number of shared
  variants.

## Statistics

Functions for statistical analysis of transmission clustering results.

- [`intra_cluster_genetic_var_analysis()`](https://theabhirath.github.io/transclust/reference/intra_cluster_genetic_var_analysis.md)
  : Analyze Intra-Cluster Genetic Variation
- [`intra_cluster_genetic_distances()`](https://theabhirath.github.io/transclust/reference/intra_cluster_genetic_distances.md)
  : Compute Genetic Distances
- [`compute_cluster_duration_metrics()`](https://theabhirath.github.io/transclust/reference/compute_cluster_duration_metrics.md)
  : Compute Cluster Duration Metrics
- [`cluster_overlap_perm_test()`](https://theabhirath.github.io/transclust/reference/cluster_overlap_perm_test.md)
  : Permutation test for cluster overlap fractions

## Cluster comparison

Functions for comparing different cluster assignments using contingency
tables, information-theoretic measures, and source concordance.

- [`cluster_contingency_table()`](https://theabhirath.github.io/transclust/reference/cluster_contingency_table.md)
  : Create Contingency Table from Two Cluster Assignments
- [`adjusted_mutual_information()`](https://theabhirath.github.io/transclust/reference/adjusted_mutual_information.md)
  : Adjusted Mutual Information from Cluster Assignments
- [`ami_from_contingency()`](https://theabhirath.github.io/transclust/reference/ami_from_contingency.md)
  : Adjusted Mutual Information from Contingency Table
- [`adjusted_rand_index()`](https://theabhirath.github.io/transclust/reference/adjusted_rand_index.md)
  : Adjusted Rand Index from Cluster Assignments
- [`ari_from_contingency()`](https://theabhirath.github.io/transclust/reference/ari_from_contingency.md)
  : Adjusted Rand Index from Contingency Table
- [`fraction_convert_same_source()`](https://theabhirath.github.io/transclust/reference/fraction_convert_same_source.md)
  : Fraction of Converts with Same Source Across Cluster Assignments
- [`fraction_convert_same_source_from_lookups()`](https://theabhirath.github.io/transclust/reference/fraction_convert_same_source_from_lookups.md)
  : Fraction of Converts with Same Source from Isolate Lookups

## Epidemiological overlap analysis

Functions for analyzing patient overlap within and between transmission
clusters, and categorizing clusters by epidemiological patterns.

- [`isolate_isolate_overlap()`](https://theabhirath.github.io/transclust/reference/isolate_isolate_overlap.md)
  : Calculate isolate-isolate overlap
- [`isolate_isolate_sequential_overlap()`](https://theabhirath.github.io/transclust/reference/isolate_isolate_sequential_overlap.md)
  : Calculate isolate-isolate sequential overlap
- [`cluster_isolate_overlap()`](https://theabhirath.github.io/transclust/reference/cluster_isolate_overlap.md)
  : Calculate isolate-isolate overlap for clusters
- [`fraction_converts_with_overlap()`](https://theabhirath.github.io/transclust/reference/fraction_converts_with_overlap.md)
  : Calculate the fraction of converts with overlap
- [`categorize_cluster_overlap()`](https://theabhirath.github.io/transclust/reference/categorize_cluster_overlap.md)
  : Categorize cluster overlap
- [`cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/cluster_patient_categorization.md)
  : Categorize patients within clusters by epidemiological status
- [`flatten_cluster_patient_categorization()`](https://theabhirath.github.io/transclust/reference/flatten_cluster_patient_categorization.md)
  : Flatten cluster-patient categorization to a data frame

## Utilities

A collection of utility functions for this package.

- [`get_snp_dist_matrix()`](https://theabhirath.github.io/transclust/reference/get_snp_dist_matrix.md)
  : Get SNP distance matrix from DNA alignment constructed using a model
  of DNA evolution.
- [`get_phylo_tree()`](https://theabhirath.github.io/transclust/reference/get_phylo_tree.md)
  : Get phylogenetic tree using neighbor-joining or maximum parsimony
  method.
- [`get_isolate_lookup()`](https://theabhirath.github.io/transclust/reference/get_isolate_lookup.md)
  : Create a lookup table for isolates along with their cluster
  assignments.
- [`get_non_single_patient_clusters()`](https://theabhirath.github.io/transclust/reference/get_non_single_patient_clusters.md)
  : Get non-single patient clusters from a vector of cluster
  assignments.
- [`remove_singleton_clusters()`](https://theabhirath.github.io/transclust/reference/remove_singleton_clusters.md)
  : Remove singleton clusters from a vector of cluster assignments.
