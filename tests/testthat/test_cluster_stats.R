load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna <- read_in_seq_aln(system.file("extdata", "example.fasta", package = "transclust"))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, function(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# Get pairwise distances
snp_dist <- get_snp_dist_matrix(dna_var)
# Create a hclust object from the distance matrix
snp_hclust <- hclust(as.dist(snp_dist))
# Create a phylogenetic tree from the hclust object
tree <- ape::as.phylo(snp_hclust)

# Get the clusters based on a SNP threshold
clusters <- get_tn_clusters_snp_thresh(snp_hclust, tree, 10)

# test intra_cluster_genetic_var_analysis
test_that("intra_cluster_genetic_var_analysis works", {
    # Get the genetic variation within each cluster
    var_df <- intra_cluster_genetic_var_analysis(clusters, dna_var)
    # var_df should be a data frame.
    expect_true(is.data.frame(var_df))
    # var_df should have the correct number of rows and columns.
    # For nrows: +1 for the overall population but -1 for the singleton, so actually no change.
    expect_equal(nrow(var_df), length(unique(clusters)))
    expect_equal(ncol(var_df), 7)
})
