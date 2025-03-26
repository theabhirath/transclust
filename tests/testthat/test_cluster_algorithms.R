load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna <- read_in_seq_aln(system.file("extdata", "example.fasta", package = "transclust"))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, function(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# Get pairwise distances
snp_dist <- get_snp_dist_matrix(dna_var)

# test get_tn_clusters_snp_thresh
test_that("get_tn_clusters_snp_thresh works", {
    # Create a hclust object from the distance matrix
    snp_hclust <- hclust(as.dist(snp_dist))
    # Create a phylogenetic tree from the hclust object
    tree <- ape::as.phylo(snp_hclust)
    clusters <- get_tn_clusters_snp_thresh(snp_hclust, tree, 10)
    # clusters should be a vector
    expect_true(is.vector(clusters))
    # clusters should be numeric
    expect_true(is.numeric(clusters))
    # clusters should have the same length as the number of sequences
    expect_equal(length(clusters), nrow(dna_var))
})

# test get_tn_clusters_sv_index
test_that("get_tn_clusters_sv_index works", {
    tree <- get_phylo_tree(dna_var, snp_dist, "pars")
    clusters <- get_tn_clusters_sv_index(
        dna_var, snp_dist, ip_seqs_3days, ip_seqs, dna_pt_labels, dates, tree
    )
    # clusters should be a vector
    expect_true(is.vector(clusters))
    # clusters should be numeric
    expect_true(is.numeric(clusters))
    # clusters should have the same length as the number of sequences
    expect_equal(length(clusters), nrow(dna_var))
})
