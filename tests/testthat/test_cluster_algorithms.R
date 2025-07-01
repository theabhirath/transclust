library(ape)
load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna_aln <- read.dna(system.file("extdata", "example.fasta", package = "transclust"), format = "fasta")

# Only keep those sequences that are in the trace matrix
dna_aln <- dna_aln[dna_pt_labels[labels(dna_aln)] %in% colnames(trace_mat), ]

# Get all variable positions in the alignment
var_pos <- apply(dna_aln, 2, function(x) sum(x == x[1]) < nrow(dna_aln))

# Only keep variable positions
dna_var <- dna_aln[, var_pos]

# Get pairwise distances
snp_dist <- get_snp_dist_matrix(dna_var)

# test get_tn_clusters_snp_thresh
test_that("get_tn_clusters_snp_thresh works", {
    clusters <- get_tn_clusters_snp_thresh(snp_dist, 10)
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
