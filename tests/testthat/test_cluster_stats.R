library(ape)
load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna_aln <- read.dna(
    system.file("extdata", "example.fasta", package = "transclust"),
    format = "fasta"
)

# Only keep those sequences that are in the trace matrix
dna_aln <- dna_aln[dna_pt_labels[labels(dna_aln)] %in% colnames(trace_mat), ]

# Get all variable positions in the alignment
var_pos <- apply(dna_aln, 2, function(x) sum(x == x[1]) < nrow(dna_aln))

# Only keep variable positions
dna_var <- dna_aln[, var_pos]

# Get pairwise distances
snp_dist <- get_snp_dist_matrix(dna_var)
# Get the clusters based on a SNP threshold
clusters <- get_tn_clusters_snp_thresh(snp_dist, 10)

# test intra_cluster_genetic_var_analysis
test_that("intra_cluster_genetic_var_analysis works", {
    # Get the genetic variation within each cluster
    var_df <- intra_cluster_genetic_var_analysis(clusters, dna_aln, var_pos)
    # var_df should be a data frame.
    expect_true(is.data.frame(var_df))
    # var_df should have the correct number of rows and columns.
    # For nrows: +1 for the overall population but -1 for the singleton, so actually no change.
    expect_equal(nrow(var_df), length(unique(clusters)))
    expect_equal(ncol(var_df), 7)
})

# test cluster_property_perm_test
test_that("cluster_property_perm_test works", {
    nperm <- 10 # 10 permutations is enough for unit tests
    perm_props <- cluster_property_perm_test(
        clusters,
        trace_mat,
        dna_pt_labels,
        ip_seqs,
        ip_seqs_3days,
        dates,
        snp_dist,
        nperm = nperm,
        floor_trace = floor_mat,
        room_trace = room_mat,
        num_cores = 1
    ) # num_cores = 1 for unit tests, because this has to run on CI (and CRAN)
    # perm_props should be an array
    expect_true(is.array(perm_props))
    # nrows should be nperm + 1
    expect_equal(nrow(perm_props), nperm + 1)
})
