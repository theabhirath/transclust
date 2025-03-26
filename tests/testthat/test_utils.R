load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna <- read_in_seq_aln(system.file("extdata", "example.fasta", package = "transclust"))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, function(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# test get_snp_dist_matrix
test_that("get_snp_dist_matrix works", {
    # Get the SNP distance matrix
    snp_dist <- get_snp_dist_matrix(dna_var)
    # Check that the function returns a matrix
    expect_true(is.matrix(snp_dist))
    # Check that the matrix is square
    expect_equal(nrow(snp_dist), ncol(snp_dist))
    # Check that the matrix has the same number of rows and columns as the input data
    expect_equal(nrow(snp_dist), nrow(dna_var))
})

# test get_phylo_tree
test_that("get_phylo_tree works", {
    # Get the SNP distance matrix
    snp_dist <- get_snp_dist_matrix(dna_var)

    # Neighbor-Joining tree
    tree <- get_phylo_tree(dna_var, snp_dist, "nj")
    # Check that the function returns a phylogenetic tree
    expect_true(inherits(tree, "phylo"))
    # Check that the tree has the same number of tips as the input data
    expect_equal(length(tree$tip.label), nrow(dna_var))

    # Parsimony tree
    tree <- get_phylo_tree(dna_var, snp_dist, "pars")
    # Check that the function returns a phylogenetic tree
    expect_true(inherits(tree, "phylo"))
    # Check that the tree has the same number of tips as the input data
    expect_equal(length(tree$tip.label), nrow(dna_var))
})
