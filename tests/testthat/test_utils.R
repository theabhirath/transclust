library(ape)
load(system.file("extdata", "example.Rdata", package = "transclust"))

# Read in sequence file
dna_aln <- read.dna(
    system.file("extdata", "example.fasta", package = "transclust"),
    format = "fasta"
)

# Only keep those sequences that are in the trace matrix
dna_aln <- dna_aln[dna_pt_labels[labels(dna_aln)] %in% colnames(trace_mat), ]

# Get variable positions in the alignment
dna_var <- dna_aln[, apply(dna_aln, 2, function(x) sum(x == x[1]) < nrow(dna_aln))]

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
