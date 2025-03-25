load(test_path("testdata", "init.Rdata"))

# Read in sequence file
aln_file <- "KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
dna <- read_in_seq_aln(test_path("testdata", aln_file))

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

# test cluster_genetic_context
test_that("cluster_genetic_context works", {
    distances <- cluster_genetic_context(clusters, dna_pt_labels, ip_seqs, snp_dist, "my_analysis")
    # distances should be a matrix.
    expect_true(is.matrix(distances))
})
