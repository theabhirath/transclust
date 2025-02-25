load(test_path("testdata", "init.RData"))

# Read in sequence file
aln_file <- "KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
dna <- read_in_seq_aln(test_path("testdata", aln_file))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, \(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# Get pairwise distances core
snp_dist <- as.matrix(ape::dist.dna(dna_var, model = "N"))
ref <- which.max(rowMeans(as.matrix(snp_dist)))

# test get_tn_clusters_snp_thresh
test_that("get_tn_clusters_snp_thresh works", {
    clusters <- get_tn_clusters_snp_thresh(dna_var, snp_dist, 10)
    # clusters should be a vector.
    expect_true(is.vector(clusters))
})
