load(test_path("testdata", "init.Rdata"))

# Read in sequence file
aln_file <- "KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
dna <- read_in_seq_aln(test_path("testdata", aln_file))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, \(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# Get pairwise distances core
snp_dist <- as.matrix(ape::dist.dna(dna_var, model = "N"))
ref <- which.max(rowMeans(as.matrix(snp_dist)))

# Get the clusters based on a SNP threshold
clusters <- get_tn_clusters_snp_thresh(dna_var, snp_dist, 10)

# test cluster_genetic_context
test_that("cluster_genetic_context works", {
    distances <- cluster_genetic_context(clusters, dna_pt_labels, ip_seqs, snp_dist, "my_analysis")
    # distances should be a matrix.
    expect_true(is.matrix(distances))
})
