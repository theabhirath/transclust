# load(test_path("testdata", "init.Rdata"))

# # Read in sequence file
# aln_file <- "KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
# dna <- read_in_seq_aln(test_path("testdata", aln_file))

# Get variable positions in the alignment
dna_var_pos <- apply(dna, 2, function(x) sum(x == x[1]) < nrow(dna))
dna_var <- dna[dna_pt_labels[labels(dna)] %in% colnames(trace_mat), dna_var_pos]

# Get pairwise distances core
snp_dist <- as.matrix(ape::dist.dna(dna_var, model = "N"))
ref <- which.max(rowMeans(as.matrix(snp_dist)))

clusters <- get_tn_clusters_MSV_SVst_index_first(
    dna_var, snp_dist, ip_seqs_3days, ip_seqs, dna_pt_labels, dates,
    pars = TRUE
)
