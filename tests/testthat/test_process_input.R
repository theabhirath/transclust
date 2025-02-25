# Unit test for read_in_seq_aln function
test_that("read_in_seq_aln processes alignment correctly", {
    # Read in sequence file
    aln_file <- "KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
    dna <- read_in_seq_aln(test_path("testdata", aln_file))

    expect_s3_class(dna, "DNAbin")
    expect_gt(length(labels(dna)), 0)
    expect_equal(length(unique(labels(dna))), length(labels(dna)),
        info = "Labels should be unique after cleaning"
    )
})
