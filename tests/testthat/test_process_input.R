# Unit test for read_in_seq_aln function
test_that("read_in_seq_aln processes alignment correctly", {
    aln_file <- "../../inst/extdata/KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta"
    dna <- read_in_seq_aln(aln_file)

    expect_s3_class(dna, "DNAbin")
    expect_gt(length(labels(dna)), 0)
    expect_equal(length(unique(labels(dna))), length(labels(dna)),
        info = "Labels should be unique after cleaning"
    )
})
