# Unit test for read_in_seq_aln function
test_that("read_in_seq_aln processes alignment correctly", {
    # Read in sequence file
    dna_aln <- read_in_seq_aln(system.file("extdata", "example.fasta", package = "transclust"))

    expect_s3_class(dna_aln, "DNAbin")
    expect_gt(length(labels(dna_aln)), 0)
    expect_equal(length(unique(labels(dna_aln))), length(labels(dna_aln)),
        info = "Labels should be unique after cleaning"
    )
})
