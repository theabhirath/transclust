# Tests for the cluster-comparison metrics with small synthetic examples

# --- cluster_contingency_table -------------------------------------------------

test_that("cluster_contingency_table counts shared isolates with sorted labels", {
    c1 <- c(a = 1, b = 1, c = 2, d = 2)
    c2 <- c(a = 1, b = 2, c = 1, d = 2)

    tab <- cluster_contingency_table(c1, c2)

    expect_equal(dim(tab), c(2, 2))
    expect_equal(rownames(tab), c("1", "2"))
    expect_equal(colnames(tab), c("1", "2"))
    # a->(1,1), b->(1,2), c->(2,1), d->(2,2): one isolate in every cell
    expect_equal(unname(tab), matrix(1L, 2, 2))
    # every isolate is accounted for exactly once
    expect_equal(sum(tab), length(c1))
})

test_that("cluster_contingency_table subsets to common isolates with a warning", {
    c1 <- c(a = 1, b = 1, c = 2)
    c2 <- c(a = 1, b = 1, d = 2)

    expect_warning(tab <- cluster_contingency_table(c1, c2), "common isolates")
    # only a and b are shared, both in cluster 1 of each assignment
    expect_equal(sum(tab), 2)
    expect_equal(unname(tab), matrix(2L, 1, 1))
})

test_that("cluster_contingency_table errors when no isolates are shared", {
    # the mismatch warning fires first, then the error; suppress the former
    expect_error(
        suppressWarnings(cluster_contingency_table(c(a = 1), c(b = 1))),
        "No common isolates"
    )
})

# --- Adjusted Rand Index -------------------------------------------------------

test_that("ari_from_contingency returns 1 for identical clusterings", {
    cl <- c(a = 1, b = 1, c = 2, d = 2, e = 3)
    expect_equal(adjusted_rand_index(cl, cl), 1)
})

test_that("ari_from_contingency matches a hand-computed value", {
    # contingency table [[2,1],[1,2]] over n = 6:
    #   sum C(nij,2) = 1+0+0+1 = 2; row/col sums all 3 -> sum C(.,2) = 6 each
    #   expected = 6*6/15 = 2.4; max = 6; ARI = (2 - 2.4)/(6 - 2.4) = -1/9
    tab <- matrix(c(2, 1, 1, 2), nrow = 2)
    expect_equal(ari_from_contingency(tab), -1 / 9)
})

test_that("adjusted_rand_index is symmetric and label-invariant", {
    c1 <- c(a = 1, b = 1, c = 2, d = 2, e = 3, f = 3)
    c2 <- c(a = 5, b = 5, c = 7, d = 7, e = 1, f = 2)

    expect_equal(adjusted_rand_index(c1, c2), adjusted_rand_index(c2, c1))
    # relabeling c1's clusters must not change the index
    relabel <- c("1" = 9, "2" = 4, "3" = 8)
    c1_relabeled <- relabel[as.character(c1)]
    names(c1_relabeled) <- names(c1)
    expect_equal(adjusted_rand_index(c1, c2), adjusted_rand_index(c1_relabeled, c2))
})

test_that("ari_from_contingency returns NA when fewer than two isolates", {
    expect_true(is.na(ari_from_contingency(matrix(1L, 1, 1))))
})

# --- Adjusted Mutual Information ----------------------------------------------

test_that("adjusted_mutual_information returns 1 for identical clusterings", {
    cl <- c(a = 1, b = 1, c = 2, d = 2, e = 3)
    expect_equal(adjusted_mutual_information(cl, cl), 1)
})

test_that("adjusted_mutual_information is symmetric and bounded above by 1", {
    c1 <- c(a = 1, b = 1, c = 2, d = 2, e = 3, f = 3)
    c2 <- c(a = 1, b = 2, c = 2, d = 3, e = 3, f = 1)

    ami <- adjusted_mutual_information(c1, c2)
    expect_equal(ami, adjusted_mutual_information(c2, c1))
    expect_lte(ami, 1)
})

test_that("wrappers agree with the contingency-table forms", {
    c1 <- c(a = 1, b = 1, c = 2, d = 2, e = 3, f = 3)
    c2 <- c(a = 1, b = 2, c = 2, d = 3, e = 3, f = 1)
    tab <- cluster_contingency_table(c1, c2)

    expect_equal(adjusted_rand_index(c1, c2), ari_from_contingency(tab))
    expect_equal(adjusted_mutual_information(c1, c2), ami_from_contingency(tab))
})
