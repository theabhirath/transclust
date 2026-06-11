# Tests for the pure cluster-bookkeeping helpers with small synthetic examples

test_that("remove_singleton_clusters drops only single-sequence clusters", {
    clusters <- c(a = 1, b = 1, c = 2, d = 3, e = 3)
    out <- remove_singleton_clusters(clusters)

    # cluster 2 is the only singleton and should be removed; names are preserved
    expect_equal(out, c(a = 1, b = 1, d = 3, e = 3))
})

test_that("remove_singleton_clusters returns empty when all clusters are singletons", {
    clusters <- c(a = 1, b = 2, c = 3)
    expect_length(remove_singleton_clusters(clusters), 0)
})

test_that("get_non_single_patient_clusters keeps clusters with >1 distinct patient", {
    lookup <- data.frame(
        cluster = c(1, 1, 2, 2, 3),
        patient_id = c("p1", "p2", "p3", "p3", "p4"),
        stringsAsFactors = FALSE
    )

    # cluster 1 has two patients; cluster 2 has one patient (sampled twice);
    # cluster 3 is a single isolate -> only cluster 1 qualifies
    expect_equal(get_non_single_patient_clusters(lookup), 1)
})

test_that("flatten_cluster_patient_categorization yields one row per cluster-patient pair", {
    categorization <- list(
        "1" = c(p1 = "index", p2 = "convert"),
        "2" = c(p3 = "index")
    )

    df <- flatten_cluster_patient_categorization(categorization)

    expect_equal(names(df), c("cluster", "patient_id", "category"))
    expect_equal(nrow(df), 3)
    expect_equal(df$cluster, c("1", "1", "2"))
    expect_equal(df$patient_id, c("p1", "p2", "p3"))
    expect_equal(df$category, c("index", "convert", "index"))
})

test_that("flatten_cluster_overlap_categorization yields one row per cluster", {
    categorization <- c("1" = "patient-to-patient", "2" = "inexplicable")

    df <- flatten_cluster_overlap_categorization(categorization)

    expect_equal(names(df), c("cluster", "category"))
    expect_equal(df$cluster, c("1", "2"))
    expect_equal(df$category, c("patient-to-patient", "inexplicable"))
})
