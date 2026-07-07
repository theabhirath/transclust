#' Compute intra-cluster SNP-distance summaries
#'
#' @description
#' Returns the mean, median and maximum pairwise SNP distance within a single cluster, taken from
#' the upper triangle of the cluster's own submatrix (excluding the diagonal). Inter-cluster and
#' inter-isolate distances are computed separately, over all clusters at once, by
#' [cluster_inter_distances()].
#'
#' @param cluster_seqs Vector of sequence IDs in the cluster.
#' @param snp_dist Matrix of SNP distances between isolates.
#'
#' @return A named numeric vector with `mean_genetic_distance`, `median_genetic_distance` and
#'   `max_genetic_distance`.
#'
#' @importFrom stats median
#' @export
cluster_pairwise_distances <- function(cluster_seqs, snp_dist) {
    # intra-cluster distances: upper triangle of the cluster's submatrix (excludes the diagonal)
    intra <- snp_dist[cluster_seqs, cluster_seqs, drop = FALSE]
    intra_vals <- intra[upper.tri(intra)]

    c(
        mean_genetic_distance = mean(intra_vals, na.rm = TRUE),
        median_genetic_distance = median(intra_vals, na.rm = TRUE),
        max_genetic_distance = max(intra_vals, na.rm = TRUE)
    )
}

#' Compute per-cluster nearest-neighbor SNP distances
#'
#' @description
#' For every non-singleton cluster, the minimum SNP distance to an isolate in another cluster
#' (`min_inter_cluster`) and to any isolate outside the cluster (`min_inter_isolate`). Unlike
#' [cluster_pairwise_distances()], this works over the whole clustering at once: it derives each
#' cluster's "other cluster" and "non cluster" comparison sets from the full assignment vector.
#' Singleton clusters are dropped via [remove_singleton_clusters()], so their sequences never form
#' a cluster of their own but still count as unclustered isolates toward `min_inter_isolate`.
#'
#' @param clusters A vector named by sequence IDs giving the cluster each sequence belongs to.
#' @param snp_dist A matrix of SNP distances between isolates. Its row/column names define the full
#'   universe of isolates, including those not assigned to any cluster.
#'
#' @return A matrix with one row per cluster (named by cluster ID) and columns `min_inter_cluster`
#'   and `min_inter_isolate`.
#'
#' @export
cluster_inter_distances <- function(clusters, snp_dist) {
    # The full universe of isolates, including any not assigned to a (multi-isolate) cluster
    all_seqs <- rownames(snp_dist)
    # Drop singleton clusters so their sequences count as unclustered isolates
    clusters <- remove_singleton_clusters(clusters)
    unique_clusters <- sort(unique(clusters))

    inter <- vapply(
        unique_clusters,
        function(cl) {
            cluster_seqs <- names(clusters[clusters == cl])
            # Isolates assigned to a different cluster
            other_cluster_seqs <- names(clusters[clusters != cl])
            # Any isolate not in this cluster, including unclustered (singleton) isolates
            non_cluster_seqs <- setdiff(all_seqs, cluster_seqs)

            c(
                min_inter_cluster = min(
                    snp_dist[cluster_seqs, other_cluster_seqs, drop = FALSE],
                    na.rm = TRUE
                ),
                min_inter_isolate = min(
                    snp_dist[cluster_seqs, non_cluster_seqs, drop = FALSE],
                    na.rm = TRUE
                )
            )
        },
        numeric(2)
    )

    inter_mat <- t(inter)
    rownames(inter_mat) <- unique_clusters
    inter_mat
}

# Stack a per-cluster reduction over every multi-patient cluster in an isolate
# lookup. `metric_fn` receives the lookup rows for one cluster and returns a
# named numeric vector; the results become a data frame with a trailing
# `cluster` column. Single-patient clusters are excluded throughout.
# @keyword internal
per_cluster_metric_df <- function(isolate_lookup, metric_fn) {
    cluster_ids <- get_non_single_patient_clusters(isolate_lookup)
    rows <- lapply(cluster_ids, function(cluster) {
        metric_fn(isolate_lookup[isolate_lookup$cluster == cluster, , drop = FALSE])
    })
    out <- as.data.frame(do.call(rbind, rows))
    out$cluster <- cluster_ids
    out
}

#' Compute genetic distances by cluster
#'
#' @description
#' Summarize intra-cluster SNP distances (mean/median/max) for every multi-patient
#' cluster in an isolate lookup, via [cluster_pairwise_distances()].
#'
#' @param isolate_lookup An isolate lookup table from [get_isolate_lookup()].
#' @param snp_dist Matrix of SNP distances between isolates.
#'
#' @return Data frame with `mean_genetic_distance`, `median_genetic_distance`,
#'   `max_genetic_distance` and `cluster`, one row per cluster.
#' @export
compute_genetic_distances_by_cluster <- function(isolate_lookup, snp_dist) {
    per_cluster_metric_df(isolate_lookup, function(cluster_lookup) {
        cluster_pairwise_distances(cluster_lookup$isolate_id, snp_dist)
    })
}

#' Compute intra-cluster duration metrics
#'
#' @description
#' Calculate temporal metrics for spread of a single cluster.
#'
#' @param seqs Vector of sequence IDs in the cluster
#' @param seq2pt Named vector mapping sequence IDs to patient IDs
#' @param dates Vector of isolate dates named by sequence IDs
#'
#' @return Named numeric vector with `cluster_start_date`,
#'   `time_to_first_acquisition`, `time_to_last_acquisition` and
#'   `median_time_to_acquisition`.
#'
#' @importFrom stats median
#' @export
intra_cluster_duration_metrics <- function(seqs, seq2pt, dates) {
    # earliest positive date for each patient in the cluster, sorted ascending
    earliest_by_pt <- sort(vapply(
        unique(seq2pt[seqs]),
        function(pt_id) min(dates[seqs[seq2pt[seqs] == pt_id]]),
        numeric(1)
    ))

    n <- length(earliest_by_pt)
    if (n == 0) {
        return(c(
            cluster_start_date = NA_real_,
            time_to_first_acquisition = NA_real_,
            time_to_last_acquisition = NA_real_,
            median_time_to_acquisition = NA_real_
        ))
    }

    # acquisition gaps relative to the first patient's date (empty when n == 1)
    since_start <- earliest_by_pt[-1] - earliest_by_pt[[1]]
    c(
        cluster_start_date = earliest_by_pt[[1]],
        time_to_first_acquisition = if (n >= 2) since_start[[1]] else 0,
        time_to_last_acquisition = earliest_by_pt[[n]] - earliest_by_pt[[1]],
        median_time_to_acquisition = if (n >= 2) median(since_start) else 0
    )
}

#' Compute duration metrics by cluster
#'
#' @description
#' Calculate temporal metrics for cluster spread for a set of cluster assignments.
#'
#' @param isolate_lookup An isolate lookup table from [get_isolate_lookup()].
#'
#' @return Data frame with duration metrics by cluster
#'
#' @importFrom stats setNames
#' @export
compute_duration_metrics_by_cluster <- function(isolate_lookup) {
    per_cluster_metric_df(isolate_lookup, function(cluster_lookup) {
        seqs <- cluster_lookup$isolate_id
        intra_cluster_duration_metrics(
            seqs,
            setNames(cluster_lookup$patient_id, seqs),
            setNames(cluster_lookup$date, seqs)
        )
    })
}
