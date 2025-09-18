#' Calculate Trace Overlap
#'
#' @description
#' Calculate the overlap between isolates in a trace matrix
#'
#' @param trace_mat A matrix with rows representing days and columns representing isolates
#'
#' @return A matrix with the overlap between isolates
#'
#' @export
calculate_trace_overlap <- function(trace_mat) {
    isolates <- row.names(trace_mat)
    overlap_mat <- matrix(0, nrow = length(isolates), ncol = length(isolates))
    rownames(overlap_mat) <- isolates
    colnames(overlap_mat) <- isolates

    for (i in 1:(length(isolates) - 1)) {
        for (j in (i + 1):length(isolates)) {
            iso1 <- isolates[i]
            iso2 <- isolates[j]

            shared_days <- sum(
                trace_mat[, iso1] == trace_mat[, iso2] &
                    !is.na(trace_mat[, iso1]) &
                    !is.na(trace_mat[, iso2])
            )
            overlap_mat[iso1, iso2] <- shared_days
            overlap_mat[iso2, iso1] <- shared_days
        }
    }
    overlap_mat
}

#' Analyze Overlap by Cluster
#'
#' @description
#' Analyze the overlap between isolates in a trace matrix by cluster
#'
#' @param clusters A named numeric vector where names are sequence IDs and values are subtrees defining the cluster.
#' @param trace_mat A matrix with rows representing days and columns representing isolates
#'
#' @return A data frame with the overlap between isolates
#'
#' @importFrom dplyr filter mutate
#' @export
analyze_overlap_by_cluster <- function(clusters, trace_mat) {
    overlap_mat <- calculate_trace_overlap(trace_mat)
    cluster_seqs <- names(clusters)
    cluster_df <- expand.grid(i = cluster_seqs, j = cluster_seqs, stringsAsFactors = FALSE) |>
        filter(i < j) |>
        filter(i %in% row.names(overlap_mat), j %in% colnames(overlap_mat)) |>
        mutate(
            same_cluster = clusters[i] == clusters[j],
            overlap = mapply(function(x, y) overlap_mat[x, y], i, j)
        )
    cluster_df
}
