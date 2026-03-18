#' Create Contingency Table from Two Cluster Assignments
#'
#' @description
#' Generates a contingency table comparing two cluster assignments. The table
#' counts the number of isolates shared between each pair of clusters from
#' the two assignments. Only isolates present in both assignments are considered.
#'
#' @param clusters1 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#' @param clusters2 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#'
#' @return A matrix where rows correspond to clusters in `clusters1`,
#'   columns correspond to clusters in `clusters2`, and cell values
#'   represent the count of shared isolates.
#'
#' @export
cluster_contingency_table <- function(clusters1, clusters2) {
    # find isolates common to both assignments
    common_isolates <- intersect(names(clusters1), names(clusters2))

    if (length(common_isolates) == 0) {
        stop("No common isolates found between the two cluster assignments")
    }

    # subset to common isolates
    c1 <- clusters1[common_isolates]
    c2 <- clusters2[common_isolates]

    # get unique cluster labels
    labels1 <- sort(unique(c1))
    labels2 <- sort(unique(c2))

    # build contingency table
    cont_table <- matrix(
        0L,
        nrow = length(labels1),
        ncol = length(labels2),
        dimnames = list(labels1, labels2)
    )

    for (isolate in common_isolates) {
        row_idx <- as.character(c1[isolate])
        col_idx <- as.character(c2[isolate])
        cont_table[row_idx, col_idx] <- cont_table[row_idx, col_idx] + 1L
    }

    cont_table
}

#' Adjusted Rand Index from Contingency Table
#'
#' @description
#' Computes the Adjusted Rand Index (ARI) from a contingency table. The ARI
#' measures the similarity between two clusterings, adjusted for chance.
#' Values range from -1 to 1, where 1 indicates perfect agreement, 0 indicates
#' random agreement, and negative values indicate less agreement than expected
#' by chance.
#'
#' @param cont_table A contingency table (matrix) where rows correspond to
#'   clusters in the first assignment and columns to clusters in the second.
#'
#' @return A numeric value representing the Adjusted Rand Index.
#'
#' @details
#' The Adjusted Rand Index is calculated as:
#' \deqn{ARI = \frac{Index - Expected}{MaxIndex - Expected}}
#'
#' where Index is the sum of \eqn{\binom{n_{ij}}{2}} over all cells,
#' Expected is derived from the row and column marginals, and MaxIndex
#' is the average of row and column sum contributions. The values
#' \eqn{n_{ij}} are from the contingency table, \eqn{a_i} are row sums,
#' \eqn{b_j} are column sums, and \eqn{n} is the total number of isolates.
#'
#' @references
#' Hubert, L. and Arabie, P. (1985). Comparing partitions. Journal of
#' Classification, 2(1), 193-218.
#'
#' @export
ari_from_contingency <- function(cont_table) {
    # total number of isolates
    n <- sum(cont_table)

    if (n < 2) {
        return(NA_real_)
    }

    # row and column sums
    a <- rowSums(cont_table)
    b <- colSums(cont_table)

    # sum of C(n_ij, 2) over all cells
    sum_nij_choose2 <- sum(choose(cont_table, 2))

    # sum of C(a_i, 2) and C(b_j, 2)
    sum_a_choose2 <- sum(choose(a, 2))
    sum_b_choose2 <- sum(choose(b, 2))

    # total combinations
    n_choose2 <- choose(n, 2)

    # expected index
    expected_index <- (sum_a_choose2 * sum_b_choose2) / n_choose2

    # max index
    max_index <- (sum_a_choose2 + sum_b_choose2) / 2

    # ARI formula
    numerator <- sum_nij_choose2 - expected_index
    denominator <- max_index - expected_index

    if (denominator == 0) {
        # all isolates in one cluster for both assignments
        return(1)
    }

    numerator / denominator
}

#' Adjusted Rand Index from Cluster Assignments
#'
#' @description
#' Computes the Adjusted Rand Index (ARI) between two cluster assignments.
#' This is a convenience wrapper around [ari_from_contingency()].
#'
#' @param clusters1 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#' @param clusters2 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#'
#' @return A numeric value representing the Adjusted Rand Index.
#'
#' @seealso [ari_from_contingency()], [cluster_contingency_table()]
#'
#' @export
adjusted_rand_index <- function(clusters1, clusters2) {
    cont_table <- cluster_contingency_table(clusters1, clusters2)
    ari_from_contingency(cont_table)
}

#' Adjusted Mutual Information from Contingency Table
#'
#' @description
#' Computes the Adjusted Mutual Information (AMI) from a contingency table.
#' The AMI measures the mutual information between two clusterings, adjusted
#' for chance. Values range from 0 to 1, where 1 indicates perfect agreement
#' and 0 indicates independence.
#'
#' @param cont_table A contingency table (matrix) where rows correspond to
#'   clusters in the first assignment and columns to clusters in the second.
#'
#' @return A numeric value representing the Adjusted Mutual Information.
#'
#' @details
#' The Adjusted Mutual Information is calculated as:
#' \deqn{AMI = \frac{MI - E[MI]}{\max(H(U), H(V)) - E[MI]}}
#'
#' where \eqn{MI} is the mutual information, \eqn{E[MI]} is the expected mutual
#' information under random permutation, and \eqn{H(U)} and \eqn{H(V)} are the
#' entropies of the two clusterings.
#'
#' @references
#' Vinh, N. X., Epps, J., and Bailey, J. (2010). Information theoretic measures
#' for clusterings comparison: Variants, properties, normalization and correction
#' for chance. Journal of Machine Learning Research, 11, 2837-2854.
#'
#' @export
ami_from_contingency <- function(cont_table) {
    # total number of isolates
    n <- sum(cont_table)

    if (n < 2) {
        return(NA_real_)
    }

    # row and column sums (unname to avoid propagating names)
    a <- unname(rowSums(cont_table))
    b <- unname(colSums(cont_table))

    # number of clusters in each assignment
    n_rows <- nrow(cont_table)
    n_cols <- ncol(cont_table)

    # compute entropy H(U) for clusters1
    p_a <- a / n
    h_u <- -sum(p_a[p_a > 0] * log(p_a[p_a > 0]))

    # compute entropy H(V) for clusters2
    p_b <- b / n
    h_v <- -sum(p_b[p_b > 0] * log(p_b[p_b > 0]))

    # compute mutual information MI(U, V)
    mi <- 0
    for (i in seq_len(n_rows)) {
        for (j in seq_len(n_cols)) {
            if (cont_table[i, j] > 0) {
                p_ij <- cont_table[i, j] / n
                mi <- mi + p_ij * log(p_ij / (p_a[i] * p_b[j]))
            }
        }
    }

    # compute expected mutual information E[MI] using the hypergeometric distribution
    emi <- 0
    for (i in seq_len(n_rows)) {
        for (j in seq_len(n_cols)) {
            # range of possible values for n_ij
            nij_min <- max(0, a[i] + b[j] - n)
            nij_max <- min(a[i], b[j])

            if (nij_min <= nij_max) {
                for (nij in nij_min:nij_max) {
                    if (nij > 0) {
                        # log probability of n_ij under hypergeometric distribution
                        log_prob <- lchoose(a[i], nij) +
                            lchoose(n - a[i], b[j] - nij) -
                            lchoose(n, b[j])

                        # contribution to expected MI
                        p_ij <- nij / n
                        p_ai <- a[i] / n
                        p_bj <- b[j] / n
                        log_ratio <- log(p_ij / (p_ai * p_bj))

                        emi <- emi + exp(log_prob) * p_ij * log_ratio
                    }
                }
            }
        }
    }

    # compute AMI
    max_entropy <- max(h_u, h_v)
    denominator <- max_entropy - emi

    if (abs(denominator) < .Machine$double.eps) {
        # perfect clustering or edge case
        if (abs(mi - emi) < .Machine$double.eps) {
            return(1)
        }
        return(0)
    }

    (mi - emi) / denominator
}

#' Adjusted Mutual Information from Cluster Assignments
#'
#' @description
#' Computes the Adjusted Mutual Information (AMI) between two cluster assignments.
#' This is a convenience wrapper around [ami_from_contingency()].
#'
#' @param clusters1 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#' @param clusters2 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#'
#' @return A numeric value representing the Adjusted Mutual Information.
#'
#' @seealso [ami_from_contingency()], [cluster_contingency_table()]
#'
#' @export
adjusted_mutual_information <- function(clusters1, clusters2) {
    cont_table <- cluster_contingency_table(clusters1, clusters2)
    ami_from_contingency(cont_table)
}

#' Fraction of Converts with Same Source Across Cluster Assignments
#'
#' @description
#' Compares two cluster assignments by examining patients categorized as "convert".
#' For each patient who is a convert in both assignments, checks whether the
#' source patient (categorized as "index" or "weak-index") in their respective
#' clusters is the same across both assignments.
#'
#' @param clusters1 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#' @param clusters2 A named vector of cluster assignments where names are isolate
#'   IDs and values are cluster numbers.
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param adm_seqs A vector of sequence IDs corresponding to admission positive patients.
#' @param dates A named vector mapping sequence IDs to dates.
#' @param surv_df A data frame with surveillance data containing columns:
#'   patient_id, genome_id, surv_date, result.
#'
#' @return A numeric value between 0 and 1 representing the fraction of convert
#'   patients whose source (index or weak-index) is the same in both cluster
#'   assignments. Returns NA if there are no common converts.
#'
#' @details
#' The function:
#' 1. Creates isolate lookups for both cluster assignments
#' 2. Categorizes patients in both using [cluster_patient_categorization()]
#' 3. Identifies patients categorized as "convert" in both assignments
#' 4. For each such convert, finds the "index" or "weak-index" patient in their
#'    cluster for each assignment
#' 5. Returns the fraction where the source patient matches
#'
#' @seealso [fraction_convert_same_source_from_lookups()], [cluster_patient_categorization()],
#'   [get_isolate_lookup()]
#'
#' @export
fraction_convert_same_source <- function(
    clusters1,
    clusters2,
    dna_aln,
    seq2pt,
    adm_seqs,
    dates,
    surv_df
) {
    # Create isolate lookups for both cluster assignments
    isolate_lookup1 <- get_isolate_lookup(
        clusters = clusters1,
        dna_aln = dna_aln,
        seq2pt = seq2pt,
        adm_seqs = adm_seqs,
        dates = dates,
        surv_df = surv_df
    )

    isolate_lookup2 <- get_isolate_lookup(
        clusters = clusters2,
        dna_aln = dna_aln,
        seq2pt = seq2pt,
        adm_seqs = adm_seqs,
        dates = dates,
        surv_df = surv_df
    )

    # Delegate to the lookup-based version
    fraction_convert_same_source_from_lookups(isolate_lookup1, isolate_lookup2, surv_df, surv_df)
}

#' Fraction of Converts with Same Source from Isolate Lookups
#'
#' @description
#' Compares two isolate lookups by examining patients categorized as "convert".
#' This is the internal implementation that works directly with isolate lookups.
#'
#' @param isolate_lookup1 A data frame from [get_isolate_lookup()] for the first
#'   cluster assignment.
#' @param isolate_lookup2 A data frame from [get_isolate_lookup()] for the second
#'   cluster assignment.
#' @param surv_df A data frame with surveillance data containing columns:
#'   patient_id, genome_id, surv_date, result.
#'
#' @return A numeric value between 0 and 1 representing the fraction of convert
#'   patients whose source (index or weak-index) is the same in both cluster
#'   assignments. Returns NA if there are no common converts.
#'
#' @seealso [fraction_convert_same_source()], [cluster_patient_categorization()]
#'
#' @export
fraction_convert_same_source_from_lookups <- function(
    isolate_lookup1,
    isolate_lookup2,
    surv_df_1,
    surv_df_2
) {
    # Helper: build convert -> sources mapping from isolate_lookup
    # Returns a named list: names are convert patient IDs, values are character
    # vectors of source patient IDs (one per cluster the convert appears in)
    get_convert_source_map <- function(isolate_lookup, surv_df) {
        categories <- cluster_patient_categorization(isolate_lookup, surv_df)
        non_single <- get_non_single_patient_clusters(isolate_lookup)
        categories <- categories[names(categories) %in% non_single]

        # Build convert -> sources list, accumulating across clusters
        convert_to_sources <- list()
        for (cl in names(categories)) {
            cats <- categories[[cl]]
            converts <- names(cats)[cats %in% c("convert", "secondary-convert")]

            if (length(converts) == 0) {
                next
            }

            # The source must be an index, weak-index, or multiply-colonized-index.
            # Skip clusters where all patients are converts or adm-pos.
            sources <- names(cats)[cats %in% c("index", "weak-index", "multiply-colonized-index")]
            if (length(sources) == 0) {
                next
            }

            source <- sources[1]
            converts <- setdiff(converts, source)
            if (length(converts) == 0) {
                next
            }

            for (conv in converts) {
                convert_to_sources[[conv]] <- c(convert_to_sources[[conv]], source)
            }
        }
        convert_to_sources
    }

    # Get convert -> sources mappings for both assignments
    map1 <- get_convert_source_map(isolate_lookup1, surv_df_1)
    map2 <- get_convert_source_map(isolate_lookup2, surv_df_2)

    # Find common converts
    common_converts <- intersect(names(map1), names(map2))

    if (length(common_converts) == 0) {
        return(NA_real_)
    }

    # Count where any source matches between the two assignments
    matches <- vapply(
        common_converts,
        function(conv) {
            length(intersect(map1[[conv]], map2[[conv]])) > 0
        },
        logical(1)
    )
    print(matches)
    print(length(common_converts))
    sum(matches) / length(common_converts)
}
