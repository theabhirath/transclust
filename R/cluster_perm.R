#' Permutation test for cluster overlap fractions
#'
#' @description
#' Performs a permutation test to assess whether the observed fraction of converts
#' with overlap is significantly different from what would be expected by chance.
#' Tests overlap at facility, floor, and room levels.
#'
#' @param clusters A named numeric vector of cluster assignments.
#' @param dna_aln A DNA alignment object (used for isolate IDs via labels()).
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param adm_seqs A vector of sequence IDs which correspond to admission positive sequences.
#' @param adm_pos_pt_seqs A vector of sequence IDs for patients who are admission positive.
#' @param dates A named vector mapping sequence IDs to dates.
#' @param surv_df A data frame with surveillance data.
#' @param facility_trace A matrix with patient IDs as row names and dates as column names.
#' @param floor_trace A matrix with floor-level trace data (same structure as facility_trace).
#' @param room_trace A matrix with room-level trace data (same structure as facility_trace).
#' @param nperm Number of permutations to perform. Default 1000.
#' @param num_cores Number of cores for parallel processing. Default is detectCores() - 1.
#'
#' @returns A list containing:
#'   \itemize{
#'     \item \code{observed}: Named list with facility, floor, room fractions for observed data
#'     \item \code{permuted}: Array of dimensions (n_clusters, 3 trace_types, nperm)
#'     \item \code{valid_clusters}: Vector of cluster IDs that have more than one patient
#'   }
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @export
cluster_overlap_perm_test <- function(
    clusters,
    dna_aln,
    seq2pt,
    adm_seqs,
    adm_pos_pt_seqs,
    dates,
    surv_df,
    facility_trace,
    floor_trace,
    room_trace,
    nperm = 1000,
    num_cores = detectCores() - 1
) {
    cluster_names <- names(clusters)

    # Create observed isolate_lookup
    observed_lookup <- get_isolate_lookup(
        clusters = clusters,
        dna_aln = dna_aln,
        seq2pt = seq2pt,
        adm_seqs = adm_seqs,
        dates = dates,
        surv_df = surv_df
    )

    # Get non-single-patient clusters to analyze
    valid_clusters <- get_non_single_patient_clusters(observed_lookup)

    if (length(valid_clusters) == 0) {
        stop("No clusters with more than one patient found.")
    }

    # =========================================================================
    # OPTIMIZATION: Precompute isolate-isolate overlaps once
    # The overlap between isolates depends only on patient locations over time,
    # not on cluster assignments. So we compute this expensive O(n^2) operation
    # once and reuse it for all permutations.
    # =========================================================================
    iso_overlap_facility <- isolate_isolate_overlap(observed_lookup, facility_trace)
    iso_overlap_floor <- isolate_isolate_overlap(observed_lookup, floor_trace)
    iso_overlap_room <- isolate_isolate_overlap(observed_lookup, room_trace)

    # Precompute a base lookup without cluster dependency (for fast updates)
    base_lookup <- observed_lookup[, c("isolate_id", "patient_id", "date", "adm_pos", "prev_surv")]

    # Calculate observed overlap fractions using precomputed overlaps
    observed_fractions <- calculate_overlap_fractions_fast(
        observed_lookup,
        valid_clusters,
        iso_overlap_facility,
        iso_overlap_floor,
        iso_overlap_room
    )

    # Create eligibility matrices for permutation
    elig_mats <- create_eligibility_matrices(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs)

    # Calculate cluster counts for each patient category
    cluster_per_index_start <- cluster_per_patient(elig_mats$index_start)
    cluster_per_index_not_start <- cluster_per_patient(elig_mats$index_not_start)
    cluster_per_convert <- cluster_per_patient(elig_mats$convert)

    # Calculate patient counts per cluster for each category
    index_start_per_cluster <- pt_per_cluster(elig_mats$index_start, clusters)
    index_not_start_per_cluster <- pt_per_cluster(elig_mats$index_not_start, clusters)
    convert_per_cluster <- pt_per_cluster(elig_mats$convert, clusters)

    # Run permutations in parallel with progress bar
    perm_results <- pbmclapply(
        seq_len(nperm),
        function(n) {
            # Create permuted cluster assignment
            perm_clust <- assign_permuted_clusters(
                cluster_names,
                elig_mats,
                list(
                    index_start = cluster_per_index_start,
                    index_not_start = cluster_per_index_not_start,
                    convert = cluster_per_convert
                ),
                list(
                    index_start = index_start_per_cluster,
                    index_not_start = index_not_start_per_cluster,
                    convert = convert_per_cluster
                )
            )

            # Fast lookup update: just update the cluster column
            perm_lookup <- base_lookup
            perm_lookup$cluster <- perm_clust[perm_lookup$isolate_id]

            # Get valid clusters for this permutation
            perm_valid_clusters <- get_non_single_patient_clusters(perm_lookup)

            # Calculate overlap fractions using precomputed overlaps
            calculate_overlap_fractions_fast(
                perm_lookup,
                perm_valid_clusters,
                iso_overlap_facility,
                iso_overlap_floor,
                iso_overlap_room
            )
        },
        mc.cores = num_cores
    )

    # Organize results into array
    trace_types <- c("facility", "floor", "room")
    n_clusters <- length(valid_clusters)

    perm_array <- array(
        dim = c(n_clusters, 3, nperm),
        dimnames = list(
            as.character(valid_clusters),
            trace_types,
            seq_len(nperm)
        )
    )

    for (i in seq_len(nperm)) {
        for (j in seq_along(trace_types)) {
            trace <- trace_types[j]
            perm_frac <- perm_results[[i]][[trace]]
            # Match clusters - permuted may have different valid clusters
            matched <- perm_frac[as.character(valid_clusters)]
            perm_array[, j, i] <- matched
        }
    }

    list(
        observed = observed_fractions,
        permuted = perm_array,
        valid_clusters = valid_clusters
    )
}

# =============================================================================
# Helper functions for permutation tests
# =============================================================================

#' Calculate overlap fractions for all three trace types (original version)
#' @noRd
calculate_overlap_fractions <- function(
    isolate_lookup,
    valid_clusters,
    facility_trace,
    floor_trace,
    room_trace
) {
    # Filter lookup to valid clusters for overlap calculation
    lookup_filtered <- isolate_lookup[isolate_lookup$cluster %in% valid_clusters, ]

    # Calculate isolate-isolate overlap for each trace type
    iso_overlap_facility <- isolate_isolate_overlap(lookup_filtered, facility_trace)
    iso_overlap_floor <- isolate_isolate_overlap(lookup_filtered, floor_trace)
    iso_overlap_room <- isolate_isolate_overlap(lookup_filtered, room_trace)

    # Calculate cluster-isolate overlap
    cluster_overlap_facility <- cluster_isolate_overlap(lookup_filtered, iso_overlap_facility)
    cluster_overlap_floor <- cluster_isolate_overlap(lookup_filtered, iso_overlap_floor)
    cluster_overlap_room <- cluster_isolate_overlap(lookup_filtered, iso_overlap_room)

    # Calculate fraction of converts with overlap
    list(
        facility = fraction_converts_with_overlap(cluster_overlap_facility, lookup_filtered),
        floor = fraction_converts_with_overlap(cluster_overlap_floor, lookup_filtered),
        room = fraction_converts_with_overlap(cluster_overlap_room, lookup_filtered)
    )
}

#' Calculate overlap fractions using precomputed isolate-isolate overlaps (fast version)
#' @noRd
calculate_overlap_fractions_fast <- function(
    isolate_lookup,
    valid_clusters,
    iso_overlap_facility,
    iso_overlap_floor,
    iso_overlap_room
) {
    # Filter lookup to valid clusters for overlap calculation
    lookup_filtered <- isolate_lookup[isolate_lookup$cluster %in% valid_clusters, ]

    # Use precomputed isolate-isolate overlaps to calculate cluster-isolate overlap
    cluster_overlap_facility <- cluster_isolate_overlap(lookup_filtered, iso_overlap_facility)
    cluster_overlap_floor <- cluster_isolate_overlap(lookup_filtered, iso_overlap_floor)
    cluster_overlap_room <- cluster_isolate_overlap(lookup_filtered, iso_overlap_room)

    # Calculate fraction of converts with overlap
    list(
        facility = fraction_converts_with_overlap(cluster_overlap_facility, lookup_filtered),
        floor = fraction_converts_with_overlap(cluster_overlap_floor, lookup_filtered),
        room = fraction_converts_with_overlap(cluster_overlap_room, lookup_filtered)
    )
}

#' Create eligibility matrices for permutation
#' @noRd
create_eligibility_matrices <- function(clusters, seq2pt, adm_seqs, adm_pos_pt_seqs) {
    cluster_names <- names(clusters)

    # Identify index patients who started clusters vs those who didn't
    index_pt_start_seqs <- unlist(sapply(adm_seqs, function(seq_id) {
        cluster_names[clusters == clusters[seq_id] & seq2pt[cluster_names] == seq2pt[seq_id]]
    }))
    index_seqs_start <- intersect(index_pt_start_seqs, cluster_names)
    index_seqs_not_start <- intersect(setdiff(adm_pos_pt_seqs, index_seqs_start), cluster_names)
    convert_seqs <- setdiff(cluster_names, c(index_seqs_start, index_seqs_not_start))

    # Create eligibility matrix
    get_elig_mat <- function(seqs) {
        cbind(
            seq = as.character(seqs),
            patient = seq2pt[seqs],
            cluster = clusters[seqs],
            comb = paste(seq2pt[seqs], clusters[seqs], sep = "-")
        )
    }

    list(
        index_start = get_elig_mat(index_seqs_start),
        index_not_start = get_elig_mat(index_seqs_not_start),
        convert = get_elig_mat(convert_seqs)
    )
}

#' Count clusters per patient from eligibility matrix
#' @noRd
cluster_per_patient <- function(elig_mat) {
    if (nrow(elig_mat) == 0) {
        return(setNames(integer(0), character(0)))
    }

    pt_unique_sorted <- sort(unique(elig_mat[, "patient"]))
    setNames(
        vapply(
            pt_unique_sorted,
            function(pt) {
                length(unique(elig_mat[elig_mat[, "patient"] == pt, "cluster"]))
            },
            integer(1)
        ),
        pt_unique_sorted
    )
}

#' Count patients per cluster from eligibility matrix
#' @noRd
pt_per_cluster <- function(elig_mat, clusters) {
    clusters_unique_sorted <- sort(unique(clusters))
    setNames(
        vapply(
            clusters_unique_sorted,
            function(clust) {
                if (nrow(elig_mat) == 0) {
                    return(0L)
                }
                length(unique(elig_mat[elig_mat[, "cluster"] == clust, "patient"]))
            },
            integer(1)
        ),
        clusters_unique_sorted
    )
}

#' Assign patients to clusters randomly while preserving structure
#' @noRd
assign_pt_clusters <- function(elig_vec, cluster_per_pt, pt_per_cluster, rand_clusters) {
    if (length(cluster_per_pt) == 0) {
        return(rand_clusters)
    }

    names_pt_per_cluster <- names(pt_per_cluster)

    for (pt in names(sort(cluster_per_pt, decreasing = TRUE))) {
        pt_clust_ids <- unique(elig_vec[elig_vec[, "patient"] == pt, "comb"])
        cluster_assign <- sample(names_pt_per_cluster[pt_per_cluster > 0], length(pt_clust_ids))
        pt_per_cluster[cluster_assign] <- pt_per_cluster[cluster_assign] - 1

        for (pt_c_i in seq_along(pt_clust_ids)) {
            seqs_to_assign <- elig_vec[elig_vec[, "comb"] == pt_clust_ids[pt_c_i], "seq"]
            rand_clusters[seqs_to_assign] <- as.numeric(cluster_assign[pt_c_i])
        }
    }
    rand_clusters
}

#' Assign permuted clusters
#' @noRd
assign_permuted_clusters <- function(cluster_names, elig_mats, cluster_per_pt, pt_per_clust) {
    perm_clust <- setNames(rep(-1, length(cluster_names)), cluster_names)

    # Assign converts first
    perm_clust <- assign_pt_clusters(
        elig_mats$convert,
        cluster_per_pt$convert,
        pt_per_clust$convert,
        perm_clust
    )

    # Then index patients who started clusters
    perm_clust <- assign_pt_clusters(
        elig_mats$index_start,
        cluster_per_pt$index_start,
        pt_per_clust$index_start,
        perm_clust
    )

    # Then index patients who didn't start clusters
    perm_clust <- assign_pt_clusters(
        elig_mats$index_not_start,
        cluster_per_pt$index_not_start,
        pt_per_clust$index_not_start,
        perm_clust
    )

    perm_clust
}
