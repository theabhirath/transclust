#' Categorize cluster overlap
#'
#' @description
#' This function categorizes the overlap of isolates in a cluster, based on the earliest collected isolate
#' from all patients deemed to be in the cluster (hereafter referred to as the "index isolate"), and overlap
#' explanations for other isolates in the cluster. The categories are:
#' \itemize{
#'   \item "patient-to-patient": if the index isolate is admission-positive and there is overlap explanation for
#'     all other converts in the cluster.
#'   \item "weak-index": if the index isolate is not admission-positive but is the first surveillance for the patient
#'     after admission and there is overlap explanation for all other isolates in the cluster.
#'   \item "missing-intermediate": if the index isolate is admission-positive but at least one other convert
#'     in the cluster has no overlap explanation.
#'   \item "false-negative-index": if the index isolate is not admission-positive but there is overlap explanation for
#'     all other converts in the cluster barring one (deemed to be the false negative index).
#'   \item "missing-source": if the index isolate is not admission-positive and there is no overlap explanation for
#'     more than one convert in the cluster.
#'   \item "multiply-colonized-index": if the index isolate is not in the cluster but is admission-positive, this is
#'     a "multiply-colonized index" if there is overlap explanation for all other isolates in the cluster.
#'   \item "multiply-colonized-index-missing-intermediate": if the index isolate is not in the cluster but is
#'     admission-positive, and at least one other convert in the cluster has no overlap explanation.
#'   \item "inexplicable": catch-all category for cases that do not fit into the other categories. Currently,
#'     index and overlap explanations are not enough to explain these clusters.
#' }
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [get_isolate_lookup()].
#' @param cluster_overlap_df A data frame with overlap information for isolate pairs. For more information, see
#'                   [cluster_isolate_overlap()].
#' @param surv_df A data frame with surveillance information for isolates.
#'
#' @return A vector with the categorization of the overlap for each cluster.
#'
#' @export
categorize_cluster_overlap <- function(isolate_lookup, cluster_overlap_df, surv_df) {
    unique_clusters <- as.character(unique(cluster_overlap_df$cluster))
    # iterate over each cluster in cluster_overlap_df
    setNames(
        sapply(unique_clusters, function(cl) {
            # overlap for cluster
            cluster_overlap_sub <- cluster_overlap_df[cluster_overlap_df$cluster == cl, ]

            # find the index isolate (earliest date, with tiebreakers)
            cluster_patients <- unique(isolate_lookup$patient_id[isolate_lookup$cluster == cl])
            lookup_sub <- isolate_lookup[isolate_lookup$patient_id %in% cluster_patients, ]
            index_lookup <- find_index_isolate(cl, cluster_patients, isolate_lookup)

            # overlap explanation for all isolates except the index
            # we ignore other admission-positive isolates
            overlap_expl_sub <- cluster_overlap_sub$overlap[
                cluster_overlap_sub$isolate_id != index_lookup$isolate_id
            ]
            # if overlap_expl_sub is all NA
            if (all(is.na(overlap_expl_sub))) {
                # check if all isolates are admission positive
                if (all(lookup_sub$adm_pos)) {
                    return("all-admission-positive")
                } else {
                    return("inexplicable")
                }
            }
            cluster_overlap_expl <- all(overlap_expl_sub, na.rm = TRUE)
            # check if the index isolate is in the cluster
            cluster_isolates <- lookup_sub$isolate_id[lookup_sub$cluster == cl]
            index_in_cluster <- index_lookup$isolate_id %in% cluster_isolates
            if (index_in_cluster) {
                # if the index isolate is admission-positive, it is a "true index"
                if (index_lookup$adm_pos) {
                    # for all other isolates, overlap should be TRUE or NA, not FALSE
                    # if this is the case, the cluster category is "patient-to-patient"
                    if (cluster_overlap_expl) {
                        "patient-to-patient"
                    } else {
                        "missing-intermediate"
                    }
                } else {
                    # if the isolate is not admission positive but is the first surveillance for the patient
                    # after admission, this is a "weak index"
                    surv_index_pt <- surv_df[surv_df$patient_id == index_lookup$patient_id, ]
                    # this is a "weak index" if the earliest surveillance date is the same as the isolate date
                    # and we have overlap explanation for all other isolates
                    if (min(surv_index_pt$surv_date) == index_lookup$date && cluster_overlap_expl) {
                        "weak-index"
                    } else {
                        # if we have overlap explanations for all isolates except one convert,
                        # this is a "false negative index"
                        # sum should be equal to the number of converts in the cluster minus one
                        converts_lookup_sub <- lookup_sub[lookup_sub$adm_pos == FALSE, ]
                        num_overlap_expl <- sum(cluster_overlap_sub$overlap[
                            cluster_overlap_sub$isolate_id %in%
                                cluster_isolates &
                                cluster_overlap_sub$isolate_id %in% converts_lookup_sub$isolate_id
                        ])
                        num_converts <- length(converts_lookup_sub$isolate_id)
                        if (num_overlap_expl == num_converts) {
                            "false-negative-index"
                        } else {
                            "missing-source"
                        }
                    }
                }
            } else {
                # if the earliest isolate is not in the cluster but is admission positive
                # this is a "multiply-colonized index"
                if (index_lookup$adm_pos) {
                    if (cluster_overlap_expl) {
                        "multiply-colonized-index"
                    } else {
                        "multiply-colonized-index-missing-intermediate"
                    }
                } else {
                    "inexplicable"
                }
            }
        }),
        unique_clusters
    )
}

#' Categorize patients within clusters by epidemiological status
#'
#' @description
#' Categorizes each patient in a cluster based on their role in transmission:
#' \itemize{
#'   \item index: admission positive, earliest isolate, isolate in cluster
#'   \item multiply-colonized-index: admission positive, earliest isolate, isolate not in cluster
#'   \item weak-index: not admission positive, but first surveillance is positive and in cluster
#'   \item convert: had surveillance before first positive
#'   \item adm-pos: first positive is in cluster and is first surveillance
#'   \item adm-pos-convert: first positive is not in cluster but is first surveillance
#'   \item secondary-convert: first positive is not in cluster and had prior surveillance
#' }
#'
#' @param isolate_lookup A lookup table from [get_isolate_lookup()].
#' @param surv_df A data frame with surveillance data (patient_id, genome_id, surv_date, result).
#'
#' @returns A named list where each element is a named character vector mapping
#'   patient_id to their category within that cluster.
#'
#' @export
cluster_patient_categorization <- function(isolate_lookup, surv_df) {
    setNames(
        lapply(unique(isolate_lookup$cluster), function(cl) {
            categorize_cluster_patients(cl, isolate_lookup, surv_df)
        }),
        as.character(unique(isolate_lookup$cluster))
    )
}

#' Categorize patients within a single cluster
#' @noRd
categorize_cluster_patients <- function(cluster_id, isolate_lookup, surv_df) {
    cluster_patients <- unique(
        isolate_lookup$patient_id[isolate_lookup$cluster == cluster_id]
    )
    results <- setNames(character(length(cluster_patients)), cluster_patients)

    # Find and categorize the index patient (earliest isolate)
    index_row <- find_index_isolate(cluster_id, cluster_patients, isolate_lookup)
    index_patient <- as.character(index_row$patient_id)
    results[index_patient] <- categorize_index_patient(index_row, cluster_id)

    # Categorize remaining patients
    other_patients <- setdiff(cluster_patients, index_row$patient_id)
    for (patient in other_patients) {
        results[as.character(patient)] <- categorize_non_index_patient(
            patient,
            cluster_id,
            isolate_lookup,
            surv_df
        )
    }

    # If all non-index patients are "adm-pos", reclassify index
    if (
        length(other_patients) > 0 &&
            all(results[as.character(other_patients)] == "adm-pos")
    ) {
        # check if the index isolate is admission positive
        results[index_patient] <- if (index_row$adm_pos) {
            "adm-pos"
        } else {
            "convert"
        }
    }

    results
}

#' Find the index isolate for a cluster (earliest date, with tiebreakers)
#' @noRd
find_index_isolate <- function(cluster_id, cluster_patients, isolate_lookup) {
    lookup_sub <- isolate_lookup[isolate_lookup$patient_id %in% cluster_patients, ]
    candidates <- lookup_sub[lookup_sub$date == min(lookup_sub$date), ]

    # Tiebreakers: prefer in-cluster, then admission positive, then lowest isolate_id
    if (nrow(candidates) > 1 && any(candidates$cluster == cluster_id)) {
        candidates <- candidates[candidates$cluster == cluster_id, ]
    }
    if (nrow(candidates) > 1 && any(candidates$adm_pos)) {
        candidates <- candidates[candidates$adm_pos, ]
    }
    if (nrow(candidates) > 1) {
        candidates <- candidates[which.min(candidates$isolate_id), ]
    }

    candidates[1, ]
}

#' Categorize the index patient based on admission status and cluster membership
#' @noRd
categorize_index_patient <- function(index_row, cluster_id) {
    is_in_cluster <- index_row$cluster == cluster_id
    is_first_surv <- index_row$prev_surv == index_row$date

    if (index_row$adm_pos) {
        if (is_in_cluster) "index" else "multiply-colonized-index"
    } else if (is_first_surv && is_in_cluster) {
        "weak-index"
    } else {
        "convert"
    }
}

#' Categorize a non-index patient based on their first positive surveillance
#' @noRd
categorize_non_index_patient <- function(patient, cluster_id, isolate_lookup, surv_df) {
    # Get patient's positive surveillances
    surv_patient <- surv_df[surv_df$patient_id == patient, ]
    pos_survs <- surv_patient[surv_patient$result == 1, ]

    # exceptional case, no positive surveillances - should not happen
    if (nrow(pos_survs) == 0) {
        return("no-pos-surv")
    }

    # Filter to positive surveillances with minimum date
    min_date <- min(pos_survs$surv_date)
    first_pos <- pos_survs[pos_survs$surv_date == min_date, ]

    # If multiple on same date, prefer one in cluster
    if (nrow(first_pos) > 1) {
        cluster_genome_ids <- isolate_lookup$isolate_id[isolate_lookup$cluster == cluster_id]
        in_cluster <- first_pos[first_pos$genome_id %in% cluster_genome_ids, ]
        if (nrow(in_cluster) > 0) {
            first_pos <- in_cluster
        }
    }

    # If still multiple, use lowest genome_id
    if (nrow(first_pos) > 1) {
        first_pos <- first_pos[which.min(first_pos$genome_id), ][1, ]
    }

    # Look up the isolate from that surveillance
    isolate_row <- isolate_lookup[isolate_lookup$isolate_id == first_pos$genome_id, ]

    if (nrow(isolate_row) == 0) {
        return("secondary-convert")
    }

    is_in_cluster <- isolate_row$cluster == cluster_id
    is_first_surv <- isolate_row$date == isolate_row$prev_surv

    if (is_in_cluster) {
        if (is_first_surv) "adm-pos" else "convert"
    } else {
        if (is_first_surv) "adm-pos-convert" else "secondary-convert"
    }
}
