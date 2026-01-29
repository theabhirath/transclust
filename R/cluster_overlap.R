#' Calculate isolate-isolate overlap
#'
#' @description
#' This function calculates the overlap between isolates in a given trace matrix.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [`get_isolate_lookup`].
#' @param trace_mat A matrix with rows representing days and columns representing patients.
#'
#' @return A data frame with columns: iso_donor, iso_recipient, overlap_days.
#'
#' @export
isolate_isolate_overlap <- function(isolate_lookup, trace_mat) {
    # only keep isolates that:
    # - have non-missing dates and prev_surv
    # - have patients present in the trace matrix
    valid <- !is.na(isolate_lookup$date) &
        !is.na(isolate_lookup$prev_surv) &
        isolate_lookup$patient_id %in% rownames(trace_mat)
    iso <- isolate_lookup[valid, ]
    num_isolates <- nrow(iso)

    # precompute the patient row indices in trace_mat (avoid repeated character lookups)
    pt_row_idx <- match(as.character(iso$patient_id), rownames(trace_mat))

    # preallocate the maximum possible size (num_isolates^2), then trim at the end
    max_pairs <- num_isolates * num_isolates
    iso_donor_vec <- character(max_pairs)
    iso_rec_vec <- character(max_pairs)
    overlap_days_vec <- integer(max_pairs)
    k <- 0L

    # iterate over all pairs of isolates
    for (i in seq_len(num_isolates)) {
        pt_i <- iso$patient_id[i]
        last_surv <- iso$prev_surv[i] # donor window start
        row_i <- pt_row_idx[i]
        iso_donor <- iso$isolate_id[i]

        for (j in seq_len(num_isolates)) {
            iso_recipient <- iso$isolate_id[j]
            pt_j <- iso$patient_id[j]
            recipient_date <- iso$date[j] # recipient collection date
            row_j <- pt_row_idx[j]
            # skip if:
            # - the donor and recipient are the same patient
            # - the donor window starts after the recipient collection date
            if (pt_i == pt_j || last_surv > recipient_date) {
                next
            }
            # get the time range between the donor and recipient
            cols <- as.integer(last_surv):as.integer(recipient_date)
            # get the trace matrix subset for the donor and recipient
            pt_sub <- trace_mat[c(row_i, row_j), cols, drop = FALSE]
            # count the number of days where the donor and recipient are both present
            overlap_days <- pt_sub[1, ] == pt_sub[2, ] & pt_sub[1, ] > 0 & pt_sub[2, ] > 0
            # increment the counter and add the donor, recipient, and overlap to the vectors
            k <- k + 1L
            iso_donor_vec[k] <- iso_donor
            iso_rec_vec[k] <- iso_recipient
            overlap_days_vec[k] <- sum(overlap_days)
        }
    }

    # if no overlaps were found, return an empty data frame
    # otherwise, return the overlaps found as a data frame
    overlap_df <- if (k == 0L) {
        data.frame(
            iso_donor = character(0),
            iso_recipient = character(0),
            overlap_days = integer(0)
        )
    } else {
        data.frame(
            iso_donor = iso_donor_vec[seq_len(k)],
            iso_recipient = iso_rec_vec[seq_len(k)],
            overlap_days = overlap_days_vec[seq_len(k)]
        )
    }

    overlap_df
}

#' Calculate isolate-isolate overlap for clusters
#'
#' @description
#' This function considers every isolate in a cluster as a recipient, and returns if there is
#' any overlap when other isolates in the cluster are the donors.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [`get_isolate_lookup`].
#' @param overlap_df A data frame with overlap information for isolate pairs. For more information, see
#'                   [`isolate_isolate_overlap`].
#'
#' @return A data frame with columns: cluster, isolate_id, overlap.
#'
#' @export
cluster_isolate_overlap <- function(isolate_lookup, iso_overlap_df) {
    # get non-single patient clusters
    clusters <- get_non_single_patient_clusters(isolate_lookup)

    # list of data frames for each cluster
    out_list <- lapply(clusters, function(cl) {
        idx_cl <- isolate_lookup$cluster == cl
        isolates <- isolate_lookup$isolate_id[idx_cl]
        adm_pos <- isolate_lookup$adm_pos[idx_cl]

        # overlaps within this cluster only
        od <- iso_overlap_df[
            iso_overlap_df$iso_donor %in% isolates & iso_overlap_df$iso_recipient %in% isolates,
            ,
            drop = FALSE
        ]

        ## TODO: check if there is a way to refine this overlap to be more accurate
        # currently all this cares about is if there is any overlap, not the actual overlap days
        od_overlap_days <- od$overlap_days
        any_overlap_by_rec <- if (nrow(od) > 0 && max(od_overlap_days) > 0) {
            # this could be min if we want to look at bare minimum overlap
            tapply(od_overlap_days[od_overlap_days > 0], od$iso_recipient[od_overlap_days > 0], any)
        } else {
            FALSE
        }

        # align to isolates; missing => FALSE; admission-positive => NA
        overlaps <- any_overlap_by_rec[isolates]
        overlaps[is.na(overlaps)] <- FALSE
        overlaps[adm_pos] <- NA

        data.frame(
            cluster = rep(cl, length(isolates)),
            isolate_id = isolates,
            overlap = overlaps,
            row.names = NULL
        )
    })

    # combine the data frames for each cluster into a single data frame
    do.call(rbind, out_list)
}

#' Categorize cluster overlap
#'
#' @description
#' This function categorizes the overlap of isolates in a cluster, based on the earliest collected isolate
#' from all patients deemed to be in the cluster (hereafter referred to as the "index isolate"), and overlap
#' explanations for other isolates in the cluster. The categories are:
#'  - "patient-to-patient": if the index isolate is admission-positive and there is overlap explanation for
#'    all other converts in the cluster.
#'  - "weak-index": if the index isolate is not admission-positive but is the first surveillance for the patient
#'    after admission and there is overlap explanation for all other isolates in the cluster.
#'  - "missing-intermediate": if the index isolate is admission-positive but at least one other convert
#'    in the cluster has no overlap explanation.
#'  - "false-negative-index": if the index isolate is not admission-positive but there is overlap explanation for
#'    all other converts in the cluster barring one (deemed to be the false negative index).
#'  - "missing-source": if the index isolate is not admission-positive and there is no overlap explanation for
#'    more than one convert in the cluster.
#'  - "multiply-colonized-index": if the index isolate is not in the cluster but is admission-positive, this is
#'    a "multiply-colonized index" if there is overlap explanation for all other isolates in the cluster.
#'  - "multiply-colonized-index-missing-intermediate": if the index isolate is not in the cluster but is
#'    admission-positive, and at least one other convert in the cluster has no overlap explanation.
#'  - "inexplicable": catch-all category for cases that do not fit into the other categories. Currently,
#'    index and overlap explanations are not enough to explain these clusters.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [`get_isolate_lookup`].
#' @param cluster_overlap_df A data frame with overlap information for isolate pairs. For more information, see
#'                   [`cluster_isolate_overlap`].
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
            # if overlap_expl_sub is all NA, the cluster has only admission positive isolates
            if (all(is.na(overlap_expl_sub))) {
                return("all-admission-positive")
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

#' Calculate the fraction of converts with overlap
#'
#' @description
#' This function calculates the fraction of converts with overlap for each cluster.
#'
#' @param cluster_overlap_df A data frame with overlap information for isolate pairs. For more
#'                            information, see [`cluster_isolate_overlap`].
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                        other relevant epidemiological information. For more information, see
#'                        [`get_isolate_lookup`].
#'
#' @return A named vector with the fraction of converts with overlap for each cluster.
#'
#' @export
fraction_converts_with_overlap <- function(cluster_overlap_df, isolate_lookup) {
    setNames(
        vapply(
            unique(cluster_overlap_df$cluster),
            function(cl) {
                lookup_sub <- isolate_lookup[isolate_lookup$cluster == cl, ]
                convert_isolates <- lookup_sub$isolate_id[
                    lookup_sub$adm_pos == FALSE
                ]
                num_converts <- length(unique(lookup_sub$isolate_id[!lookup_sub$adm_pos]))
                num_converts_with_overlap <- sum(
                    cluster_overlap_df$overlap[
                        cluster_overlap_df$isolate_id %in% convert_isolates
                    ],
                    na.rm = TRUE
                )
                if (num_converts > 0) {
                    num_converts_with_overlap / num_converts
                } else {
                    NA
                }
            },
            numeric(1)
        ),
        unique(cluster_overlap_df$cluster)
    )
}
