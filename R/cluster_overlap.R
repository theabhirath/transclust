#' Calculate isolate-isolate overlap
#'
#' @description
#' This function calculates the overlap between isolates in a given trace matrix.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [get_isolate_lookup()].
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

#' Calculate isolate-isolate sequential overlap
#'
#' @description
#' This function calculates the sequential overlap between isolates in a given trace matrix.
#' Sequential overlap counts the number of days where the recipient is at a location that the
#' donor had previously been at (on a strictly earlier day), capturing potential indirect
#' transmission via shared locations.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [get_isolate_lookup()].
#' @param trace_mat A matrix with rows representing days and columns representing patients.
#'
#' @return A data frame with columns: iso_donor, iso_recipient, overlap_days.
#'
#' @export
isolate_isolate_sequential_overlap <- function(isolate_lookup, trace_mat) {
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
            # if there is ANY spatiotemporal (concurrent) overlap, this pair is
            # invalid for sequential overlap — set to 0
            donor_locs <- pt_sub[1, ]
            recip_locs <- pt_sub[2, ]
            has_concurrent <- any(
                donor_locs == recip_locs & donor_locs > 0 & recip_locs > 0
            )
            seq_overlap <- 0L
            n_days <- length(donor_locs)
            if (!has_concurrent && n_days > 1L) {
                # count days where the recipient is at a location that the donor
                # had previously been at (on a strictly earlier day)
                donor_present <- donor_locs > 0
                unique_donor_locs <- unique(donor_locs[donor_present])
                for (loc in unique_donor_locs) {
                    first_donor_day <- which(donor_locs == loc)[1L]
                    recip_at_loc <- which(recip_locs == loc)
                    seq_overlap <- seq_overlap + sum(recip_at_loc > first_donor_day)
                }
            }
            # increment the counter and add the donor, recipient, and overlap to the vectors
            k <- k + 1L
            iso_donor_vec[k] <- iso_donor
            iso_rec_vec[k] <- iso_recipient
            overlap_days_vec[k] <- seq_overlap
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
#'                       [get_isolate_lookup()].
#' @param iso_overlap_df A data frame with overlap information for isolate pairs. For more information, see
#'                   [isolate_isolate_overlap()].
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

#' Calculate the fraction of converts with overlap
#'
#' @description
#' This function calculates the fraction of converts with overlap for each cluster.
#'
#' @param cluster_overlap_df A data frame with overlap information for isolate pairs. For more
#'                            information, see [cluster_isolate_overlap()].
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                        other relevant epidemiological information. For more information, see
#'                        [get_isolate_lookup()].
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
