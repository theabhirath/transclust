#' Analyze Intra-Cluster Genetic Variation
#'
#' @description
#' This function analyzes genetic diversity within clusters by identifying variant sites and calculating
#' mutational spectra (e.g. transition and transversion rates). It uses a DNA alignment and a vector
#' of cluster assignments (named by sequence IDs) and returns a data frame summarizing the diversity
#' measures for each cluster as well as for the overall population.
#'
#' @param clusters A named vector of cluster assignments (cluster IDs) with names corresponding to sequence IDs.
#' @param dna_aln A matrix representing the DNA sequence alignment, where rows are sequences (named by sequence IDs)
#'                and columns represent nucleotide positions.
#' @param var_pos A logical vector indicating which positions in the alignment are variable.
#'
#' @return A data frame with one row per cluster (plus a row for the overall population) containing the
#'   number of variable sites and the rates for six mutation types. Mutation rates are rounded to three decimals.
#'
#' @export
intra_cluster_genetic_var_analysis <- function(clusters, dna_aln, var_pos) {
    # Pre-compute a character matrix of the alignment
    dna_aln_char <- matrix(as.character(dna_aln), nrow = nrow(dna_aln), dimnames = dimnames(dna_aln))

    # Identify cluster IDs (excluding a potential singleton "1" if present)
    cluster_ids <- setdiff(sort(unique(clusters)), 1)

    # Compute major alleles across the entire alignment
    major_alleles_all <- apply(dna_aln_char, 2, function(pos_aln) {
        allele_tab <- table(pos_aln)
        names(allele_tab)[which.max(allele_tab)]
    })

    # Filter the alignment to only include variable positions
    dna_aln_char <- dna_aln_char[, var_pos, drop = FALSE]
    major_alleles_all <- major_alleles_all[var_pos]

    # Pre-allocate vectors for results
    n_iter <- length(cluster_ids) + 1 # +1 for the overall population
    num_var_sites <- numeric(n_iter)
    gc_at_transition_rate <- numeric(n_iter)
    at_gc_transition_rate <- numeric(n_iter)
    gc_ta_transversion_rate <- numeric(n_iter)
    at_cg_transversion_rate <- numeric(n_iter)
    at_ta_transversion_rate <- numeric(n_iter)
    gc_cg_transversion_rate <- numeric(n_iter)

    # Define row names for the output: overall population and each cluster
    cluster_labels <- c("Pop. freq.", cluster_ids)

    iter <- 1
    # Loop over overall population (represented by -1) and each cluster
    for (cluster in c(-1, cluster_ids)) {
        if (cluster == -1) {
            seq_ids <- names(clusters)
        } else {
            seq_ids <- names(clusters)[clusters == cluster]
        }

        # Subset the alignment for the current group
        aln_subset <- dna_aln_char[seq_ids, , drop = FALSE]
        # Determine variable positions within the current group
        var_pos_cluster <- apply(aln_subset, 2, function(pos_aln) {
            sum(pos_aln == pos_aln[1]) < length(seq_ids)
        })
        var_idx <- which(var_pos_cluster)

        # Compute minor allele for each variable position
        minor_alleles <- vapply(var_idx, function(j) {
            pos_aln_full <- dna_aln_char[, j]
            allele_tab <- table(pos_aln_full)
            major_allele <- major_alleles_all[j]
            cluster_alleles <- unique(aln_subset[, j])
            allele_count <- length(cluster_alleles)
            # Check if the major allele is present in the cluster
            is_major_present <- major_allele %in% cluster_alleles
            if (allele_count == 1) {
                # Cluster has a major or minor allele
                if (is_major_present) "X" else "Y"
            } else if (allele_count == 2 && is_major_present) {
                # Cluster has a variable minor allele
                setdiff(cluster_alleles, major_allele)
            } else {
                "Z" # Cluster has multiple minor alleles
            }
        }, character(1))
        # Compute major allele for each variable position
        major_alleles <- major_alleles_all[var_idx]

        # Filter positions to those with valid nucleotide calls (both alleles among a, c, t, g)
        valid_pos <- which(
            minor_alleles != major_alleles &
                major_alleles %in% c("a", "c", "t", "g") &
                minor_alleles %in% c("a", "c", "t", "g")
        )

        num_sites <- length(valid_pos)
        num_var_sites[iter] <- num_sites

        # If no positions are variable or no valid sites are found, assign NA to mutation rates
        if (num_sites == 0 || !any(var_pos_cluster)) {
            gc_at_transition_rate[iter] <- NA
            at_gc_transition_rate[iter] <- NA
            gc_ta_transversion_rate[iter] <- NA
            at_cg_transversion_rate[iter] <- NA
            at_ta_transversion_rate[iter] <- NA
            gc_cg_transversion_rate[iter] <- NA
            iter <- iter + 1
            next
        }

        major_valid <- major_alleles[valid_pos]
        minor_valid <- minor_alleles[valid_pos]

        # Calculate mutation rates
        gc_at_transition_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "a") |
                (major_valid == "c" & minor_valid == "t")
        ) / num_sites
        at_gc_transition_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "g") |
                (major_valid == "t" & minor_valid == "c")
        ) / num_sites
        gc_ta_transversion_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "t") |
                (major_valid == "c" & minor_valid == "a")
        ) / num_sites
        at_cg_transversion_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "c") |
                (major_valid == "t" & minor_valid == "g")
        ) / num_sites
        at_ta_transversion_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "t") |
                (major_valid == "t" & minor_valid == "a")
        ) / num_sites
        gc_cg_transversion_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "c") |
                (major_valid == "c" & minor_valid == "g")
        ) / num_sites

        iter <- iter + 1
    }

    data.frame(
        num_var_sites = num_var_sites,
        GC_AT_transition_rate = round(gc_at_transition_rate, 3),
        AT_GC_transition_rate = round(at_gc_transition_rate, 3),
        GC_TA_transversion_rate = round(gc_ta_transversion_rate, 3),
        AT_TA_transversion_rate = round(at_ta_transversion_rate, 3),
        GC_CG_transversion_rate = round(gc_cg_transversion_rate, 3),
        AT_CG_transversion_rate = round(at_cg_transversion_rate, 3),
        row.names = cluster_labels,
        check.names = FALSE
    )
}

#' Compute various cluster-level properties
#'
#' @description
#' Calculates various epidemiological and genetic characteristics for a given
#' set of sequence IDs in a single cluster. This includes counts of patients,
#' convert vs. index statuses, genetic distances, cluster duration, and
#' timing-based metrics such as time to first/last acquisition.
#'
#' @param cluster_seqs A character vector of sequence IDs belonging to the cluster.
#' @param pt_trace A data frame or matrix of patient-level trace data (rows are days,
#'   columns are patient IDs).
#' @param seq2pt A named character vector mapping sequence IDs to patient IDs
#'   (names must be sequence IDs, values must be patient IDs).
#' @param ip_pt_seqs A character vector of sequence IDs corresponding to patients
#'   who tested positive on admission (\"intake positive\").
#' @param ip_seqs A character vector of sequence IDs presumed to be imported
#'   from outside (a subset of \code{ip_pt_seqs} or similar).
#' @param snp_dist A matrix of SNP distances between isolates (row and column
#'   names should match sequence IDs).
#' @param dates A named numeric or Date vector of isolate dates by sequence ID. Often
#'   these are numeric if read from Excel-style dates, in which case an origin
#'   should be considered (e.g. \code{\"1899-12-30\"}).
#' @param floor_trace (Optional) A data frame or matrix for floor-level tracing.
#' @param room_trace (Optional) A data frame or matrix for room-level tracing.
#'
#' @return A named numeric vector containing various cluster-level properties:
#'   \itemize{
#'     \item Number_of_patients
#'     \item Number_of_start_indexes
#'     \item Number_of_non-start_indexes
#'     \item Number_of_study_start_indexes
#'     \item Number_of_converts
#'     \item Number_of_initial_converts
#'     \item Number_of_later_converts
#'     \item Cluster_duration
#'     \item Mean_genetic_distance
#'     \item Median_genetic_distance
#'     \item Max_genetic_distance
#'     \item Number_of_converts_after_index
#'     \item Number_of_initial_converts_after_index
#'     \item Number_of_converts_after_index_seq
#'     \item Number_of_initial_converts_after_index_seq
#'     \item Number_of_converts_with_source
#'     \item Number_of_initial_converts_with_source
#'     \item Number_of_converts_with_floor_source
#'     \item Number_of_initial_converts_with_floor_source
#'     \item Number_of_converts_with_room_source
#'     \item Number_of_initial_converts_with_room_source
#'     \item Time_to_first_acquisition
#'     \item Time_to_last_acquisition
#'     \item Median_time_to_acquisition
#'     \item Time_from_index_to_first_convert
#'     \item Time_from_index_to_last_convert
#'   }
#' @keywords internal
cluster_properties <- function(cluster_seqs, pt_trace, seq2pt, ip_pt_seqs, ip_seqs, snp_dist,
                               dates, floor_trace = NULL, room_trace = NULL) {
    # Define names for all output properties
    attr_names <- c(
        "Number_of_patients",
        "Number_of_start_indexes",
        "Number_of_non-start_indexes",
        "Number_of_study_start_indexes",
        "Number_of_converts",
        "Number_of_initial_converts",
        "Number_of_later_converts",
        "Cluster_duration",
        "Mean_genetic_distance",
        "Median_genetic_distance",
        "Max_genetic_distance",
        "Number_of_converts_after_index",
        "Number_of_initial_converts_after_index",
        "Number_of_converts_after_index_seq",
        "Number_of_initial_converts_after_index_seq",
        "Number_of_converts_with_source",
        "Number_of_initial_converts_with_source",
        "Number_of_converts_with_floor_source",
        "Number_of_initial_converts_with_floor_source",
        "Number_of_converts_with_room_source",
        "Number_of_initial_converts_with_room_source",
        "Time_to_first_acquisition",
        "Time_to_last_acquisition",
        "Median_time_to_acquisition",
        "Time_from_index_to_first_convert",
        "Time_from_index_to_last_convert"
    )

    # Initialize output vector for cluster properties
    cluster_prop <- setNames(numeric(length(attr_names)), attr_names)

    # 1. Basic cluster-level counts #####
    # Number of unique patients in this cluster
    cluster_prop["Number_of_patients"] <- length(unique(seq2pt[cluster_seqs]))

    # Number of start-index patients (imported cases) who belong to this cluster
    # i.e. those with sequence IDs in ip_seqs
    cluster_prop["Number_of_start_indexes"] <- length(unique(seq2pt[cluster_seqs[cluster_seqs %in% ip_seqs]]))

    # Number of non-start index patients (any intake-positive but not the 'imported' subset)
    cluster_prop["Number_of_non-start_indexes"] <- length(setdiff(
        unique(seq2pt[cluster_seqs[cluster_seqs %in% ip_pt_seqs]]),
        unique(seq2pt[cluster_seqs[cluster_seqs %in% ip_seqs]])
    ))

    # Number of study start indexes
    cluster_prop["Number_of_study_start_indexes"] <- length(setdiff(
        unique(seq2pt[cluster_seqs[cluster_seqs %in% ip_seqs]]),
        unique(seq2pt[cluster_seqs[
            as.Date(dates[cluster_seqs], origin = "1899-12-30") >
                (as.Date(row.names(pt_trace)[1]) + 3)
        ]])
    ))

    # 2. Convert (non-index) patient counts #####
    # Number of converts: patients who are not in ip_pt_seqs
    # i.e. those who became positive while in the facility (not on admission).
    cluster_prop["Number_of_converts"] <- length(unique(
        seq2pt[cluster_seqs[!(cluster_seqs %in% ip_pt_seqs)]]
    ))

    # 3. Genetic distances within the cluster #####
    # Subset the SNP distance matrix to only the cluster sequences
    snp_dist_subset <- snp_dist[cluster_seqs, cluster_seqs, drop = FALSE]

    # Get the upper triangle of the SNP distance matrix (excluding the diagonal)
    seq_vals <- snp_dist_subset[upper.tri(snp_dist_subset)]

    cluster_prop["Mean_genetic_distance"] <- mean(seq_vals)
    cluster_prop["Median_genetic_distance"] <- median(seq_vals)
    cluster_prop["Max_genetic_distance"] <- max(seq_vals)

    # 4. Cluster duration metrics #####
    # Calculate earliest positive date for each patient in the cluster
    # Then compute the cluster duration = difference between min and max
    # of those earliest positive dates.
    earliest_pos_by_pt <- sapply(unique(seq2pt[cluster_seqs]), function(pt_id) {
        min(dates[cluster_seqs[seq2pt[cluster_seqs] == pt_id]])
    })
    earliest_pos_by_pt <- sort(earliest_pos_by_pt)

    cluster_prop["Cluster_duration"] <- max(earliest_pos_by_pt) - min(earliest_pos_by_pt)

    # Time to first acquisition = difference between the earliest two patients' dates
    cluster_prop["Time_to_first_acquisition"] <- earliest_pos_by_pt[2] - earliest_pos_by_pt[1]

    # Time to last acquisition = difference between the earliest and latest
    # in the sorted vector of earliest_pos_by_pt
    cluster_prop["Time_to_last_acquisition"] <- earliest_pos_by_pt[length(earliest_pos_by_pt)] - earliest_pos_by_pt[1]

    # Median time to acquisition (the median difference from the cluster's
    # earliest overall positive to each subsequent patient's earliest positive).
    cluster_prop["Median_time_to_acquisition"] <- median(sapply(
        earliest_pos_by_pt[2:length(earliest_pos_by_pt)], function(x) x - earliest_pos_by_pt[1]
    ))

    # 5. Detailed convert-patient calculations #####
    convert_pos <- numeric(0)
    conversion_day <- character(0)
    initial_convert_pts <- character(0)
    convert_pts <- character(0)

    # If we have at least one convert in the cluster:
    if (cluster_prop["Number_of_converts"] > 0) {
        # For each convert patient, find the earliest date in 'dates' for the cluster
        # and the day index from pt_trace where that patient first tested 1.5
        convert_pt_ids <- unique(seq2pt[cluster_seqs[!(cluster_seqs %in% ip_pt_seqs)]])
        convert_pos <- sapply(convert_pt_ids, function(pt_id) min(dates[cluster_seqs[seq2pt[cluster_seqs] == pt_id]]))

        # In pt_trace, 1.5 presumably indicates a positive test day
        # conversion_day is the row name of pt_trace for that day index
        # If NA, we set it to the last row in pt_trace
        conversion_day_indices <- sapply(convert_pt_ids, function(pt_id) min(which(pt_trace[, pt_id] == 1.5)))
        conversion_day_indices[is.infinite(conversion_day_indices)] <- NA

        conversion_day <- row.names(pt_trace)[conversion_day_indices]
        conversion_day[is.na(conversion_day)] <- row.names(pt_trace)[nrow(pt_trace)]

        # Distinguish 'initial' converts from 'later' converts based on a
        # 7-day difference threshold between the earliest cluster date and that day in trace
        day_diff <- as.Date(convert_pos, origin = "1899-12-30") - as.Date(conversion_day)
        initial_mask <- day_diff <= 7

        initial_convert_pts <- convert_pt_ids[initial_mask]
        convert_pts <- convert_pt_ids

        cluster_prop["Number_of_initial_converts"] <- sum(initial_mask)
        cluster_prop["Number_of_later_converts"] <- sum(!initial_mask)

        # Mark those that fail the 7-day threshold as Inf in convert_pos
        # to indicate no in-cluster first positive
        convert_pos[day_diff > 7] <- Inf
    }

    # 6. Conversions after an index patient is in the cluster #####
    # Count how many converts happen after the earliest known index's date
    if (cluster_prop["Number_of_start_indexes"] > 0 && cluster_prop["Number_of_converts"] > 0) {
        # The earliest index date among ip_seqs in this cluster
        index_pos <- min(dates[ip_seqs[ip_seqs %in% cluster_seqs]])

        # Number of all converts that appear on or after that index date
        cluster_prop["Number_of_converts_after_index"] <- sum(convert_pos >= index_pos)
        cluster_prop["Number_of_initial_converts_after_index"] <-
            sum(convert_pos >= index_pos & convert_pos != Inf)

        # Time from earliest index to first and last convert (for initial converts only)
        valid_idx <- convert_pos[convert_pos != Inf]
        valid_idx <- valid_idx[valid_idx >= index_pos]

        if (length(valid_idx) > 0) {
            cluster_prop["Time_from_index_to_first_convert"] <- min(valid_idx) - index_pos
            cluster_prop["Time_from_index_to_last_convert"] <- max(valid_idx) - index_pos
        } else {
            cluster_prop["Time_from_index_to_first_convert"] <- NA
            cluster_prop["Time_from_index_to_last_convert"] <- NA
        }
    } else {
        cluster_prop["Number_of_converts_after_index"] <- 0
        cluster_prop["Number_of_initial_converts_after_index"] <- 0
        cluster_prop["Time_from_index_to_first_convert"] <- NA
        cluster_prop["Time_from_index_to_last_convert"] <- NA
    }

    # Conversions after ANY intake-positive isolate in the cluster
    # Might matter if a patient has multiple isolates, etc.
    if ((sum(cluster_seqs %in% ip_pt_seqs) > 0) && cluster_prop["Number_of_converts"] > 0) {
        # Find earliest day index in pt_trace for the intake positives
        # row.names(pt_trace) is used to get the day label
        index_pos_indices <- sapply(unique(seq2pt[cluster_seqs[cluster_seqs %in% ip_pt_seqs]]), function(pt_id) {
            min(which(pt_trace[, pt_id] == 1.5))
        })
        index_pos_indices[is.infinite(index_pos_indices)] <- NA

        # Convert the earliest index position to a row name date
        earliest_index_pos <- min(row.names(pt_trace)[index_pos_indices], na.rm = TRUE)
        if (is.infinite(earliest_index_pos)) earliest_index_pos <- NA

        # Count how many convert_pos are >= that date
        if (!is.na(earliest_index_pos)) {
            # Convert earliest_index_pos to Date
            earliest_index_date <- as.Date(earliest_index_pos)
            convert_dates <- as.Date(convert_pos, origin = "1899-12-30")

            cluster_prop["Number_of_converts_after_index_seq"] <- sum(convert_dates >= earliest_index_date)
            cluster_prop["Number_of_initial_converts_after_index_seq"] <-
                sum(convert_dates >= earliest_index_date & convert_pos != Inf)
        } else {
            cluster_prop["Number_of_converts_after_index_seq"] <- 0
            cluster_prop["Number_of_initial_converts_after_index_seq"] <- 0
        }
    } else {
        cluster_prop["Number_of_converts_after_index_seq"] <- 0
        cluster_prop["Number_of_initial_converts_after_index_seq"] <- 0
    }

    # 7. Location overlap checks #####
    # Helper function: given a donor pt_donor, a recipient pt_recipient, earliest_pos_by_pt,
    # and a trace data set, check if there's an overlap in location
    # after pt_donor has the strain and on or before pt_recipient's first positive.
    pt_overlap <- function(pt_donor, pt_recipient, first_pos_vec, trace) {
        donor_date <- first_pos_vec[pt_donor]
        recipient_date <- first_pos_vec[pt_recipient]
        # If the donor's earliest date is <= the recipient's earliest date:
        if (first_pos_vec[pt_donor] <= first_pos_vec[pt_recipient]) {
            time_range <- seq(from = donor_date, to = recipient_date, by = 1)
            time_range <- as.character(as.Date(time_range, origin = "1899-12-30"))
            # Check if there's a day where pt_donor and pt_recipient share a location
            # (i.e., same floor/room/patient unit) in 'trace'
            any(
                floor(trace[time_range, pt_donor]) == floor(trace[time_range, pt_recipient]) &
                    (floor(trace[time_range, pt_recipient]) > 0)
            )
        } else {
            FALSE
        }
    }

    find_source <- function(pt_convert, trace) {
        all_others <- setdiff(unique(seq2pt[cluster_seqs]), pt_convert)
        any(sapply(all_others, pt_overlap, pt_recipient = pt_convert,
                   first_pos_vec = earliest_pos_by_pt, trace = trace))
    }

    # Count how many converts can be assigned a source patient in the cluster
    # based on shared location prior to convert's first positive.
    if (cluster_prop["Number_of_converts"] > 0 &&
            cluster_prop["Number_of_patients"] > 1 &&
            !any(is.infinite(earliest_pos_by_pt))) {
        # convert_pts was captured above as all convert patient IDs
        # for each convert, see if any other patient in the cluster can be a source
        cluster_prop["Number_of_converts_with_source"] <-
            sum(sapply(convert_pts, find_source, trace = pt_trace))
        # If we have floor_trace, do the same overlap check
        if (!is.null(floor_trace)) {
            cluster_prop["Number_of_converts_with_floor_source"] <-
                sum(sapply(convert_pts, find_source, trace = floor_trace))
        }
        # If we have room_trace, do the same overlap check
        if (!is.null(room_trace)) {
            cluster_prop["Number_of_converts_with_room_source"] <-
                sum(sapply(convert_pts, find_source, trace = room_trace))
        }
    } else {
        cluster_prop["Number_of_converts_with_source"] <- 0
        cluster_prop["Number_of_converts_with_floor_source"] <- 0
        cluster_prop["Number_of_converts_with_room_source"] <- 0
    }

    # If we have any initial converts, do the same overlap checks specifically for those patients
    if (cluster_prop["Number_of_initial_converts"] > 0 &&
            cluster_prop["Number_of_patients"] > 1 &&
            !any(is.infinite(earliest_pos_by_pt))) {
        cluster_prop["Number_of_initial_converts_with_source"] <-
            sum(sapply(initial_convert_pts, find_source, trace = pt_trace))
        if (!is.null(floor_trace)) {
            cluster_prop["Number_of_initial_converts_with_floor_source"] <-
                sum(sapply(initial_convert_pts, find_source, trace = floor_trace))
        }
        if (!is.null(room_trace)) {
            cluster_prop["Number_of_initial_converts_with_room_source"] <-
                sum(sapply(initial_convert_pts, find_source, trace = room_trace))
        }
    } else {
        cluster_prop["Number_of_initial_converts_with_source"] <- 0
        cluster_prop["Number_of_initial_converts_with_floor_source"] <- 0
        cluster_prop["Number_of_initial_converts_with_room_source"] <- 0
    }

    # Replace any infinite values with NA
    cluster_prop[is.infinite(cluster_prop)] <- NA

    # Return the final numeric vector of cluster properties
    cluster_prop
}

#' Permutation Test for Cluster Properties
#'
#' @description
#' This function compares the properties of actual clusters to properties of randomly
#' permuted clusters with the same size distribution. It computes various summary statistics.
#'
#' @param clusters A named numeric vector where names are sequence IDs and values are subtrees defining the cluster.
#' @param pt_trace A matrix or data frame with rows representing days and columns representing patients.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_pt_seqs A vector of sequence IDs corresponding to intake positive patients.
#' @param ip_seqs A vector of sequence IDs corresponding to intake patient sequences presumed to be imported.
#' @param dates A vector of isolate dates named by sequence IDs.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param prefix A descriptor string (used to name output figures).
#' @param nperm The number of permutations to perform (default is 1000).
#' @param floor_trace An optional floor trace (rows: days, columns: patients).
#' @param room_trace An optional room trace (rows: days, columns: patients).
#' @param mc.cores The number of cores to use for parallel processing (default is one less than the total
#'                 number of cores). This is only for *nix systems and does not work on Windows - for Windows,
#'                 set mc.cores = 1.
#'
#' @return A matrix of permutation statistics. Rows correspond to each permutation (with the last row
#' containing the observed statistic) and columns correspond to the computed properties.
#'
#' @importFrom parallel detectCores mclapply
#' @export
cluster_property_perm_test <- function(clusters, pt_trace, seq2pt, ip_pt_seqs, ip_seqs, dates, snp_dist,
                                       prefix, nperm = 1000, floor_trace = NULL, room_trace = NULL,
                                       mc.cores = detectCores() - 1) {
    # Process clusters to set single-patient clusters to 1
    unique_clusters <- sort(unique(clusters))
    cluster_size <- sapply(unique_clusters, function(x) {
        length(unique(seq2pt[names(clusters)[clusters == x]]))
    })
    single_clusters <- unique_clusters[cluster_size == 1]
    clusters[clusters %in% single_clusters] <- 1

    # Exclude patients in default cluster 1
    clusters <- clusters[clusters != 1]
    cluster_names <- names(clusters)

    # Compute properties of the actual (observed) clusters
    valid_clusters <- sort(unique(clusters))
    cluster_props <- t(sapply(valid_clusters, function(c) {
        cluster_properties(
            cluster_names[clusters == c],
            pt_trace, seq2pt, ip_pt_seqs, ip_seqs, snp_dist,
            dates, floor_trace, room_trace
        )
    }))
    row.names(cluster_props) <- valid_clusters

    # Create list of eligible index patients and convert patients to be added to clusters
    # Note that this accounts for patients being in multiple clusters (those patients are
    # eligible to go into multiple clusters) and multiple isolates from the same patient
    # being in a single cluster (those isolates always stay together)
    index_pt_start_seqs <- unlist(sapply(ip_seqs, function(seq_id) {
        cluster_names[clusters == clusters[seq_id] & seq2pt[cluster_names] == seq2pt[seq_id]]
    }))
    index_seqs_start <- intersect(index_pt_start_seqs, cluster_names)
    index_seqs_not_start <- intersect(setdiff(ip_pt_seqs, index_seqs_start), cluster_names)
    convert_seqs <- setdiff(cluster_names, c(index_seqs_start, index_seqs_not_start))

    # Create eligibility matrices
    get_elig_mat <- function(seqs) {
        cbind(
            seq = as.character(seqs),
            patient = seq2pt[seqs],
            cluster = clusters[seqs],
            comb = paste(seq2pt[seqs], clusters[seqs], sep = "-")
        )
    }
    elig_index_start_pts <- get_elig_mat(index_seqs_start)
    elig_index_not_start_pts <- get_elig_mat(index_seqs_not_start)
    elig_convert_pts <- get_elig_mat(convert_seqs)

    # Determine the number of clusters each patient is in
    cluster_per_patient <- function(elig_mat) {
        pt_unique_sorted <- sort(unique(elig_mat[, "patient"]))
        setNames(sapply(pt_unique_sorted, function(pt) {
            length(unique(elig_mat[elig_mat[, "patient"] == pt, "cluster"]))
        }), pt_unique_sorted)
    }
    cluster_per_index_start <- cluster_per_patient(elig_index_start_pts)
    cluster_per_index_not_start <- cluster_per_patient(elig_index_not_start_pts)
    cluster_per_convert <- cluster_per_patient(elig_convert_pts)

    # Determine the number of index and convert patients in each cluster
    pt_per_cluster <- function(elig_mat) {
        clusters_unique_sorted <- sort(unique(clusters))
        setNames(sapply(clusters_unique_sorted, function(clust) {
            length(unique(elig_mat[elig_mat[, "cluster"] == clust, "patient"]))
        }), clusters_unique_sorted)
    }
    index_start_per_cluster <- pt_per_cluster(elig_index_start_pts)
    index_not_start_per_cluster <- pt_per_cluster(elig_index_not_start_pts)
    convert_per_cluster <- pt_per_cluster(elig_convert_pts)

    # Permute clusters and compute their properties
    n_props <- ncol(cluster_props)
    prop_array <- array(
        dim = c(nrow(cluster_props), n_props, nperm + 1),
        dimnames = list(
            seq_len(nrow(cluster_props)),
            colnames(cluster_props),
            seq_len(nperm + 1)
        )
    )

    # Assign patients to clusters, starting with patients in multiple clusters
    assign_pt_clusters <- function(elig_vec, cluster_per_pt, pt_per_cluster, rand_clusters) {
        names_pt_per_cluster <- names(pt_per_cluster)
        # For each patient, assign them to a cluster
        for (pt in names(sort(cluster_per_pt, decreasing = TRUE))) {
            # Get unique cluster assignments for the current patient
            pt_clust_ids <- unique(elig_vec[elig_vec[, "patient"] == pt, "comb"])
            # Pick correct number of eligible clusters to assign to them
            cluster_assign <- sample(names_pt_per_cluster[pt_per_cluster > 0], length(pt_clust_ids))
            # Remove one from available openings in assigned
            pt_per_cluster[cluster_assign] <- pt_per_cluster[cluster_assign] - 1
            # For each patient-cluster, assign all associated isolates to new cluster
            for (pt_c_i in seq_along(pt_clust_ids)) {
                # Get all sequences for this patient-cluster combination
                seqs_to_assign <- elig_vec[elig_vec[, "comb"] == pt_clust_ids[pt_c_i], "seq"]
                # Assign all sequences to the target cluster
                rand_clusters[seqs_to_assign] <- as.numeric(cluster_assign[pt_c_i])
            }
        }
        rand_clusters
    }

    # Sequentially assign patients to clusters from each category
    rand_clusters <- assign_pt_clusters(
        elig_convert_pts, cluster_per_convert,
        convert_per_cluster, rand_clusters
    ) # Converts
    rand_clusters <- assign_pt_clusters(
        elig_index_start_pts, cluster_per_index_start,
        index_start_per_cluster, rand_clusters
    ) # Index patients who started clusters
    rand_clusters <- assign_pt_clusters(
        elig_index_not_start_pts, cluster_per_index_not_start,
        index_not_start_per_cluster, rand_clusters
    ) # Index patients whose isolate was not an index isolate

    rand_clusters
}

#' Permutation Test for Cluster Properties
#'
#' This function compares the properties of actual clusters to properties of randomly
#' permuted clusters with the same size distribution. It computes various summary statistics,
#' calculates permutation-based p-values, and generates histogram plots of the statistics.
#'
#' @param clusters A named numeric vector where names are sequence IDs and values are subtrees defining the cluster.
#' @param pt_trace A matrix or data frame with rows representing days and columns representing patients.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_pt_seqs A vector of sequence IDs corresponding to intake positive patients.
#' @param ip_seqs A vector of sequence IDs corresponding to intake patient sequences presumed to be imported.
#' @param dates A vector of isolate dates named by sequence IDs.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param prefix A descriptor string (used to name output figures).
#' @param nperm The number of permutations to perform (default is 1000).
#' @param floor_trace An optional floor trace (rows: days, columns: patients).
#' @param room_trace An optional room trace (rows: days, columns: patients).
#' @param mc.cores The number of cores to use for parallel processing (default is one less than the total
#'                 number of cores). This is only for *nix systems and does not work on Windows - for Windows,
#'                 set mc.cores = 1.
#'
#' @return A matrix of permutation statistics. Rows correspond to each permutation (with the last row
#' containing the observed statistic) and columns correspond to the computed properties.
#'
#' @importFrom parallel detectCores mclapply
#' @export
cluster_property_perm_test <- function(clusters, pt_trace, seq2pt, ip_pt_seqs, ip_seqs, dates, snp_dist,
                                       prefix, nperm = 1000, floor_trace = NULL, room_trace = NULL,
                                       mc.cores = detectCores() - 1) {
    # Process clusters to set single-patient clusters to 1
    unique_clusters <- sort(unique(clusters))
    cluster_size <- sapply(unique_clusters, function(x) {
        length(unique(seq2pt[names(clusters)[clusters == x]]))
    })
    names(cluster_size) <- unique_clusters
    single_clusters <- names(cluster_size)[cluster_size == 1]
    clusters[clusters %in% single_clusters] <- 1

    # Compute properties of the actual (observed) clusters
    valid_clusters <- sort(unique(clusters[clusters != 1]))
    cluster_props <- t(sapply(valid_clusters, function(c) {
        cluster_properties(
            names(clusters)[clusters == c],
            pt_trace, seq2pt, ip_pt_seqs, ip_seqs, snp_dist,
            dates, floor_trace, room_trace
        )
    }))
    row.names(cluster_props) <- valid_clusters

    # Permute clusters and compute their properties
    n_props <- ncol(cluster_props)
    prop_array <- array(
        dim = c(nrow(cluster_props), n_props, nperm + 1),
        dimnames = list(
            seq_len(nrow(cluster_props)),
            colnames(cluster_props),
            seq_len(nperm + 1)
        )
    )
    perm_results <- mclapply(seq_len(nperm), function(n) {
        perm_clust <- permute_clusters(clusters, seq2pt, ip_seqs, ip_pt_seqs, all = FALSE)
        perm_valid <- sort(unique(perm_clust))
        perm_cluster_props <- t(sapply(perm_valid, function(c) {
            cluster_properties(
                names(perm_clust)[perm_clust == c],
                pt_trace, seq2pt, ip_pt_seqs, ip_seqs,
                snp_dist, dates, floor_trace, room_trace
            )
        }))
        perm_cluster_props
    }, mc.cores = mc.cores)

    # Store the permuted cluster properties in the array
    for (n in seq_len(nperm)) {
        prop_array[, , n] <- perm_results[[n]]
    }
    # Store observed cluster properties in the final slice of the array
    prop_array[, , nperm + 1] <- cluster_props

    # Compute summary statistics for each permutation
    median_stat_names <- paste("Median", colnames(cluster_props))

    per_conv_cols <- c(
        "Number_of_converts_with_source", "Number_of_converts_with_floor_source",
        "Number_of_converts_with_room_source", "Number_of_converts_after_index"
    )
    per_conv_stat_names <- paste("Per convert", per_conv_cols)

    per_init_conv_cols <- c(
        "Number_of_initial_converts_with_source", "Number_of_initial_converts_with_floor_source",
        "Number_of_initial_converts_with_room_source", "Number_of_initial_converts_after_index"
    )
    per_init_conv_stat_names <- paste("Per convert", per_init_conv_cols)

    cluster_cols <- c(
        "Number_of_converts_with_source", "Number_of_converts_with_floor_source",
        "Number_of_converts_with_room_source", "Number_of_converts_after_index"
    )
    cluster_stat_names <- paste("Fraction of clusters", cluster_cols)

    prop_stat_names <- c(
        median_stat_names, per_conv_stat_names,
        per_init_conv_stat_names, cluster_stat_names
    )

    perm_props <- matrix(
        nrow = nperm + 1, ncol = length(prop_stat_names),
        dimnames = list(seq_len(nperm + 1), prop_stat_names)
    )

    for (i in seq_len(nperm + 1)) {
        # Median statistics
        perm_props[i, median_stat_names] <- apply(prop_array[, , i], 2, median, na.rm = TRUE)

        # Per-convert statistics (as a fraction of the total converts in observed clusters)
        perm_props[i, per_conv_stat_names] <- sapply(per_conv_cols, function(col) {
            round(sum(prop_array[, col, i]) / sum(cluster_props[, "Number_of_converts"]), 2)
        })

        # Per-initial-convert statistics
        perm_props[i, per_init_conv_stat_names] <- sapply(per_init_conv_cols, function(col) {
            round(sum(prop_array[, col, i]) / sum(cluster_props[, "Number_of_initial_converts"]), 2)
        })

        # Fraction of clusters statistic
        perm_props[i, cluster_stat_names] <- sapply(cluster_cols, function(col) {
            round(sum(prop_array[, col, i] == cluster_props[, "Number_of_converts"]) / nrow(cluster_props), 2)
        })
    }

    perm_props
}
