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
#'
#' @return A data frame with one row per cluster (plus a row for the overall population) containing the
#'   number of variable sites and the rates for six mutation types. Mutation rates are rounded to three decimals.
#'
#' @export
intra_cluster_genetic_var_analysis <- function(clusters, dna_aln) {
    # pre-compute a character matrix of the alignment
    dna_aln_char <- matrix(
        as.character(dna_aln),
        nrow = nrow(dna_aln),
        dimnames = dimnames(dna_aln)
    )

    # identify cluster IDs (excluding a potential singleton "1" if present)
    cluster_ids <- setdiff(sort(unique(clusters)), 1)

    # compute major alleles across the entire alignment
    major_alleles_all <- apply(dna_aln_char, 2, function(pos_aln) {
        allele_tab <- table(pos_aln)
        names(allele_tab)[which.max(allele_tab)]
    })

    # get all variable positions in the alignment
    var_pos <- apply(dna_aln_char, 2, function(x) sum(x == x[1]) < nrow(dna_aln_char))

    # filter the alignment to only include variable positions
    dna_aln_char <- dna_aln_char[, var_pos, drop = FALSE]
    major_alleles_all <- major_alleles_all[var_pos]

    # pre-allocate vectors for results
    n_iter <- length(cluster_ids) + 1 # +1 for the overall population
    num_var_sites <- numeric(n_iter)
    gc_at_transition_rate <- numeric(n_iter)
    at_gc_transition_rate <- numeric(n_iter)
    gc_ta_transversion_rate <- numeric(n_iter)
    at_cg_transversion_rate <- numeric(n_iter)
    at_ta_transversion_rate <- numeric(n_iter)
    gc_cg_transversion_rate <- numeric(n_iter)

    # define row names for the output: overall population and each cluster
    cluster_labels <- c("Pop. freq.", cluster_ids)

    iter <- 1
    # loop over overall population (represented by -1) and each cluster
    for (cluster in c(-1, cluster_ids)) {
        seq_ids <- if (cluster == -1) {
            names(clusters)
        } else {
            names(clusters)[clusters == cluster]
        }

        # subset the alignment for the current group
        aln_subset <- dna_aln_char[seq_ids, , drop = FALSE]
        # determine variable positions within the current group
        var_pos_cluster <- apply(aln_subset, 2, function(pos_aln) {
            sum(pos_aln == pos_aln[1]) < length(seq_ids)
        })
        var_idx <- which(var_pos_cluster)

        # compute minor allele for each variable position
        minor_alleles <- vapply(
            var_idx,
            function(j) {
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
            },
            character(1)
        )
        # compute major allele for each variable position
        major_alleles <- major_alleles_all[var_idx]

        # filter positions to those with valid nucleotide calls (both alleles among a, c, t, g)
        valid_pos <- which(
            minor_alleles != major_alleles &
                major_alleles %in% c("a", "c", "t", "g") &
                minor_alleles %in% c("a", "c", "t", "g")
        )

        num_sites <- length(valid_pos)
        num_var_sites[iter] <- num_sites

        # if no positions are variable or no valid sites are found, assign NA to mutation rates
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

        # calculate mutation rates
        gc_at_transition_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "a") |
                (major_valid == "c" & minor_valid == "t")
        ) /
            num_sites
        at_gc_transition_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "g") |
                (major_valid == "t" & minor_valid == "c")
        ) /
            num_sites
        gc_ta_transversion_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "t") |
                (major_valid == "c" & minor_valid == "a")
        ) /
            num_sites
        at_cg_transversion_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "c") |
                (major_valid == "t" & minor_valid == "g")
        ) /
            num_sites
        at_ta_transversion_rate[iter] <- sum(
            (major_valid == "a" & minor_valid == "t") |
                (major_valid == "t" & minor_valid == "a")
        ) /
            num_sites
        gc_cg_transversion_rate[iter] <- sum(
            (major_valid == "g" & minor_valid == "c") |
                (major_valid == "c" & minor_valid == "g")
        ) /
            num_sites

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

#' Compute Genetic Distances
#'
#' @description
#' Calculate genetic distance statistics within a cluster
#'
#' @param seqs Vector of sequence IDs in the cluster
#' @param snp_dist Matrix of SNP distances between isolates
#'
#' @return List with genetic distance properties
#' @export
intra_cluster_genetic_distances <- function(seqs, snp_dist) {
    # subset the SNP distance matrix to only the cluster sequences
    snp_dist_subset <- snp_dist[seqs, seqs, drop = FALSE]

    # get the upper triangle of the SNP distance matrix (excluding the diagonal)
    seq_vals <- snp_dist_subset[upper.tri(snp_dist_subset)]

    genetic_stats <- numeric(3)
    names(genetic_stats) <- c(
        "mean_genetic_distance",
        "median_genetic_distance",
        "max_genetic_distance"
    )

    genetic_stats["mean_genetic_distance"] <- mean(seq_vals)
    genetic_stats["median_genetic_distance"] <- median(seq_vals)
    genetic_stats["max_genetic_distance"] <- max(seq_vals)

    genetic_stats
}

#' Compute Cluster Duration Metrics
#'
#' @description
#' Calculate temporal metrics for cluster spread
#'
#' @param seqs Vector of sequence IDs in the cluster
#' @param seq2pt Named vector mapping sequence IDs to patient IDs
#' @param dates Vector of isolate dates named by sequence IDs
#'
#' @return List with cluster duration properties and earliest positive date for each patient
#' @export
compute_cluster_duration_metrics <- function(seqs, seq2pt, dates) {
    # calculate earliest positive date for each patient in the cluster
    earliest_pos_by_pt <- vapply(
        unique(seq2pt[seqs]),
        function(pt_id) {
            patient_seqs <- seqs[seq2pt[seqs] == pt_id]
            if (length(patient_seqs) > 0) {
                min(dates[patient_seqs])
            } else {
                Inf
            }
        },
        numeric(1)
    )
    earliest_pos_by_pt <- sort(setNames(earliest_pos_by_pt, unique(seq2pt[seqs])))

    duration_stats <- numeric(3)
    names(duration_stats) <- c(
        "time_to_first_acquisition",
        "time_to_last_acquisition",
        "median_time_to_acquisition"
    )

    # check if we have any patients
    if (length(earliest_pos_by_pt) > 0) {
        # time to first acquisition = difference between the earliest two patients' dates
        if (length(earliest_pos_by_pt) >= 2) {
            duration_stats["time_to_first_acquisition"] <- earliest_pos_by_pt[2] -
                earliest_pos_by_pt[1]
        } else {
            duration_stats["time_to_first_acquisition"] <- 0
        }

        # time to last acquisition = difference between the earliest and latest
        duration_stats["time_to_last_acquisition"] <-
            earliest_pos_by_pt[length(earliest_pos_by_pt)] - earliest_pos_by_pt[1]

        # median time to acquisition
        if (length(earliest_pos_by_pt) >= 2) {
            duration_stats["median_time_to_acquisition"] <- median(vapply(
                earliest_pos_by_pt[2:length(earliest_pos_by_pt)],
                function(x) x - earliest_pos_by_pt[1],
                numeric(1)
            ))
        } else {
            duration_stats["median_time_to_acquisition"] <- 0
        }
    } else {
        # no patients in cluster
        duration_stats["time_to_first_acquisition"] <- 0
        duration_stats["time_to_last_acquisition"] <- 0
        duration_stats["median_time_to_acquisition"] <- 0
    }

    list(stats = duration_stats, earliest_pos_by_pt = earliest_pos_by_pt)
}
