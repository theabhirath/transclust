#' Analyze Intra-Cluster Genetic Variation
#'
#' This function analyzes genetic diversity within clusters by identifying variant sites and calculating
#' mutational spectra (e.g. transition and transversion rates). It uses a DNA alignment and a vector
#' of cluster assignments (named by sequence IDs) and returns a data frame summarizing the diversity
#' measures for each cluster as well as for the overall population.
#'
#' @param clusters A named vector of cluster assignments (cluster IDs) with names corresponding to sequence IDs.
#' @param dna_aln A matrix (or data frame) representing the DNA sequence alignment, where rows are sequences
#'   (named by sequence IDs) and columns represent nucleotide positions.
#'
#' @return A data frame with one row per cluster (plus a row for the overall population) containing the
#'   number of variable sites and the rates for six mutation types. Mutation rates are rounded to two decimals.
intra_cluster_genetic_var_analysis <- function(clusters, dna_aln) {
    # Identify cluster IDs (excluding a potential singleton "1" if present)
    cluster_ids <- setdiff(sort(unique(clusters)), 1)

    # Compute major alleles across the entire alignment (as characters)
    aln_chars <- as.character(dna_aln)
    major_alleles_all <- apply(aln_chars, 2, function(pos_aln) {
        allele_tab <- table(pos_aln)
        names(sort(allele_tab, decreasing = TRUE))[1]
    })

    # Identify positions that are variable among the population sequences
    var_pos_all <- apply(aln_chars, 2, function(pos_aln) {
        sum(pos_aln == pos_aln[1]) < nrow(dna_aln)
    })
    major_alleles_all <- major_alleles_all[var_pos_all]

    # Prepare pre-allocated vectors for the results.
    n_iter <- length(cluster_ids) + 1 # +1 for the overall population
    num_var_sites <- numeric(n_iter)
    gc_at_transition_rate <- numeric(n_iter)
    at_gc_transition_rate <- numeric(n_iter)
    gc_ta_transversion_rate <- numeric(n_iter)
    at_cg_transversion_rate <- numeric(n_iter)
    at_ta_transversion_rate <- numeric(n_iter)
    gc_cg_transversion_rate <- numeric(n_iter)

    # Define row names for the output: overall population and each cluster.
    cluster_labels <- c("Population frequency", cluster_ids)

    iter <- 1
    # Loop over overall population (represented by -1) and each cluster
    for (cluster in c(-1, cluster_ids)) {
        if (cluster == -1) {
            seq_ids <- names(clusters)
        } else {
            seq_ids <- names(clusters)[clusters == cluster]
        }

        # Subset the alignment for the current group and determine variable positions within this group.
        aln_subset <- as.character(dna_aln[seq_ids, , drop = FALSE])
        var_pos_cluster <- apply(aln_subset, 2, function(pos_aln) {
            sum(pos_aln == pos_aln[1]) < length(seq_ids)
        })


        # Get indices of variable positions in this cluster.
        var_idx <- which(var_pos_cluster)

        # For each variable position, determine the "minor allele" call for the cluster.
        minor_alleles <- sapply(var_idx, function(j) {
            pos_aln <- as.character(dna_aln[, j])
            allele_tab <- table(pos_aln)
            major_allele <- names(sort(allele_tab, decreasing = TRUE))[1]
            cluster_alleles <- unique(aln_subset[, j])

            if (length(cluster_alleles) == 1 && major_allele %in% cluster_alleles) {
                "X" # Cluster has the major allele
            } else if (length(cluster_alleles) == 1) {
                "Y" # Cluster has a minor allele shared by all members
            } else if (length(cluster_alleles) == 2 && major_allele %in% cluster_alleles) {
                setdiff(cluster_alleles, major_allele)  # Cluster has a variable minor allele
            } else {
                "Z" # Cluster has multiple minor alleles
            }
        })

        # Get the corresponding major alleles from the overall set.
        maj_alleles <- major_alleles_all[var_idx]

        # Filter positions to those with valid nucleotide calls (both alleles among a, c, t, g)
        valid_pos <- which((minor_alleles != maj_alleles) & 
                            (maj_alleles %in% c("a", "c", "t", "g")) &
                            (minor_alleles %in% c("a", "c", "t", "g")))

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

        # Extract valid alleles for computation.
        maj_valid <- maj_alleles[valid_pos]
        minor_valid <- minor_alleles[valid_pos]

        # Calculate mutation rates.
        gc_at_transition_rate[iter] <- sum(
            (maj_valid == "g" & minor_valid == "a") | (maj_valid == "c" & minor_valid == "t")
        ) / num_sites
        at_gc_transition_rate[iter] <- sum(
            (maj_valid == "a" & minor_valid == "g") | (maj_valid == "t" & minor_valid == "c")
        ) / num_sites
        gc_ta_transversion_rate[iter] <- sum(
            (maj_valid == "g" & minor_valid == "t") | (maj_valid == "c" & minor_valid == "a")
        ) / num_sites
        at_cg_transversion_rate[iter] <- sum(
            (maj_valid == "a" & minor_valid == "c") | (maj_valid == "t" & minor_valid == "g")
        ) / num_sites
        at_ta_transversion_rate[iter] <- sum(
            (maj_valid == "a" & minor_valid == "t") | (maj_valid == "t" & minor_valid == "a")
        ) / num_sites
        gc_cg_transversion_rate[iter] <- sum(
            (maj_valid == "g" & minor_valid == "c") | (maj_valid == "c" & minor_valid == "g")
        ) / num_sites

        iter <- iter + 1
    }

    # Combine the computed metrics and return the result data frame
    data.frame(
        num_var_sites = num_var_sites,
        GC_AT_transition_rate = round(gc_at_transition_rate, 2),
        AT_GC_transition_rate = round(at_gc_transition_rate, 2),
        GC_TA_transversion_rate = round(gc_ta_transversion_rate, 2),
        AT_TA_transversion_rate = round(at_ta_transversion_rate, 2),
        GC_CG_transversion_rate = round(gc_cg_transversion_rate, 2),
        AT_CG_transversion_rate = round(at_cg_transversion_rate, 2),
        row.names = cluster_labels,
        check.names = FALSE
    )
}
