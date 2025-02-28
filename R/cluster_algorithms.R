#' get_tn_clusters_snpthresh
#'
#' Performs clustering of isolates using a hard SNP distance cutoff.
#'
#' @param dna_aln A DNA object for sequences of interest.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param snp_thresh A threshold for defining clusters.
#'
#' @return A vector indicating the cluster that each isolate belongs to.
#'
#' @export
#'
get_tn_clusters_snp_thresh <- function(dna_aln, snp_dist, snp_thresh) {
    # Perform hierarchical clustering based on SNP distance matrix
    snp_clust <- hclust(as.dist(snp_dist))
    clusters <- cutree(snp_clust, h = snp_thresh)

    # Convert hierarchical clustering into a phylogenetic tree
    upgma_tree <- ape::as.phylo(snp_clust)
    # Extract all subtrees from the phylogenetic tree
    upgma_sub_trees <- ape::subtrees(upgma_tree)

    # Pre-compute the size and labels of each subtree to avoid recomputing
    sub_tree_data <- lapply(upgma_sub_trees, function(st) list(size = length(st$tip.label), labels = st$tip.label))

    # Sequentially assign each cluster to the best matching subtree
    st_clusters <- vapply(unique(clusters), FUN = function(c) {
        # Count how many sequences in the cluster are in each subtree
        sub_tree_match <- vapply(sub_tree_data, function(st) {
            sum(names(clusters)[clusters == c] %in% st$labels)
        }, numeric(1))
        # Get index of best-matching subtree
        sub_tree_match_i <- which(sub_tree_match == max(sub_tree_match))
        # If the cluster has more than one sequence, assign it to the smallest matching subtree
        # else assign it to the default subtree
        if (sum(clusters == c) > 1) {
            smallest_subtree <- which.min(vapply(sub_tree_data[sub_tree_match_i], function(st) st$size, numeric(1)))
            sub_tree_match_i[smallest_subtree]
        } else {
            1
        }
    }, numeric(1))

    # Each cluster is associated with a subtree
    setNames(st_clusters[match(clusters, unique(clusters))], names(clusters))
}

#' Get Transmission Clusters
#'
#' This function identifies transmission clusters based on the number of shared variants.
#'
#' @param dna_aln A DNA object for sequences of interest.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param ip_seqs A vector of sequence IDs which correspond to intake patient sequences
#'                presumed to be imported.
#' @param ip_pt_seqs A vector of sequence IDs which correspond to intake-positive patients.
#' @param seq2pt A named vector linking sequence IDs to patient IDs.
#' @param dates A vector of isolate dates named by sequence IDs.
#' @param pars Logical, if TRUE then use parsimony tree instead of the default neighbor joining tree.
#' @param tree A tree to use if provided instead of the default neighbor joining tree.
#'
#' @return A vector indicating the cluster that each isolate belongs to.
#'
#' @details Clustering is performed to identify the maximal clusters containing a single intake-positive patient
#' that occurs before all cluster converts. The clustering metric is the number of shared variants, and clusters
#' can have multiple intake-positive patients if they share an identical number of variants with other cluster
#' members or intake-positive patients occur after converts. This clustering also requires that clusters be
#' defined by at least one shared variant that other isolates don't have.
#'
#' @export
get_tn_clusters_MSV_SVst_index_first <- function(dna_aln, snp_dist, ip_seqs, ip_pt_seqs,
                                                 seq2pt, dates, pars = FALSE, tree = NULL) {

    #####################################################################################
    # 1. Tree construction and subtree extraction #####
    #   - If a tree is provided, we simply extract its subtrees. This will ignore the pars argument.
    #   - If not provided, we construct a neighbor-joining tree and reroot it using the out-group.
    #   - If pars is TRUE, we create a parsimony tree instead of the default neighbor joining tree.
    #####################################################################################
    # Determine out-group: the isolate with the maximum average SNP distance
    out_group <- which.max(rowMeans(snp_dist))
    if (is.null(tree)) {
        nj_tree <- ape::nj(snp_dist)
        nj_tree <- phytools::reroot(nj_tree, which(nj_tree$tip.label == row.names(snp_dist)[out_group]))
        # If pars is TRUE, create a parsimony tree
        tree <- if (pars) phangorn::optim.parsimony(nj_tree, phangorn::as.phyDat(dna_aln), all = TRUE) else nj_tree
    }
    sub_trees <- ape::subtrees(tree)
    # log completion of phase to standard output
    message("Phase 1 complete: Tree construction and subtree extraction")

    #####################################################################################
    # 2. Compute the shared variant matrix #####
    # For each pair of isolates, we compute the number of positions where:
    #   1. the isolate's base differs from the out-group,
    #   2. the two isolates share the same base,
    #   3. the base is not a gap ('-') or ambiguous ('n'),
    # and then we convert this similarity count to a distance-like measure.
    #####################################################################################
    n_isolates <- nrow(dna_aln)
    isolate_names <- rownames(dna_aln)

    # Initialize an empty matrix to store the computed values
    dna_shared_mat <- matrix(0, nrow = n_isolates, ncol = n_isolates, dimnames = list(isolate_names, isolate_names))

    # Get all unique pairs of isolates
    pair_indices <- combn(n_isolates, 2)
    # Preallocate vector to hold shared variant counts
    shared_counts <- numeric(ncol(pair_indices))

    # Pre-calculate valid positions for the out-group to avoid repeated work
    # The out-group's base must be neither '-' nor 'n' i.e. not ambiguous
    out_valid <- (dna_aln[out_group, ] != "-" & dna_aln[out_group, ] != "n")

    # Loop over each pair to compute the count of shared variants
    for (k in seq_len(ncol(pair_indices))) {
        i <- pair_indices[1, k]
        j <- pair_indices[2, k]
        # Conditions:
        #   1. The base of isolate i must differ from the out-group
        #   2. Isolates i and j must have the same base
        #   3. The base of isolate i must be neither '-' nor 'n' i.e. not ambiguous
        #   4. The out-group's base must be valid
        cond <- (dna_aln[i, ] != dna_aln[out_group, ]) &
            (dna_aln[i, ] == dna_aln[j, ]) &
            (dna_aln[i, ] != "-" & dna_aln[i, ] != "n") &
            out_valid
        shared_counts[k] <- sum(cond)
    }

    # Transform counts into a distance-like metric
    # We subtract each count from the maximum observed count
    max_shared <- max(shared_counts)
    for (k in seq_len(ncol(pair_indices))) {
        i <- pair_indices[1, k]
        j <- pair_indices[2, k]
        value <- max_shared - shared_counts[k]
        dna_shared_mat[i, j] <- value
        dna_shared_mat[j, i] <- value
    }
    # Diagonals are set to infinity
    diag(dna_shared_mat) <- Inf

    # log completion of phase to standard output
    message("Phase 2 complete: Shared variant matrix computation")

    ##############################################################################
    # 3. Identify defining variants in each subtree #####
    # For every subtree, we loop over each alignment column (site) to check:
    #   - All tips in the subtree have the same base (and that base is not 'n').
    #   - Outside the subtree, there is at least one different base.
    #   - No outside isolate has the same base as the subtree.
    # The total count of such sites is stored for each subtree.
    ##############################################################################
    sub_trees_sv <- vapply(sub_trees, function(st) {
        # For each column in the alignment, test if it is a "defining variant" for the subtree.
        sv_pos <- apply(dna_aln, 2, function(column) {
            # Check if all tips in the subtree share the same base and that base is not 'n'.
            same_in_subtree <- (length(unique(column[st$tip.label])) == 1) &&
                (column[st$tip.label[1]] != "n")
            # Get bases for isolates not in the subtree.
            outside <- setdiff(isolate_names, st$tip.label)
            # Outside should have at least one base (ignoring 'n') that is different.
            outside_diff <- length(unique(setdiff(as.character(column[outside]), "n"))) > 0
            # Ensure that none of the outside isolates match the subtree’s base.
            none_match <- sum(as.character(column[outside]) == as.character(column[st$tip.label[1]])) == 0
            same_in_subtree && outside_diff && none_match
        })
        sum(sv_pos > 0)
    }, numeric(1))

    # log completion of phase to standard output
    message("Phase 3 complete: Defining variant identification")

    #####################################################################################
    # 4. Validate and score subtrees as clusters #####
    # For each subtree, we:
    #   - Identify intake-positive sequences that occur before any "convert" in the subtree.
    #   - Choose one representative sequence per patient (the one with the minimum shared variant distance).
    #   - Verify that the subtree has at least one defining variant and that index isolates are not
    #     overly differentiated (i.e. their shared variant counts are nearly identical).
    #   - Compute a cluster score that rewards early convert sequences and penalizes intake-positives
    #     that could not have started a cluster.
    # This gives us subtrees that are valid clusters. Valid clusters are those which have
    # no more than one index that occurs before converts.
    #####################################################################################
    subtree_metrics <- lapply(seq_along(sub_trees), function(st_i) {
        st <- sub_trees[[st_i]]
        # Identify intake-positive sequences in the subtree
        ip_in_sub <- intersect(ip_seqs, st$tip.label)
        # Determine the earliest date among convert sequences (non-intake-positives)
        convert_dates <- dates[setdiff(st$tip.label, ip_pt_seqs)]
        min_convert_date <- if (length(convert_dates) > 0) min(convert_dates) else Inf
        # Filter intake-positives to those that occur on/before the cutoff
        valid_ip <- ip_in_sub[dates[ip_in_sub] <= min_convert_date]
        # Select one representative sequence per unique patient among valid intake-positives
        ip_rep <- vapply(unique(seq2pt[valid_ip]), function(pt) {
            seqs <- intersect(names(seq2pt)[seq2pt == pt], st$tip.label)
            other <- setdiff(st$tip.label, seqs)
            # If there are no other sequences, return the first one
            if (length(other) == 0) {
                seqs[1]
            } else {
                # Otherwise, find the sequence with the minimum shared variant distance to all other sequences
                min_vals <- sapply(seqs, function(seq) min(dna_shared_mat[other, seq]))
                seqs[which.min(min_vals)]
            }
        }, character(1))
        # Get the distance from the representative to all other isolates in the subtree
        shared_counts_rep <- vapply(ip_rep, function(ip) {
            min(dna_shared_mat[setdiff(st$tip.label, ip), ip])
        }, numeric(1))
        # Compute cluster score if there are defining variants and index isolates are not overly
        # distant from the representative sequence (i.e. their shared variant counts are nearly identical)
        score <- if (sub_trees_sv[st_i] > 0 && length(unique(shared_counts_rep)) <= 1) {
            # Reward intake-positives that could have started a cluster
            ip_convert_seq_count <- length(setdiff(
                st$tip.label,
                ip_pt_seqs[seq2pt[ip_pt_seqs] %in% setdiff(unique(seq2pt), seq2pt[unlist(ip_rep)])]
            ))
            # Penalize intake-positives that could not have started a cluster
            ip_pt_count <- sum(st$tip.label %in% ip_pt_seqs) / 10^6
            ip_convert_seq_count - ip_pt_count
        } else {
            0
        }
        # Compute index count: subtrees that have an index before all converts
        num_converts <- length(setdiff(st$tip.label, ip_pt_seqs))
        num_ip_pts <- length(unique(seq2pt[intersect(st$tip.label, ip_seqs)]))
        index_count <- if (length(valid_ip) > 0 && num_converts > 0) num_ip_pts else 0
        # Return both metrics
        list(score = score, index_count = index_count)
    })

    # Extract both metrics in a single step using a matrix conversion, and unlist
    metrics_mat <- do.call(rbind, subtree_metrics)
    sub_trees_valid <- unlist(metrics_mat[, "score"])
    sub_trees_index_first <- unlist(metrics_mat[, "index_count"])

    # log completion of phase to standard output
    message("Phase 4 complete: Cluster validation and scoring")

    #####################################################################################
    # 5. Assign clusters to each isolate #####
    # For every isolate:
    #  - Check all valid subtrees (i.e. those with positive scores) that contain it.
    #  - Ensure that a candidate subtree is not "nested" within another subtree that
    #    has an earlier index (i.e. a lower index_count).
    #  - Assign the cluster based on the subtree with the highest score.
    #  - If no valid subtree is found, default to cluster 1.
    #  - Return a vector of clusters for each isolate.
    #####################################################################################
    clusters <- vapply(isolate_names, function(x) {
        st_scores <- vapply(which(sub_trees_valid > 0), function(st_i) {
            st <- sub_trees[[st_i]]
            # Skip if the isolate is not in the subtree
            if (!(x %in% st$tip.label)) {
                return(0)
            }
            # Check for nested clusters:
            nested <- vapply(setdiff(which(sub_trees_valid > 0), st_i), function(sub_i) {
                as.integer(all(sub_trees[[sub_i]]$tip.label %in% st$tip.label) &&
                               sub_trees_index_first[sub_i] < sub_trees_index_first[st_i] &&
                               sub_trees_index_first[sub_i] > 0)
            }, integer(1))
            if (sum(nested) == 0) sub_trees_valid[[st_i]] else 0
        }, numeric(1))
        # If no valid subtrees are found, default to cluster 1
        # Otherwise, assign the cluster based on the subtree with the highest score
        if (length(st_scores) == 0 || max(st_scores) <= 0) 1 else which(sub_trees_valid > 0)[which.max(st_scores)]
    }, numeric(1))

    # log completion of phase to standard output
    message("Phase 5 complete: Cluster assignment")

    # Return the clusters
    clusters
}
