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
#' @importFrom ape as.phylo subtrees
#' @importFrom stats as.dist hclust cutree setNames
#'
#' @export
get_tn_clusters_snp_thresh <- function(dna_aln, snp_dist, snp_thresh) {
    # Perform hierarchical clustering based on SNP distance matrix
    snp_clust <- hclust(as.dist(snp_dist))
    clusters <- cutree(snp_clust, h = snp_thresh)

    # Convert hierarchical clustering into a phylogenetic tree using the UPGMA method
    upgma_tree <- as.phylo(snp_clust)
    # Extract all subtrees from the phylogenetic tree
    upgma_sub_trees <- subtrees(upgma_tree)

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
#' @param method A string indicating the method to use for tree construction. Options are
#'               "nj" (neighbor-joining) or "pars" (maximum parsimony). Default is "pars".
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
#' @importFrom ape subtrees
#' @export
get_tn_clusters_MSV_SVst_index_first <- function(dna_aln, snp_dist, ip_seqs, ip_pt_seqs, seq2pt,
                                                 dates, method = "pars", tree = NULL) {
    #####################################################################################
    # 1. Construct the phylogenetic tree #####
    # Get subtrees from the phylogenetic tree
    #####################################################################################
    # If a tree is provided, use it. Otherwise, construct a new tree using the DNA alignment.
    tree <- if (is.null(tree)) get_phylo_tree(dna_aln, snp_dist, method)
    sub_trees <- subtrees(tree)
    # Log completion of phase to standard output
    message("Phase 1 complete: Tree construction and subtree extraction")

    #####################################################################################
    # 2. Compute the shared variant matrix #####
    # For each pair of isolates, we compute the number of positions where:
    #   a. the isolate's base differs from the out-group,
    #   b. the two isolates share the same base,
    #   c. the base is not a gap ('-') or ambiguous ('n'),
    # and then we convert this similarity count to a distance-like measure.
    ####################################################################################
    # Call the Rcpp function to compute the shared variant matrix
    out_group <- which.max(rowMeans(snp_dist))
    dna_shared_mat <- computeSharedMatrix(dna_aln, out_group)

    isolate_names <- row.names(dna_aln)
    rownames(dna_shared_mat) <- isolate_names
    colnames(dna_shared_mat) <- isolate_names

    # Log completion of phase to standard output
    message("Phase 2 complete: Shared variant matrix computation")

    ##############################################################################
    # 3. Identify defining variants in each subtree #####
    # For every subtree, we loop over each alignment column (site) to check:
    #   - All tips in the subtree have the same base (and that base is not 'n').
    #   - Outside the subtree, there is at least one different base.
    #   - No outside isolate has the same base as the subtree.
    # The total count of such sites is stored for each subtree.
    ##############################################################################
    # Call the Rcpp function to compute defining variants
    sub_trees_sv <- computeDefiningVariants(dna_aln, isolate_names, sub_trees)

    # Log completion of phase to standard output
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
    patient_seq_map <- split(names(seq2pt), seq2pt)
    subtree_metrics <- lapply(seq_along(sub_trees), function(st_i) {
        st <- sub_trees[[st_i]]
        tip_labels <- st$tip.label
        # Identify intake-positive sequences in the subtree
        ip_in_sub <- intersect(ip_seqs, tip_labels)
        # Determine the earliest date among convert sequences (non-intake-positives)
        convert_dates <- dates[setdiff(tip_labels, ip_pt_seqs)]
        min_convert_date <- if (length(convert_dates) > 0) min(convert_dates) else Inf
        # Filter intake-positives to those that occur on/before the cutoff
        valid_ip <- ip_in_sub[dates[ip_in_sub] <= min_convert_date]
        # Select one representative sequence per unique patient among valid intake-positives
        ip_patients <- unique(seq2pt[valid_ip])
        ip_rep <- vapply(ip_patients, function(pt) {
            seqs <- intersect(patient_seq_map[[pt]], tip_labels)
            other <- setdiff(tip_labels, seqs)
            # If there are no other sequences, return the first one
            # Otherwise, find the sequence with the minimum shared variant distance to all other sequences
            if (length(other) > 0) {
                min_vals <- sapply(seqs, function(seq) min(dna_shared_mat[other, seq]))
                seqs[which.min(min_vals)]
            } else {
                seqs[1]
            }
        }, character(1))
        # Get the distance from the representative to all other isolates in the subtree
        shared_counts_rep <- vapply(ip_rep, function(ip) {
            min(dna_shared_mat[setdiff(tip_labels, ip), ip])
        }, numeric(1))
        # Compute cluster score if there are defining variants and index isolates are not overly
        # distant from the representative sequence (i.e. their shared variant counts are nearly identical)
        score <- if (sub_trees_sv[st_i] > 0 && length(unique(shared_counts_rep)) <= 1) {
            # Reward intake-positives that could have started a cluster
            ip_convert_seq_count <- length(setdiff(
                tip_labels,
                ip_pt_seqs[seq2pt[ip_pt_seqs] %in% setdiff(unique(seq2pt), seq2pt[unlist(ip_rep)])]
            ))
            # Penalize intake-positives that could not have started a cluster
            ip_pt_count <- sum(tip_labels %in% ip_pt_seqs) / 1e6
            ip_convert_seq_count - ip_pt_count
        } else {
            0
        }
        # Compute index count: subtrees that have an index before all converts
        num_converts <- length(setdiff(tip_labels, ip_pt_seqs))
        num_ip_pts <- length(unique(seq2pt[intersect(tip_labels, ip_seqs)]))
        index_count <- if (length(valid_ip) > 0 && num_converts > 0) num_ip_pts else 0
        # Return both metrics
        list(score = score, index_count = index_count)
    })

    # Extract both metrics in a single step using a matrix conversion, and unlist
    metrics_mat <- do.call(rbind, subtree_metrics)
    sub_trees_valid <- unlist(metrics_mat[, "score"])
    sub_trees_index_first <- unlist(metrics_mat[, "index_count"])

    # Log completion of phase to standard output
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
        valid_indices <- which(sub_trees_valid > 0)
        st_scores <- vapply(valid_indices, function(st_i) {
            st <- sub_trees[[st_i]]
            # Skip and early return if the isolate is not in the subtree
            if (!(x %in% st$tip.label)) return(0)
            # Check for nested clusters
            nested <- vapply(setdiff(valid_indices, st_i), function(sub_i) {
                as.integer(all(sub_trees[[sub_i]]$tip.label %in% st$tip.label) &&
                               sub_trees_index_first[sub_i] < sub_trees_index_first[st_i] &&
                               sub_trees_index_first[sub_i] > 0)
            }, integer(1))
            # If not nested, return score
            if (sum(nested) == 0) sub_trees_valid[[st_i]] else 0
        }, numeric(1))
        # If no valid subtrees are found, default to cluster 1
        # Otherwise, assign the cluster based on the subtree with the highest score
        if (length(st_scores) == 0 || max(st_scores) <= 0) 1 else valid_indices[which.max(st_scores)]
    }, numeric(1))

    # Log completion of phase to standard output
    message("Phase 5 complete: Cluster assignment")

    # Return the clusters
    clusters
}
