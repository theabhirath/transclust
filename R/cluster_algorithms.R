#' Perform clustering of isolates using a hard SNP distance cutoff.
#'
#' @description
#' This function uses a hard SNP distance cutoff to define clusters naïvely, without using
#' the phylogenetic tree.
#'
#' @param snp_dist A matrix of SNP distances between isolates constructed using a model of DNA evolution.
#'                 See [`get_snp_dist_matrix`] for a useful function to generate this.
#' @param snp_thresh A threshold for defining clusters.
#'
#' @return A numeric vector indicating the cluster that each isolate belongs to.
#'
#' @importFrom stats hclust as.dist cutree
#' @export
get_tn_clusters_snp_thresh <- function(snp_dist, snp_thresh) {
    # convert distance matrix to hclust object
    snp_hclust <- hclust(as.dist(snp_dist))
    # cut at the specified threshold to create clusters and return the clusters
    cutree(snp_hclust, h = snp_thresh)
}

#' Identify transmission clusters based on the number of shared variants.
#'
#' @description
#' Clustering is performed to identify the maximal clusters containing a single intake-positive patient that occurs
#' before all cluster converts. The clustering metric is the number of shared variants, and clusters can have multiple
#' intake-positive patients if they share an identical number of variants with other cluster members or intake-positive
#' patients occur after converts. This clustering also requires that clusters be defined by at least one shared variant
#' that other isolates don't have.
#'
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param snp_dist A matrix of SNP distances between isolates constructed using a model of DNA evolution.
#'                 See [`get_snp_dist_matrix`] for a useful function to generate this.
#' @param ip_seqs A vector of sequence IDs which correspond to intake patient sequences presumed to be imported.
#' @param ip_pt_seqs A vector of all sequence IDs which correspond to intake-positive patients, imported or
#'                   collected later. This will be a superset of `ip_seqs` by definition.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param dates A vector of isolate dates named by sequence IDs.
#' @param tree A phylogenetic tree object of class `phylo` constructed from the DNA alignment. This can be constructed
#'             using the [`get_phylo_tree`] or can be any other tree object constructed from the same isolates.
#'
#' @return A numeric vector indicating the cluster that each isolate belongs to.
#'
#' @references Hawken, S. E., Yelin, R. D., Lolans, K., Pirani, A., Weinstein, R. A., Lin, M. Y., Hayden, M. K., &
#'             Snitkin, E. S. (2022). Threshold-free genomic cluster detection to track transmission pathways in
#'             health-care settings: A genomic epidemiology analysis. The Lancet Microbe, 3(9), e652–e662.
#'             \doi{10.1016/S2666-5247(22)00115-X}
#'
#' @importFrom ape subtrees
#' @export
get_tn_clusters_sv_index <- function(dna_aln, snp_dist, ip_seqs, ip_pt_seqs, seq2pt, dates, tree) {
    #####################################################################################
    # 1. Compute the shared variant matrix #####
    # For each pair of isolates, we compute the number of positions where:
    #   a. the isolate's base differs from the out-group,
    #   b. the two isolates share the same base,
    #   c. the base is not a gap ('-') or ambiguous ('n'),
    # and then we convert this similarity count to a distance-like measure.
    ####################################################################################
    # get root of the tree, which is the out-group
    parents <- unique(tree$edge[, 1])
    children <- unique(tree$edge[, 2])
    out_group <- setdiff(parents, children)
    # convert the DNA alignment object to a character matrix
    dna_char <- as.character(dna_aln)
    # call the Rcpp function to compute the shared variant matrix
    dna_shared_mat <- computeSharedMatrix(dna_char, out_group)

    isolate_names <- row.names(dna_aln)
    row.names(dna_shared_mat) <- isolate_names
    colnames(dna_shared_mat) <- isolate_names

    # Log completion of phase to standard output
    message("Phase 1 complete: Shared variant matrix computation.")

    ##############################################################################
    # 2. Identify defining variants in each subtree #####
    # For every subtree, we loop over each alignment column (site) to check:
    #   - All tips in the subtree have the same base (and that base is not 'n').
    #   - Outside the subtree, there is at least one different base.
    #   - No outside isolate has the same base as the subtree.
    # The total count of such sites is stored for each subtree.
    ##############################################################################
    # get subtrees from the provided phylogenetic tree
    sub_trees <- subtrees(tree)
    # call the Rcpp function to compute defining variants
    sub_trees_dv <- computeDefiningVariants(dna_char, isolate_names, sub_trees)

    # log completion of phase to standard output
    message("Phase 2 complete: Defining variant identification.")

    #####################################################################################
    # 3. Validate and score subtrees as clusters #####
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
        # identify intake-positive sequences in the subtree
        ip_in_sub <- intersect(ip_seqs, tip_labels)
        # determine the earliest date among convert sequences (non-intake-positive patients)
        convert_dates <- dates[setdiff(tip_labels, ip_pt_seqs)]
        min_convert_date <- if (length(convert_dates) > 0) min(convert_dates) else Inf
        # filter intake-positives to those that occur on/before the cutoff
        valid_ip <- ip_in_sub[dates[ip_in_sub] <= min_convert_date]
        # select one representative sequence per unique patient among valid intake-positives
        ip_patients <- unique(seq2pt[valid_ip])
        ip_rep <- vapply(
            ip_patients,
            function(pt) {
                seqs <- intersect(patient_seq_map[[as.character(pt)]], tip_labels)
                other <- setdiff(tip_labels, seqs)
                # if there are no other sequences, return the first one
                # otherwise, find the sequence with the minimum shared variant distance to all other sequences
                if (length(other) > 0) {
                    min_vals <- vapply(
                        seqs,
                        function(seq) min(dna_shared_mat[other, seq]),
                        numeric(1)
                    )
                    seqs[which.min(min_vals)]
                } else {
                    seqs[1]
                }
            },
            character(1)
        )
        # get the distance from the representative to all other isolates in the subtree
        shared_counts_rep <- vapply(
            ip_rep,
            function(ip) {
                min(dna_shared_mat[setdiff(tip_labels, ip), ip])
            },
            numeric(1)
        )
        # compute cluster score if there are defining variants and index isolates are not overly
        # distant from the representative sequence (i.e. their shared variant counts are nearly identical)
        score <- if (sub_trees_dv[st_i] > 0 && length(unique(shared_counts_rep)) <= 1) {
            # all patients that are not represented in this subtree by an intake-positive
            patients_non_rep <- setdiff(unique(seq2pt), seq2pt[unlist(ip_rep)])
            # reward intake-positives that could have started a cluster
            ip_convert_seq_count <- length(setdiff(
                tip_labels,
                ip_pt_seqs[seq2pt[ip_pt_seqs] %in% patients_non_rep]
            ))
            # penalize intake-positives that could not have started a cluster
            ip_pt_count <- sum(tip_labels %in% ip_pt_seqs) / 1e6
            ip_convert_seq_count - ip_pt_count
        } else {
            0
        }
        # compute index count: subtrees that have an index before all converts
        num_converts <- length(setdiff(tip_labels, ip_pt_seqs))
        num_ip_pts <- length(unique(seq2pt[intersect(tip_labels, ip_seqs)]))
        index_count <- if (length(valid_ip) > 0 && num_converts > 0) num_ip_pts else 0
        # return both metrics
        list(score = score, index_count = index_count)
    })

    # extract both metrics in a single step using a matrix conversion, and unlist
    metrics_mat <- do.call(rbind, subtree_metrics)
    sub_trees_valid <- unlist(metrics_mat[, "score"])
    sub_trees_index_first <- unlist(metrics_mat[, "index_count"])

    # log completion of phase to standard output
    message("Phase 3 complete: Cluster validation and scoring.")

    #####################################################################################
    # 4. Assign clusters to each isolate #####
    # For every isolate:
    #  - Check all valid subtrees (i.e. those with positive scores) that contain it.
    #  - Ensure that a candidate subtree is not "nested" within another subtree that
    #    has an earlier index (i.e. a lower index_count).
    #  - Assign the cluster based on the subtree with the highest score.
    #  - If no valid subtree is found, default to cluster 1.
    #  - Remap cluster values: non-1 values are mapped to 1:n, each 1 is treated as separate.
    #####################################################################################
    clusters <- vapply(
        isolate_names,
        function(x) {
            valid_indices <- which(sub_trees_valid > 0)
            st_scores <- vapply(
                valid_indices,
                function(st_i) {
                    st <- sub_trees[[st_i]]
                    # skip and early return if the isolate is not in the subtree
                    if (!(x %in% st$tip.label)) {
                        return(0)
                    }
                    # check for nested clusters
                    nested <- vapply(
                        setdiff(valid_indices, st_i),
                        function(sub_i) {
                            as.integer(
                                all(sub_trees[[sub_i]]$tip.label %in% st$tip.label) &&
                                    sub_trees_index_first[sub_i] < sub_trees_index_first[st_i] &&
                                    sub_trees_index_first[sub_i] > 0
                            )
                        },
                        integer(1)
                    )
                    # if not nested, return score
                    if (sum(nested) == 0) sub_trees_valid[[st_i]] else 0
                },
                numeric(1)
            )
            # if no valid subtrees are found, default to cluster 0
            # otherwise, assign the cluster based on the subtree with the highest score
            if (length(st_scores) == 0 || max(st_scores) <= 0) {
                0
            } else {
                valid_indices[which.max(st_scores)]
            }
        },
        numeric(1)
    )

    # log completion of phase to standard output
    message("Phase 4 complete: Cluster assignment.")

    # remap cluster values to ensure that clusters are numbered sequentially starting from 1
    remap_cluster_values(clusters)
}
