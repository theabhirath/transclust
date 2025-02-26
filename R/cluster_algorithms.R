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
