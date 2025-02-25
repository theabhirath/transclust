#' Compare Clusters
#'
#' Compares the content of clusters created with two different methods.
#'
#' @param clusters1 A vector named by sequence IDs with values being subtrees defining the cluster.
#' @param clusters2 A vector named by the same sequence IDs as clusters1 with values being subtrees
#' defining the cluster.
#' @param prefix Prefix to use for naming figure output files.
#'
#' @return A comparison of the content of clusters.
#' @export
#'
#' @examples
#' \dontrun{
#' compare_clusters(clusters1, clusters2, "prefix")
#' }
#'
compare_clusters <- function(clusters1, clusters2, prefix) {
    # Unique cluster labels
    clusters1_unique <- unique(clusters1)
    clusters2_unique <- unique(clusters2)

    # Create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        ncol = length(clusters2_unique), nrow = length(clusters1_unique),
        dimnames = list(sort(clusters1_unique), sort(clusters2_unique))
    )

    # Compute cluster overlap
    for (c1 in as.character(sort(clusters1_unique))) {
        for (c2 in as.character(sort(clusters2_unique))) {
            cluster_overlap[c1, c2] <- length(intersect(
                names(clusters1)[clusters1 == c1],
                names(clusters2)[clusters2 == c2]
            )) / length(union(
                names(clusters1)[clusters1 == c1],
                names(clusters2)[clusters2 == c2]
            ))
        }
    }

    # Plot heatmap of cluster overlap
    file <- paste("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_member_compare.pdf", sep = "")
    pheatmap::pheatmap(cluster_overlap, filename = file)
}

#' Cluster Genetic Context
#'
#' Produces summary plots regarding genetic distances within and between transmission clusters.
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_seqs A vector of sequence IDs corresponding to intake positive patients.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param prefix A descriptor indicating how the clusters were generated (used to name figure outputs).
#'
#' @return A matrix containing intra-cluster and inter-cluster genetic distances.
#' @export
#'
#' @examples
#' \dontrun{
#' # distances <- cluster_genetic_context(clusters, seq2pt, ip_seqs, snp_dist, "my_analysis")
#' }
#'
cluster_genetic_context <- function(clusters, seq2pt, ip_seqs, snp_dist, prefix) {
    # Subset of isolates assigned to clusters
    clusters_subset <- clusters[clusters != 1]
    unique_clusters <- sort(unique(clusters_subset))

    # Compute maximum genetic distance within cluster, minimum genetic distance to another cluster,
    # and minimum genetic distance to an isolate not in the same cluster
    cluster_distances <- vapply(unique_clusters, function(c) {
        cluster_seqs <- names(clusters[clusters == c])
        non_cluster_seqs <- names(clusters[clusters != c])

        max_intra <- max(snp_dist[cluster_seqs, cluster_seqs], na.rm = TRUE)
        min_inter_cluster <- min(snp_dist[cluster_seqs, names(clusters_subset[clusters_subset != c])], na.rm = TRUE)
        min_inter_isolate <- min(snp_dist[cluster_seqs, non_cluster_seqs], na.rm = TRUE)

        c(max_intra, min_inter_cluster, min_inter_isolate)
    }, numeric(3))

    # Define function for generating plots
    generate_plot <- function(x, y, xlab, ylab, file_name) {
        # Create directory if it doesn't exist
        dir.create("figures", showWarnings = FALSE)
        file_path <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_", file_name, ".pdf")
        pdf(file_path)
        plot(jitter(x), jitter(y),
            xlab = xlab, ylab = ylab,
            xlim = c(0, max(x, na.rm = TRUE) + 1),
            ylim = c(0, max(y, na.rm = TRUE) + 1)
        )
        par(new = TRUE)
        upper_right <- c(0, min(max(x, na.rm = TRUE), max(y, na.rm = TRUE)))
        plot(upper_right, upper_right,
            xlab = "", ylab = "", type = "l",
            xlim = c(0, max(x, na.rm = TRUE) + 1),
            ylim = c(0, max(y, na.rm = TRUE) + 1)
        )
        dev.off()
    }

    # Generate plots
    generate_plot(
        cluster_distances[1, ], cluster_distances[2, ],
        "Max genetic distance within cluster",
        "Min genetic distance to another cluster",
        "intra_vs_inter_cluster_gen_dist"
    )
    generate_plot(
        cluster_distances[1, ], cluster_distances[3, ],
        "Max genetic distance within cluster",
        "Min genetic distance to another isolate",
        "intra_vs_inter_isolate_gen_dist"
    )

    # Create and return result matrix
    intra_inter_cluster_dist_mat <- t(cluster_distances)
    rownames(intra_inter_cluster_dist_mat) <- unique_clusters
    colnames(intra_inter_cluster_dist_mat) <- c("Max intra-clust", "Min inter-clust", "Min inter-isolate")

    intra_inter_cluster_dist_mat
}
