#' Maps clusters onto a phylogenetic tree and visualizes them
#'
#' @author [Aryan Singh](mailto:aryansin@umich.edu) -
#'         ORCID ID: [0000-0002-1850-5598](https://orcid.org/0000-0002-1850-5598)
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint
#' @importFrom hues iwanthue
#' @importFrom dplyr .data left_join
#' @importFrom ggplot2 aes scale_color_manual ggtitle theme element_text unit guides guide_legend
#' @importFrom stats setNames
#' @export
plot_clusters_phylo <- function(tree, clusters) {
    # Convert phylo object to ggtree object
    tree <- ggtree(tree)

    # Format the clusters into a dataframe
    cluster_df <- data.frame(isolate = names(clusters), clust_id = factor(clusters))

    # Use iwanthue to generate colors for clusters
    cluster_colors <- setNames(
        iwanthue(length(unique(cluster_df$clust_id))),
        levels(cluster_df$clust_id)
    )

    # Add cluster information to the tree
    tree$data <- tree$data |> left_join(cluster_df, by = c("label" = "isolate"))

    # Return the tree labelled with clusters
    tree +
        geom_tippoint(aes(color = .data$clust_id), size = 1.5) +
        scale_color_manual(values = cluster_colors) +
        guides(
            color = guide_legend(
                title = "Cluster ID",
                override.aes = list(size = 5, shape = "square")
            )
        ) +
        theme(
            plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm")
        )
}

#' Compares the content of clusters created with two different methods.
#'
#' @param clusters1 A vector named by sequence IDs with values being subtrees defining the cluster.
#' @param clusters2 A vector named by the same sequence IDs as clusters1 with values being subtrees
#' defining the cluster.
#' @param width The width of the heatmap plot.
#' @param height The height of the heatmap plot.
#'
#' @return A ggplot plot object showing the overlap between clusters.
#'
#' @importFrom ggalign ggheatmap
#' @importFrom ggplot2 theme element_text scale_fill_gradient
#' @export
compare_clusters <- function(clusters1, clusters2, width = 10, height = 10) {
    # Unique cluster labels excluding 1 ("default" cluster)
    clusters1_unique <- setdiff(sort(unique(clusters1)), 1)
    clusters2_unique <- setdiff(sort(unique(clusters2)), 1)

    # Create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        nrow = length(clusters1_unique),
        ncol = length(clusters2_unique),
        dimnames = list(clusters1_unique, clusters2_unique)
    )

    clusters1_char <- as.character(clusters1_unique)
    clusters2_char <- as.character(clusters2_unique)

    # Compute cluster overlap
    for (cluster1 in clusters1_char) {
        for (cluster2 in clusters2_char) {
            cluster_overlap[cluster1, cluster2] <- length(intersect(
                names(clusters1)[clusters1 == cluster1],
                names(clusters2)[clusters2 == cluster2]
            )) /
                length(union(
                    names(clusters1)[clusters1 == cluster1],
                    names(clusters2)[clusters2 == cluster2]
                ))
        }
    }

    # Return the ggheatmap plot
    ggheatmap(cluster_overlap, width = width, height = height)
}

#' Produces summary plots regarding genetic distances within and between transmission clusters.
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_seqs A vector of sequence IDs corresponding to intake positive patients.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param prefix A descriptor indicating how the clusters were generated (used to name figure outputs).
#'
#' @return A matrix containing intra-cluster and inter-cluster genetic distances.
#'
#' @export
cluster_genetic_context <- function(clusters, seq2pt, ip_seqs, snp_dist, prefix) {
    # Subset of isolates assigned to clusters
    clusters_subset <- clusters[clusters != 1]
    unique_clusters <- sort(unique(clusters_subset))

    # Compute maximum genetic distance within cluster, minimum genetic distance to another cluster,
    # and minimum genetic distance to an isolate not in the same cluster
    cluster_distances <- vapply(
        unique_clusters,
        function(c) {
            cluster_seqs <- names(clusters[clusters == c])
            non_cluster_seqs <- names(clusters[clusters != c])

            max_intra <- max(snp_dist[cluster_seqs, cluster_seqs], na.rm = TRUE)
            min_inter_cluster <- min(
                snp_dist[cluster_seqs, names(clusters_subset[clusters_subset != c])],
                na.rm = TRUE
            )
            min_inter_isolate <- min(snp_dist[cluster_seqs, non_cluster_seqs], na.rm = TRUE)

            c(max_intra, min_inter_cluster, min_inter_isolate)
        },
        numeric(3)
    )

    # Define function for generating plots
    generate_plot <- function(x, y, xlab, ylab, file_name) {
        # Create directory if it doesn't exist
        dir.create("figures", showWarnings = FALSE)
        file_path <- paste0(
            "figures/",
            format(Sys.time(), "%Y-%m-%d"),
            "_",
            prefix,
            "_",
            file_name,
            ".pdf"
        )
        pdf(file_path)
        plot(
            jitter(x),
            jitter(y),
            xlab = xlab,
            ylab = ylab,
            xlim = c(0, max(x, na.rm = TRUE) + 1),
            ylim = c(0, max(y, na.rm = TRUE) + 1)
        )
        par(new = TRUE)
        upper_right <- c(0, min(max(x, na.rm = TRUE), max(y, na.rm = TRUE)))
        plot(
            upper_right,
            upper_right,
            xlab = "",
            ylab = "",
            type = "l",
            xlim = c(0, max(x, na.rm = TRUE) + 1),
            ylim = c(0, max(y, na.rm = TRUE) + 1)
        )
        dev.off()
    }

    # Generate plots
    generate_plot(
        cluster_distances[1, ],
        cluster_distances[2, ],
        "Max genetic distance within cluster",
        "Min genetic distance to another cluster",
        "intra_vs_inter_cluster_gen_dist"
    )
    generate_plot(
        cluster_distances[1, ],
        cluster_distances[3, ],
        "Max genetic distance within cluster",
        "Min genetic distance to another isolate",
        "intra_vs_inter_isolate_gen_dist"
    )

    # Create and return result matrix
    intra_inter_cluster_dist_mat <- t(cluster_distances)
    rownames(intra_inter_cluster_dist_mat) <- unique_clusters
    colnames(intra_inter_cluster_dist_mat) <- c(
        "Max intra-clust",
        "Min inter-clust",
        "Min inter-isolate"
    )

    intra_inter_cluster_dist_mat
}
