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
#'
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
#'
#' @export
compare_clusters <- function(clusters1, clusters2, width = 10, height = 10) {
    # Unique cluster labels
    clusters1_unique_labels <- as.character(sort(unique(clusters1)))
    clusters2_unique_labels <- as.character(sort(unique(clusters2)))

    # Create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        nrow = length(clusters1_unique_labels),
        ncol = length(clusters2_unique_labels),
        dimnames = list(clusters1_unique_labels, clusters2_unique_labels)
    )

    # Compute cluster overlap
    for (cluster1 in clusters1_unique_labels) {
        for (cluster2 in clusters2_unique_labels) {
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

    # get clusters that map exactly 1-1
    one_to_one_mapping <- rowSums(cluster_overlap == 1) == 1
    # get isolates that are in one-to-one mapping
    isolates_in_one_to_one_mapping <- names(clusters1)[one_to_one_mapping]
    # percentage of isolates that are in one-to-one mapping
    overlap_percentage <- length(isolates_in_one_to_one_mapping) / length(clusters1)

    # Return the ggheatmap plot
    ggheatmap(cluster_overlap, width = width, height = height) +
        ggtitle(paste0(
            "Comparison of clusters. Isolate overlap: ",
            round(overlap_percentage * 100, 2),
            "%"
        ))
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

#' @importFrom ape keep.tip
#' @importFrom ggtree gheatmap vexpand
#' @importFrom dplyr group_by slice_min ungroup rename left_join select mutate arrange mutate_all
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom paletteer paletteer_d
#' @importFrom ggplot2 geom_tile scale_fill_manual scale_fill_gradientn ggtitle
#' @importFrom ggtree geom_tiplab
plot_trace <- function(tree, clusters, snp_dist, dates, ip_seqs, seq2pt, trace_data) {
    # prepare a data frame for all isolates
    df <- data.frame(
        id = names(clusters),
        cluster = as.character(clusters),
        patient = seq2pt[names(clusters)],
        date = dates[names(clusters)]
    )

    # get cluster order
    cluster_order <- names(table(clusters))

    # for each cluster-patient group, find earliest isolate in time
    # then map each isolate to that earliest isolate
    df_earliest <- df |>
        group_by(cluster, patient) |>
        slice_min(date, with_ties = FALSE) |>
        ungroup() |>
        rename(earliest_id = id)

    # add earliest isolate info back to df
    df <- left_join(
        df,
        select(df_earliest, cluster, patient, earliest_id),
        by = c("cluster", "patient")
    )

    # identify isolate order: cluster_order, then by earliest isolate's date
    df_earliest <- df_earliest |>
        mutate(cluster_factor = factor(cluster, levels = cluster_order)) |>
        arrange(cluster_factor, date)
    isolate_order <- df_earliest$earliest_id

    # identify if each earliest isolate is in the user-supplied index set
    df_earliest <- df_earliest |> mutate(is_index = earliest_id %in% ip_seqs)

    # create a named vector to map each isolate to the earliest isolate of its cluster–patient group
    rep_isolate <- setNames(df$earliest_id, df$id)

    # get max SNP distance for all isolates for a single patient or for an entire cluster
    get_max_snp <- function(iso_id, patient_only = TRUE) {
        clust <- df$cluster[df$id == iso_id]
        ptval <- df$patient[df$id == iso_id]
        sub_ids <- if (patient_only) {
            df$id[df$cluster == clust & df$patient == ptval] # same cluster, same patient
        } else {
            df$id[df$cluster == clust] # entire cluster
        }
        max(snp_dist[sub_ids, sub_ids])
    }

    # for each earliest isolate, compute max distances
    max_dist <- sapply(isolate_order, get_max_snp, patient_only = TRUE)
    max_clust_dist <- sapply(isolate_order, get_max_snp, patient_only = FALSE)

    # subset the trace_data matrix
    pt_in_trace <- df$patient %in% colnames(trace_data)
    seq_ids_in_trace <- df$id[pt_in_trace]
    trace_sub <- t(trace_data[, as.character(df$patient[pt_in_trace]), drop = FALSE])
    rownames(trace_sub) <- seq_ids_in_trace

    # mark the actual isolate day with 1.75
    trace_sub[cbind(rep_isolate[seq_ids_in_trace], as.character(df$date[pt_in_trace]))] <- 1.75

    # reorder rows by earliest isolates in cluster/date order
    keep_rows <- intersect(isolate_order, rownames(trace_sub))
    trace_sub <- trace_sub[keep_rows, , drop = FALSE]

    # identify convert isolates
    is_convert_isolate <- sapply(keep_rows, function(id) {
        row_vals <- trace_sub[id, ]
        i_exposed <- which(row_vals %in% c(1.25))
        i_positive <- which(row_vals %in% c(1.5, 1.75))
        if (length(i_exposed) == 0 || length(i_positive) == 0) {
            return(FALSE)
        }
        is_convert <- min(i_exposed) < min(i_positive)
        trace_date <- as.numeric(colnames(trace_sub)[min(i_positive)])
        iso_date <- df$date[df$id == id]
        (iso_date - trace_date < 7) && is_convert && !(id %in% ip_seqs)
    })

    # index isolates
    is_index_isolate <- keep_rows %in% ip_seqs

    # build annotation data
    patient_labels <- c("Secondary convert", "Index patient", "Convert patient")
    convert_class <- rep("Secondary convert", length(keep_rows))
    convert_class[is_index_isolate] <- "Index patient"
    convert_class[!is_index_isolate & is_convert_isolate] <- "Convert patient"

    # construct the annotation data frame
    annotation_row <- data.frame(
        id = keep_rows,
        Cluster = df$cluster[match(keep_rows, df$id)],
        Convert = convert_class,
        Intra_pt_dist = max_dist[keep_rows] + 1,
        Intra_clust_dist = max_clust_dist[keep_rows] + 1
    )

    # annotation colors
    cluster_colors <- setNames(iwanthue(length(unique(df$cluster))), levels(df$cluster))
    bluescale <- colorRampPalette(c("white", "blue"))
    greenscale <- colorRampPalette(c("white", "darkgreen"))
    ann_colors <- list(
        Cluster = structure(
            cluster_colors[seq_along(unique(df$cluster))],
            names = sort(unique(df$cluster))
        ),
        Convert = structure(c("gray95", "gray", "red"), names = patient_labels),
        Intra_pt_dist = bluescale(max(max_dist, 1)),
        Intra_clust_dist = greenscale(max(max_clust_dist, 1))
    )

    # prepare the column labels
    col_lab <- colnames(trace_sub)
    idx_keep_labels <- seq(1, length(col_lab), 14)
    col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

    # set up breaks and colors
    # get unique values from trace_sub and sort them
    custom_breaks <- sort(unique(as.numeric(as.matrix(trace_sub))))

    # calculate number of colors needed (excluding white for 0)
    n_colors <- length(custom_breaks) - 1

    # Generate colors. All colors must be distinct.
    # If n_colors > 9, recycle colors.
    base_colors <- paletteer_d("ggthemes::Classic_Cyclic")
    if (n_colors > 13) {
        base_colors <- rep(base_colors, length.out = n_colors)
        warning("More than 13 colors needed for trace. Recycled colors.")
    }
    # three high-contrast colors which clash with the color scheme, cannot be white as that is the base color
    high_contrast_colors <- c("#000000", "red", "darkgray")

    # create final color vector with white for 0 value
    custom_colors <- c(
        "white",
        base_colors[1],
        high_contrast_colors,
        base_colors[2:length(base_colors)]
    )

    tree_sub <- keep.tip(tree, keep_rows)
    # tree_sub$tip.label <- paste0(tree_sub$tip.label, "(", seq2pt[tree_sub$tip.label], ")")
    trace_df <- as.data.frame(trace_sub)
    trace_df[] <- lapply(trace_df, function(x) {
        factor(as.numeric(as.character(x)), levels = custom_breaks)
    })

    tree_plot <- ggtree(tree_sub, branch.length = "none") +
        # add labels for each tip
        geom_tiplab(size = 0.75, align = TRUE, offset = 3, linetype = NULL)

    gheatmap(
        tree_plot,
        trace_df,
        offset = 6,
        width = 4,
        font.size = 1,
        custom_column_labels = col_lab,
        colnames_angle = 90,
        colnames_offset_y = -8
    ) +
        vexpand(0.1, -1) +
        scale_fill_manual(
            name = "Trace",
            values = custom_colors,
            breaks = custom_breaks,
            labels = custom_breaks,
            na.value = "white",
            drop = FALSE
        ) +
        new_scale_fill() +
        # Cluster annot
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Cluster),
            offset = -0.2,
            width = 1
        ) +
        scale_fill_manual(
            name = "Cluster",
            values = ann_colors$Cluster,
            drop = FALSE,
            guide = "none"
        ) +
        new_scale_fill() +
        # Convert annot
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Convert),
            offset = -0.1825,
            width = 1
        ) +
        scale_fill_manual(
            name = "Convert",
            values = ann_colors$Convert,
            labels = names(ann_colors$Convert),
            drop = FALSE
        ) +
        new_scale_fill() +
        # Intra_clust_dist annot
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Intra_clust_dist),
            offset = -0.1725,
            width = 1
        ) +
        scale_fill_gradientn(name = "Intra_clust_dist", colors = ann_colors$Intra_clust_dist) +
        new_scale_fill() +
        # Intra_pt_dist annot
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Intra_pt_dist),
            offset = -0.1625,
            width = 1
        ) +
        scale_fill_gradientn(name = "Intra_pt_dist", colors = ann_colors$Intra_pt_dist)
}

#' Plots the patient and floor trace heatmaps for each cluster separately.
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param snp_dist A matrix of SNP distances between isolates
#' @param ip_seqs Vector of sequence IDs presumed to be imported
#' @param ip_pt_seqs Vector of sequence IDs for intake positive patients
#' @param seq2pt A named vector mapping sequence IDs to patient IDs
#' @param dates A named vector of isolate dates by sequence ID
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param prefix A string descriptor used for naming figure outputs.
#' @param pt_trace A data frame or matrix representing daily patient trace (rows) by patient ID (columns)
#' @param floor_trace Optional data for floor-level tracing.
#'
#' @importFrom ggplot2 ggsave
#' @export
#' @keywords internal
cluster_context <- function(
    tree,
    clusters,
    snp_dist,
    ip_seqs,
    ip_pt_seqs,
    seq2pt,
    dates,
    pt_trace,
    floor_trace = NULL,
    by = c("isolate", "patient"),
    prefix = "plot"
) {
    by <- match.arg(by)
    # remove singleton clusters
    clusters <- remove_singleton_clusters(clusters)
    # Iterate over each cluster
    for (cluster in setdiff(unique(clusters), 1)) {
        # Get patients corresponding to cluster members
        cluster_members <- names(clusters)[clusters == cluster]
        cluster_pts <- seq2pt[as.character(cluster_members)]
        # Plot trace
        ggsave(
            paste0("figures/", prefix, "_pt_trace_", cluster, ".pdf"),
            plot_trace(
                tree,
                clusters[clusters == cluster],
                snp_dist,
                dates,
                ip_seqs,
                seq2pt,
                pt_trace
            )
        )
        # Plot floor trace if given
        if (!is.null(floor_trace) && length(unique(cluster_pts)) > 1) {
            ggsave(
                paste0("figures/", prefix, "_floor_trace_", cluster, ".pdf"),
                plot_trace(
                    tree,
                    clusters[clusters == cluster],
                    snp_dist,
                    dates,
                    ip_seqs,
                    seq2pt,
                    floor_trace
                )
            )
        }
    }
}

#' Plots the patient and floor trace heatmaps for all clusters belonging to a single sequence type.
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param snp_dist A matrix of SNP distances between isolates
#' @param ip_seqs Vector of sequence IDs presumed to be imported
#' @param ip_pt_seqs Vector of sequence IDs for intake positive patients
#' @param seq2pt A named vector mapping sequence IDs to patient IDs
#' @param dates A named vector of isolate dates by sequence ID
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param pt_trace A data frame or matrix representing daily patient trace (rows) by patient ID (columns)
#' @param floor_trace Optional data for floor-level tracing.
#' @param prefix A string descriptor used for naming figure outputs.
#'
#' @importFrom ggplot2 ggsave
#' @export
st_context <- function(
    tree,
    clusters,
    snp_dist,
    ip_seqs,
    seq2pt,
    dates,
    pt_trace,
    floor_trace = NULL,
    prefix = "plot"
) {
    # TODO: handle this so that only one patient is plotted per cluster
    # plot trace for entire cluster
    ggsave(
        paste0("figures/", prefix, "_pt_trace_", "all", ".pdf"),
        plot_trace(tree, clusters, snp_dist, dates, ip_seqs, seq2pt, pt_trace)
    )
    # plot floor trace if given
    if (!is.null(floor_trace)) {
        ggsave(
            paste0("figures/", prefix, "_floor_trace_", "all", ".pdf"),
            plot_trace(tree, clusters, snp_dist, dates, ip_seqs, seq2pt, floor_trace)
        )
    }
}
