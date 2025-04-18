#' Maps clusters onto a phylogenetic tree and visualizes them
#'
#' @author [Aryan Singh](mailto:aryansin@umich.edu) -
#'         ORCID ID: [0000-0002-1850-5598](https://orcid.org/0000-0002-1850-5598)
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param tree A phylogenetic tree object of class `phylo`.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint
#' @importFrom hues iwanthue
#' @importFrom dplyr .data left_join
#' @importFrom ggplot2 aes scale_color_manual theme element_text unit guides guide_legend
#' @importFrom stats setNames
#' @export
plot_clusters_phylo <- function(clusters, tree) {
    # Convert phylo object to ggtree object
    tree <- ggtree(tree)

    # Format the clusters into a dataframe
    cluster_df <- data.frame(isolate = names(clusters), clust_id = factor(clusters))

    # Use iwanthue to generate colors for clusters
    cluster_colors <- setNames(iwanthue(length(unique(cluster_df$clust_id))), levels(cluster_df$clust_id))

    # Add cluster information to the tree
    tree$data <- tree$data |> left_join(cluster_df, by = c("label" = "isolate"))

    # Return the tree labelled with clusters
    tree +
        geom_tippoint(aes(color = .data$clust_id), size = 1.5) +
        scale_color_manual(values = cluster_colors) +
        guides(color = guide_legend(title = "Cluster ID", override.aes = list(size = 5, shape = "square"))) +
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
#' @importFrom ggalign ggheatmap
#' @importFrom ggplot2 theme element_text scale_fill_gradient
#'
#' @export
compare_clusters <- function(clusters1, clusters2, width = 10, height = 10) {
    # Unique cluster labels excluding 1 ("default" cluster)
    clusters1_unique <- setdiff(sort(unique(clusters1)), 1)
    clusters2_unique <- setdiff(sort(unique(clusters2)), 1)

    # Create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        nrow = length(clusters1_unique), ncol = length(clusters2_unique),
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
            )) / length(union(
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

#' @importFrom ape keep.tip
#' @importFrom ggtree gheatmap vexpand
#' @importFrom dplyr group_by slice_min ungroup rename left_join select mutate arrange mutate_all
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales colour_ramp
#' @importFrom ggplot2 geom_tile scale_fill_manual scale_fill_gradientn ggtitle
plot_trace_map <- function(sequence_ids, ip_seqs, seq2pt, snp_dist, dates, tree,
                           trace_data, clusters, # required for floor trace plotting
                           prefix = "plot", cluster_class = NULL, # optional
                           plot_type = c("pt_trace", "floor_trace")) {
    # Check that we have a valid plot type
    plot_type <- match.arg(plot_type)

    # Build a data frame for all sequences
    df <- data.frame(
        id = sequence_ids,
        pt = seq2pt[sequence_ids],
        date = dates[sequence_ids]
    )

    if (plot_type == "pt_trace") {
        # Identify which sequences/patients are present in pt_trace's columns
        df$in_trace <- df$pt %in% colnames(trace_data)
        df_sub <- df[df$in_trace, ] # Subset df to keep only those in pt_trace

        # Build the trace_sub matrix by transposing only the necessary columns
        trace_sub <- t(trace_data[, as.character(df_sub$pt), drop = FALSE])
        rownames(trace_sub) <- df_sub$id

        # Mark each isolate's date with a distinct value (2.5) in trace_sub
        # cbind() pairs the row name (isolate ID) and the date as a character
        trace_sub[cbind(df_sub$id, as.character(df_sub$date))] <- 2.5

        # If tr provided, keep only the relevant tips; else, nj tree
        if (!is.null(tree)) {
            tree_sub <- keep.tip(tree, df_sub$id)
            tree_plot <- ggtree(tree_sub)
        } else {
            tree_sub <- nj(snp_dist[df_sub$id, df_sub$id])
            tree_plot <- ggtree(tree_sub, size = 0)
        }

        # if no branch lengths create default lengths
        ## if (is.null(tree_sub$edge.length)) {
        ##  tree_sub$edge.length <- rep(1, nrow(tree_sub$edge))
        ## }

        # Column labeling: only label every 14th column i.e. every two weeks
        col_lab <- rownames(trace_data)
        idx_labeled <- seq(1, length(col_lab), 14)
        col_lab[setdiff(seq_along(col_lab), idx_labeled)] <- ""

        # Row labels: ID and patient ID in parentheses
        labels <- paste0(df_sub$id, "(", df_sub$pt, ")")

        # Row coloring: white for index patient isolates, lightblue otherwise
        row_colors <- ifelse(df_sub$id %in% ip_seqs, "white", "lightblue")

        # Create plot: tr + heatmap of trace + annote rows w the patient/origin color
        gheatmap(tree_plot, as.data.frame(trace_sub) |> mutate_all(as.factor),
            offset = 1, width = 4, font.size = 2,
            custom_column_labels = col_lab, colnames_angle = 90
        ) + vexpand(.1, -1) + # adjust for long vertical column labels
            scale_fill_manual(values = c("white", "lightgray", "blue", "red", "yellow"), name = "Trace") +
            new_scale_fill() + # color the tiles
            geom_fruit(
                data = data.frame(id = df_sub$id, color = row_colors),
                geom = geom_tile,
                mapping = aes(y = id, x = 1, fill = color),
                offset = -0.2, width = 1
            ) +
            scale_fill_manual(
                name = "Isolate Origin",
                values = c("white" = "white", "lightblue" = "lightblue"),
                labels = c("Index Patient", "Non-Index Patient")
            ) + ggtitle(paste0(prefix, " - Patient Trace Map"))
    } else if (plot_type == "floor_trace") {
        if (is.null(clusters)) {
            stop("Clusters vector must be provided for floor_trace plot type")
        }

        # Update df to include cluster info
        df$cluster <- as.character(clusters)

        # Ensure cluster_class has correct names if provided; otherwise create a simple vector
        if (is.null(cluster_class)) {
            cluster_order <- names(table(clusters))
            cluster_class <- structure(rep("NA", length(unique(clusters))), names = unique(clusters))
        } else {
            # Only keep cluster names that actually appear in 'clusters'
            cluster_order <- intersect(names(cluster_class), unique(df$cluster))
        }

        # For each cluster-pt group, find earliest isolate (in time)
        # Then map each isolate to that earliest isolate
        df_earliest <- df |>
            group_by(cluster, pt) |>
            slice_min(date, with_ties = FALSE) |>
            ungroup() |>
            rename(earliest_id = id)

        # Add earliest isolate info back to df
        df <- left_join(df, select(df_earliest, cluster, pt, earliest_id), by = c("cluster", "pt"))

        # Identify isolate order: cluster_order, then by earliest isolate's date
        df_earliest <- df_earliest |>
            mutate(cluster_factor = factor(cluster, levels = cluster_order)) |>
            arrange(cluster_factor, date)
        isolate_order <- df_earliest$earliest_id

        # Identify if each earliest isolate is in the user-supplied index set
        df_earliest <- df_earliest |> mutate(is_index = earliest_id %in% ip_seqs)

        # Create a named vector to map each isolate to the earliest isolate of its cluster–patient group.
        rep_isolate <- setNames(df$earliest_id, df$id)

        # Helper function to get max SNP distance within cluster or within cluster–patient
        get_max_snp <- function(iso_id, use_entire_cluster = FALSE) {
            clust <- df$cluster[df$id == iso_id]
            ptval <- df$pt[df$id == iso_id]
            sub_ids <- if (!use_entire_cluster) {
                df$id[df$cluster == clust & df$pt == ptval] # same cluster, same patient
            } else {
                df$id[df$cluster == clust] # entire cluster
            }
            max(snp_dist[sub_ids, sub_ids])
        }

        # For each earliest isolate, compute max distances
        max_dist <- sapply(isolate_order, get_max_snp, use_entire_cluster = FALSE)
        max_clust_dist <- sapply(isolate_order, get_max_snp, use_entire_cluster = TRUE)

        # Subset the floor_trace matrix
        pt_in_trace <- seq2pt[sequence_ids] %in% colnames(trace_data)
        seq_ids_in_trace <- sequence_ids[pt_in_trace]
        trace_sub <- t(trace_data[, as.character(seq2pt[seq_ids_in_trace]), drop = FALSE])
        rownames(trace_sub) <- seq_ids_in_trace

        # Mark the actual isolate day with 1.75
        these_dates <- dates[seq_ids_in_trace]
        trace_sub[cbind(rep_isolate[seq_ids_in_trace], as.character(these_dates))] <- 1.75

        # Reorder rows by earliest isolates in cluster/date order
        keep_rows <- intersect(isolate_order, rownames(trace_sub))
        trace_sub <- trace_sub[keep_rows, , drop = FALSE]

        # Identify convert isolates
        is_convert_isolate <- sapply(keep_rows, function(id) {
            row_vals <- trace_sub[id, ]
            i_exposed <- which(row_vals %in% c(1.25))
            i_positive <- which(row_vals %in% c(1.5, 1.75))
            if (length(i_exposed) == 0 || length(i_positive) == 0) {
                return(FALSE)
            }
            is_convert <- min(i_exposed) < min(i_positive)
            trace_date <- as.numeric(colnames(trace_sub)[min(i_positive)])
            iso_date <- dates[as.character(id)]
            (iso_date - trace_date < 7) && is_convert && !(id %in% ip_seqs)
        })

        # Index isolates
        is_index_isolate <- keep_rows %in% ip_seqs

        # Build annotation data
        patient_labels <- c("Secondary convert", "Index patient", "Convert patient")
        convert_class <- rep("Secondary convert", length(keep_rows))
        convert_class[is_index_isolate] <- "Index patient"
        convert_class[!is_index_isolate & is_convert_isolate] <- "Convert patient"

        # Construct the annotation data frame.
        annotation_row <- data.frame(
            Convert          = convert_class,
            Intra_pt_dist    = max_dist[keep_rows] + 1,
            Intra_clust_dist = max_clust_dist[keep_rows] + 1,
            row.names        = keep_rows,
            stringsAsFactors = FALSE
        )

        # Annot colors
        bluescale <- colorRampPalette(c("white", "blue"))
        greenscale <- colorRampPalette(c("white", "darkgreen"))

        ann_colors <- list(
            Convert          = structure(c("gray95", "gray", "red"), names = patient_labels),
            Intra_pt_dist    = bluescale(max(max_dist, 1)),
            Intra_clust_dist = greenscale(max(max_clust_dist, 1))
        )

        annotation_row$id <- rownames(annotation_row)

        # Gaps between clusters
        ctab <- table(annotation_row$Cluster)
        ctab <- ctab[as.character(cluster_order)]
        row_gaps <- cumsum(ctab)
        row_gaps <- row_gaps[row_gaps < nrow(annotation_row)]

        # Dynamic font sizing based on the dimensions of the matrix
        fsize_row <- 2 * (1 - (nrow(trace_sub) / ncol(trace_data)))

        # Prepare the column labels
        col_lab <- rownames(trace_data)
        idx_keep_labels <- seq(1, length(col_lab), 14)
        col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

        labels <- paste0(keep_rows, "(", seq2pt[keep_rows], ")") # row labels

        # Set up breaks and colors
        custom_breaks <- c(0, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6)
        custom_colors <- brewer.pal(7, "Pastel2")
        custom_colors <- c("white", custom_colors[1], "blue", "red", "yellow", custom_colors[c(2:5, 7)])

        tree_sub <- keep.tip(tree, keep_rows)
        trace_df <- as.data.frame(trace_sub)
        trace_df[] <- lapply(trace_df, function(x) {
            factor(as.numeric(as.character(x)), levels = custom_breaks)
        })

        tree_plot <- ggtree(tree_sub, branch.length = "none")

        gheatmap(tree_plot, trace_df, offset = 4, width = 4, font.size = 3,
            custom_column_labels = col_lab, colnames_angle = 90, colnames_offset_y = -8
        ) +
            vexpand(0.1, -1) +
            scale_fill_manual(
                name = "Trace", values = custom_colors, breaks = custom_breaks, labels = custom_breaks,
                na.value = "white", drop = FALSE
            ) + new_scale_fill() +
            # Convert annot
            geom_fruit(
                data = annotation_row, geom = geom_tile, mapping = aes(y = id, x = 1, fill = Convert),
                offset = -0.1825, width = 1
            ) +
            scale_fill_manual(
                name = "Convert", values = ann_colors$Convert, labels = names(ann_colors$Convert), drop = FALSE
            ) + new_scale_fill() +
            # Intra_clust_dist annot
            geom_fruit(
                data = annotation_row, geom = geom_tile, mapping = aes(y = id, x = 1, fill = Intra_clust_dist),
                offset = -0.1725, width = 1
            ) +
            scale_fill_gradientn(name = "Intra_clust_dist", colors = ann_colors$Intra_clust_dist) + new_scale_fill() +
            # Intra_pt_dist annot
            geom_fruit(
                data = annotation_row, geom = geom_tile, mapping = aes(y = id, x = 1, fill = Intra_pt_dist),
                offset = -0.1625, width = 1
            ) +
            scale_fill_gradientn(
                name = "Intra_pt_dist", colors = ann_colors$Intra_pt_dist
            ) + ggtitle(paste0(prefix, " - Floor Trace Map"))
    }
}

