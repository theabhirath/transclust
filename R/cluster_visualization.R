#' Maps clusters onto a phylogenetic tree and visualizes them
#'
#' @author [Aryan Singh](mailto:aryansin@umich.edu) -
#'         ORCID ID: [0000-0002-1850-5598](https://orcid.org/0000-0002-1850-5598)
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param seq2pt OPTIONAL. A named vector mapping sequence IDs to patient IDs. Default is NULL – no patient labels will
#'               be added to the tree.
#' @param patient_label OPTIONAL. A logical value indicating whether to add patient labels to the tree. Default is FALSE
#'                      – no patient labels will be added to the tree.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint geom_tiplab
#' @importFrom hues iwanthue
#' @importFrom rlang .data
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes scale_color_manual theme element_text unit guides guide_legend
#' @importFrom stats setNames
#' @export
plot_clusters_phylo <- function(tree, clusters, seq2pt = NULL, patient_label = FALSE) {
    # convert phylo object to ggtree object
    tree <- ggtree(tree)

    # format the clusters into a dataframe
    cluster_df <- data.frame(isolate = names(clusters), clust_id = factor(clusters))

    # use iwanthue to generate colors for clusters
    cluster_colors <- setNames(
        iwanthue(length(unique(cluster_df$clust_id))),
        levels(cluster_df$clust_id)
    )

    # add cluster information to the tree
    tree$data <- tree$data |> left_join(cluster_df, by = c("label" = "isolate"))

    # build the plot
    p <- tree + geom_tippoint(aes(color = .data$clust_id), size = 2, alpha = 0.8)

    # add patient labels if requested
    if (!is.null(seq2pt) && patient_label) {
        p <- p +
            geom_tiplab(
                aes(label = paste0(label, " (", seq2pt[label], ")")),
                size = 2,
                offset = 0.5
            )
    }

    # return the final plot
    p +
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
    # unique cluster labels
    clusters1_unique_labels <- as.character(sort(unique(clusters1)))
    clusters2_unique_labels <- as.character(sort(unique(clusters2)))

    # create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        nrow = length(clusters1_unique_labels),
        ncol = length(clusters2_unique_labels),
        dimnames = list(clusters1_unique_labels, clusters2_unique_labels)
    )

    # compute cluster overlap
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
    # subset of isolates assigned to clusters
    clusters_subset <- clusters[clusters != 1]
    unique_clusters <- sort(unique(clusters_subset))

    # compute maximum genetic distance within cluster, minimum genetic distance to another cluster,
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

    # define function for generating plots
    generate_plot <- function(x, y, xlab, ylab, file_name) {
        # create directory if it doesn't exist
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

    # generate plots
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

    # create and return result matrix
    intra_inter_cluster_dist_mat <- t(cluster_distances)
    row.names(intra_inter_cluster_dist_mat) <- unique_clusters
    colnames(intra_inter_cluster_dist_mat) <- c(
        "Max intra-clust",
        "Min inter-clust",
        "Min inter-isolate"
    )

    intra_inter_cluster_dist_mat
}

#' Plots a trace heatmap for all patients for a sequence type.
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param dates A named vector of isolate dates by sequence ID.
#' @param ip_seqs A vector of sequence IDs which correspond to intake patient sequences presumed to be imported.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param trace_data A data frame or matrix representing daily patient trace (rows) by patient ID (columns).
#' @param trace_colors A vector of colors for the trace, one for each trace value.
#' @param surv_colors A vector of two colors for the surveillance, negative and positive.
#'
#' @return A `ggplot` object with the trace heatmap for all patients for a sequence type.
#'
#' @importFrom ape keep.tip
#' @importFrom dplyr group_by slice_min ungroup select mutate arrange pull
#' @importFrom paletteer paletteer_d
#' @importFrom ggplot2 geom_tile scale_fill_manual scale_fill_gradientn theme unit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggtree ggtree geom_tiplab gheatmap vexpand
#' @importFrom ggtreeExtra geom_fruit
#' @export
plot_st_patient_trace <- function(
    tree,
    clusters,
    snp_dist,
    dates,
    ip_seqs,
    seq2pt,
    trace_data,
    trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
    surv_colors = c("darkblue", "red")
) {
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
    isolate_order <- df |>
        group_by(patient) |>
        mutate(cluster_factor = factor(cluster, levels = cluster_order)) |>
        arrange(cluster_factor, date) |>
        slice_min(date, with_ties = FALSE) |>
        ungroup() |>
        pull(id)

    # get max SNP distance between all isolates for a single patient in that cluster
    get_max_snp <- function(iso_id) {
        clust <- df$cluster[df$id == iso_id]
        ptval <- df$patient[df$id == iso_id]
        sub_ids <- df$id[df$cluster == clust & df$patient == ptval]
        max(snp_dist[sub_ids, sub_ids])
    }
    max_dist <- vapply(isolate_order, get_max_snp, numeric(1))

    # subset the trace_data matrix
    pt_in_trace <- df$patient %in% row.names(trace_data)
    seq_ids_in_trace <- df$id[pt_in_trace]
    trace_sub <- trace_data[as.character(df$patient[pt_in_trace]), , drop = FALSE]
    row.names(trace_sub) <- seq_ids_in_trace

    # reorder rows by earliest isolates in cluster/date order
    keep_rows <- intersect(isolate_order, seq_ids_in_trace)
    trace_sub <- trace_sub[keep_rows, , drop = FALSE]

    # identify convert isolates
    is_convert_isolate <- sapply(keep_rows, function(id) {
        row_vals <- trace_sub[id, ]
        i_negative <- which(row_vals %in% c(1.25))
        i_positive <- which(row_vals %in% c(1.5))
        if (length(i_negative) == 0 || length(i_positive) == 0) {
            return(FALSE)
        }
        is_convert <- min(i_negative) < min(i_positive)
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
    cluster_levels <- sort(unique(df$cluster))
    annotation_row <- data.frame(
        id = keep_rows,
        Cluster = factor(df$cluster[match(keep_rows, df$id)], levels = cluster_levels),
        Convert = factor(convert_class, levels = patient_labels),
        Intra_pt_dist = max_dist[keep_rows] + 1
    )

    bluescale <- colorRampPalette(c("white", "blue"))
    ann_colors <- list(
        Cluster = setNames(iwanthue(length(cluster_levels)), cluster_levels),
        Convert = setNames(c("gray95", "red", "gray"), patient_labels),
        Intra_pt_dist = bluescale(max(max_dist, 1))
    )

    # prepare the column labels
    col_lab <- colnames(trace_sub)
    idx_keep_labels <- seq(1, length(col_lab), 14)
    col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

    # set up breaks
    custom_breaks <- sort(unique(as.numeric(as.matrix(trace_sub))))

    # calculate number of colors needed (excluding white for 0)
    n_colors <- length(custom_breaks) - 1

    # if n_colors > 13, recycle colors
    if (n_colors > length(trace_colors)) {
        trace_colors <- rep(trace_colors, length.out = n_colors)
        warning(
            paste0(
                "More than ",
                length(trace_colors),
                " colors needed for trace. Consider increasing the ",
                "number of colors in the palette."
            )
        )
    }

    # create final color vector with white for 0 value
    trace_custom_colors <- c(
        "white",
        trace_colors[1],
        surv_colors,
        trace_colors[2:length(trace_colors)]
    )

    tree_sub <- keep.tip(tree, keep_rows)
    trace_df <- as.data.frame(trace_sub)
    trace_df[] <- lapply(trace_df, function(x) {
        factor(as.numeric(as.character(x)), levels = custom_breaks)
    })

    # Create tip label mapping with patient IDs
    tip_label_map <- sapply(tree_sub$tip.label, function(tip) {
        # convert tip to character for consistent lookup
        tip_char <- as.character(tip)

        if (tip_char %in% names(seq2pt)) {
            paste0(seq2pt[tip_char], " (", tip, ")")
        } else {
            tip # use original tip label if no mapping exists
        }
    })

    # update all relevant objects with new tip labels
    tree_sub$tip.label <- tip_label_map
    annotation_row$id <- tip_label_map[annotation_row$id]
    row.names(trace_df) <- tip_label_map[row.names(trace_df)]

    # validate mappings and anchor y-order to tree tips
    if (any(is.na(annotation_row$id))) {
        stop(
            "Some annotation IDs could not be mapped to tree tip labels. Check seq2pt and input IDs."
        )
    }
    if (!all(row.names(trace_df) %in% tree_sub$tip.label)) {
        stop("Trace row names are not all present in tree tip labels after relabeling.")
    }
    annotation_row$id <- factor(annotation_row$id, levels = tree_sub$tip.label)
    # align heatmap rows to tree tip order to lock row positions
    trace_df <- trace_df[tree_sub$tip.label, , drop = FALSE]

    tree_plot <- ggtree(tree_sub) +
        geom_tiplab(size = 0.75, align = TRUE, offset = 1.5, linetype = NULL)

    gheatmap(
        tree_plot,
        trace_df,
        offset = 4,
        width = 4,
        font.size = 1,
        custom_column_labels = col_lab,
        colnames_angle = 90,
        colnames_offset_y = -2
    ) +
        vexpand(0.1, -1) +
        # trace annotation
        scale_fill_manual(
            name = "Trace",
            values = trace_custom_colors,
            breaks = custom_breaks,
            labels = custom_breaks,
            na.value = "white",
            drop = FALSE
        ) +
        new_scale_fill() +
        # cluster annotation
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
            limits = levels(annotation_row$Cluster),
            guide = "none",
            drop = FALSE
        ) +
        new_scale_fill() +
        # convert annotation
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Convert),
            offset = -0.19,
            width = 1
        ) +
        scale_fill_manual(
            name = "Convert",
            values = ann_colors$Convert,
            labels = names(ann_colors$Convert),
            limits = levels(annotation_row$Convert),
            drop = FALSE
        ) +
        new_scale_fill() +
        # intra_pt_dist annotation
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Intra_pt_dist),
            offset = -0.18,
            width = 1
        ) +
        scale_fill_gradientn(name = "Intra_pt_dist", colors = ann_colors$Intra_pt_dist) +
        theme(
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.position = "right",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm")
        )
}

#' Plots a trace heatmap for all patients for a cluster.
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param cluster_seqs A vector of sequence IDs which correspond to the cluster.
#' @param dates A named vector of isolate dates by sequence ID.
#' @param ip_seqs A vector of sequence IDs which correspond to intake patient sequences presumed to be imported.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param trace_data A data frame or matrix representing daily patient trace (rows) by patient ID (columns).
#' @param trace_colors A vector of colors for the trace, one for each trace value.
#' @param surv_colors A vector of two colors for the surveillance, negative and positive.
#'
#' @return A `ggplot` object with the trace heatmap for all patients for a cluster.
#'
#' @importFrom ape keep.tip
#' @importFrom dplyr group_by slice_min ungroup select mutate arrange pull
#' @importFrom paletteer paletteer_d
#' @importFrom ggplot2 geom_tile scale_fill_manual scale_fill_gradientn theme unit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggtree ggtree geom_tiplab gheatmap vexpand
#' @importFrom ggtreeExtra geom_fruit
#' @export
plot_cluster_patient_trace <- function(
    tree,
    cluster_seqs,
    dates,
    ip_seqs,
    seq2pt,
    trace_data,
    trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
    surv_colors = c("darkblue", "red", "yellow")
) {
    # prepare a data frame for all isolates in the cluster
    df <- data.frame(
        id = names(cluster_seqs),
        patient = seq2pt[names(cluster_seqs)],
        date = dates[names(cluster_seqs)]
    )

    # subset the trace_data matrix
    pt_in_trace <- df$patient %in% row.names(trace_data)
    seq_ids_in_trace <- df$id[pt_in_trace]
    trace_sub <- trace_data[as.character(df$patient[pt_in_trace]), , drop = FALSE]
    row.names(trace_sub) <- seq_ids_in_trace

    # for each patient, find earliest isolate in time
    # then map each isolate to that earliest isolate
    isolate_order <- df |>
        group_by(patient) |>
        arrange(date) |>
        slice_min(date, with_ties = FALSE) |>
        ungroup() |>
        pull(id)

    # reorder rows by earliest isolates in cluster/date order
    keep_rows <- intersect(isolate_order, seq_ids_in_trace)
    trace_sub <- trace_sub[keep_rows, , drop = FALSE]

    # for each isolate for each patient in this cluster, set the trace value to 1.75
    for (id in keep_rows) {
        # find patient for this id
        pt_id <- df$patient[df$id == id]
        for (isolate in df$id[df$patient == pt_id]) {
            trace_sub[id, as.character(df$date[df$id == isolate])] <- 1.75
        }
    }

    # identify convert isolates
    is_convert_isolate <- sapply(keep_rows, function(id) {
        row_vals <- trace_sub[id, ]
        i_negative <- which(row_vals %in% c(1.25))
        i_positive <- which(row_vals %in% c(1.5, 1.75))
        if (length(i_negative) == 0 || length(i_positive) == 0) {
            return(FALSE)
        }
        is_convert <- min(i_negative) < min(i_positive)
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
        Convert = factor(convert_class, levels = patient_labels)
    )

    # annotation colors
    ann_colors <- list(Convert = setNames(c("gray95", "red", "gray"), patient_labels))

    # prepare the column labels
    col_lab <- colnames(trace_sub)
    idx_keep_labels <- seq(1, length(col_lab), 14)
    col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

    # set up breaks
    custom_breaks <- sort(unique(as.numeric(as.matrix(trace_sub))))

    # calculate number of colors needed (excluding white for 0)
    n_colors <- length(custom_breaks) - 1

    # if n_colors > length(trace_colors), error out
    if (n_colors > length(trace_colors)) {
        warning(
            paste0(
                "More than ",
                length(trace_colors),
                " colors needed for trace. Consider increasing the ",
                "number of colors in the palette."
            )
        )
    }

    # create final color vector with white for 0 value
    trace_custom_colors <- c(
        "white",
        trace_colors[1],
        surv_colors,
        trace_colors[2:length(trace_colors)]
    )

    tree_sub <- keep.tip(tree, keep_rows)
    trace_df <- as.data.frame(trace_sub)
    trace_df[] <- lapply(trace_df, function(x) {
        factor(as.numeric(as.character(x)), levels = custom_breaks)
    })

    # create tip label mapping with patient IDs
    tip_label_map <- sapply(tree_sub$tip.label, function(tip) {
        # convert tip to character for consistent lookup
        tip_char <- as.character(tip)

        if (tip_char %in% names(seq2pt)) {
            paste0(seq2pt[tip_char], " (", tip, ")")
        } else {
            tip # use original tip label if no mapping exists
        }
    })

    # Update all relevant objects with new tip labels
    tree_sub$tip.label <- tip_label_map
    annotation_row$id <- tip_label_map[annotation_row$id]
    row.names(trace_df) <- tip_label_map[row.names(trace_df)]

    # validate mappings and anchor y-order to tree tips
    if (any(is.na(annotation_row$id))) {
        stop(
            "Some annotation IDs could not be mapped to tree tip labels. Check seq2pt and input IDs."
        )
    }
    if (!all(row.names(trace_df) %in% tree_sub$tip.label)) {
        stop("Trace row names are not all present in tree tip labels after relabeling.")
    }
    annotation_row$id <- factor(annotation_row$id, levels = tree_sub$tip.label)

    tree_plot <- ggtree(tree_sub) +
        geom_tiplab(size = 0.75, align = TRUE, offset = 1.5, linetype = NULL)

    gheatmap(
        tree_plot,
        trace_df,
        offset = 2,
        width = 4,
        font.size = 1,
        custom_column_labels = col_lab,
        colnames_angle = 90
    ) +
        vexpand(0.1, -1) +
        # trace annotation
        scale_fill_manual(
            name = "Trace",
            values = trace_custom_colors,
            breaks = custom_breaks,
            labels = custom_breaks,
            na.value = "white",
            drop = FALSE
        ) +
        new_scale_fill() +
        # convert annotation
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Convert),
            offset = -0.19,
            width = 1
        ) +
        scale_fill_manual(
            name = "Convert",
            values = ann_colors$Convert,
            labels = names(ann_colors$Convert),
            limits = levels(annotation_row$Convert),
            drop = FALSE
        ) +
        theme(
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.position = "right",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm")
        )
}
