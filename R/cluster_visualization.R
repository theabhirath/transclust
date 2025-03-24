#' Maps clusters onto a phylogenetic tree and visualizes them
#'
#' @author [Aryan Singh](mailto:aryansin@umich.edu)
#'
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param tree A phylogenetic tree object of class `phylo`.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint
#' @importFrom hues iwanthue
#' @importFrom dplyr .data left_join
#' @importFrom ggplot2 aes scale_color_manual ggtitle theme element_text unit guides guide_legend
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

#' Produces summary figures and information regarding the clusters.
#'
#' @param clusters A named vector indicating the cluster (subtree) for each sequence.
#' @param sub_trees A list of phylogenetic subtrees.
#' @param dna_aln A DNA object corresponding to sequence IDs.
#' @param pt_trace A data frame or matrix representing daily patient trace (rows) by patient ID (columns).
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_pt_seqs A vector of sequence IDs which correspond to intake positive patients.
#' @param ip_seqs A vector of sequence IDs which correspond to intake positive patients presumed to be imported.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param dates A named vector of isolate dates by sequence ID.
#' @param prefix A descriptor for naming figure outputs.
#' @param floor_trace  An optional floor trace where the rows are days and the columns are patients.
#'
#' @return No return value, invoked for side effects.
#'
#' @importFrom stats heatmap
#' @export
cluster_summary <- function(clusters, sub_trees, dna_aln, pt_trace, seq2pt, ip_pt_seqs,
                            ip_seqs, snp_dist, dates, prefix, floor_trace = NULL) {
    #####################################################################################
    # Plot the size distribution of clusters #####
    #####################################################################################
    cluster_sizes <- table(clusters)
    # Remove singleton clusters
    cluster_sizes <- cluster_sizes[names(cluster_sizes) != 1]
    # Early return if no clusters with more than one member
    if (length(cluster_sizes) == 0) {
        print("No clusters with more than one member!")
        return()
    }
    file <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_size_dist.pdf")
    pdf(file)
    hist(cluster_sizes, 20, main = "Cluster Size Distribution")
    dev.off()

    #####################################################################################
    # Bar plot of cluster sizes showing separate counts of: #####
    #   - Converts (non-IP)
    #   - Intake-positive isolate (ip_seq_count)
    #   - Intake-positive patient (ip_pt_seq_count)
    #####################################################################################
    # Helper function to count how many sequences from ip_vector fall into each cluster
    # By default this counts the number of isolates, but can be set to count the number of patients
    count_ip_in_cluster <- function(cluster_id, ip_vector, patient = FALSE) {
        seqs_in_cluster <- names(clusters)[clusters == cluster_id]
        found_seqs <- intersect(ip_vector, seqs_in_cluster)
        length(if (patient) unique(seq2pt[found_seqs]) else found_seqs)
    }
    ip_pt_seq_count <- sapply(names(cluster_sizes), count_ip_in_cluster, ip_vector = ip_pt_seqs)
    ip_seq_count <- sapply(names(cluster_sizes), count_ip_in_cluster, ip_vector = ip_seqs)

    # Plot bar plot of clusters showing number of intake positive/convert isolates
    cluster_sort <- sort.int(cluster_sizes, index.return = TRUE)
    sorted_ix <- cluster_sort$ix
    # cluster_w_ip_count organizes row data:
    #  1st row: total cluster size - # of ip_pt_seq_count = converts?
    #  2nd row: intake-positive isolate count
    #  3rd row: difference between ip_pt_seq_count and ip_seq_count (i.e. patient-level vs isolate-level)
    cluster_w_ip_count <- rbind(
        cluster_sort$x - ip_pt_seq_count[sorted_ix],
        ip_seq_count[sorted_ix],
        ip_pt_seq_count[sorted_ix] - ip_seq_count[sorted_ix]
    )
    file <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_size_vs_ip_count_barplot.pdf")
    pdf(file)
    barplot(cluster_w_ip_count,
        xlab = "Clusters", ylab = "Number of isolates in cluster",
        legend.text = c("Converts", "Admit positive isolate", "Admit positive patient"),
        col = c("lightblue", "white", "lightgray"),
        names.arg = rep("", ncol(cluster_w_ip_count)),
        args.legend = list(x = "topleft")
    )
    dev.off()

    #####################################################################################
    # Bar plot of clusters showing number of intake-positive/convert patients. #####
    # This considers the unique patients rather than the number of isolates.
    #####################################################################################
    cluster_sizes_pts <- vapply(names(cluster_sizes), function(c) {
        length(unique(sort(seq2pt[names(clusters[clusters == c])])))
    }, integer(1))
    ip_seq_count_pts <- sapply(names(cluster_sizes), count_ip_in_cluster,
        ip_vector = ip_seqs, patient = TRUE
    )
    ip_pt_seq_count_pts <- sapply(names(cluster_sizes), count_ip_in_cluster,
        ip_vector = ip_pt_seqs, patient = TRUE
    )
    cluster_sort_pts <- sort.int(cluster_sizes_pts, index.return = TRUE)
    sorted_ix_pts <- cluster_sort_pts$ix
    cluster_w_ip_count_pts <- rbind(
        cluster_sort_pts$x - ip_pt_seq_count_pts[sorted_ix_pts],
        ip_seq_count_pts[sorted_ix_pts],
        ip_pt_seq_count_pts[sorted_ix_pts] - ip_seq_count_pts[sorted_ix_pts]
    )
    file <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_size_vs_ip_count_pts_barplot.pdf")
    pdf(file)
    barplot(cluster_w_ip_count_pts,
        xlab = "Clusters", ylab = "Number of patients in cluster",
        legend.text = c("Converts", "Admit positive isolate", "Admit positive patient"),
        col = c("lightblue", "white", "lightgray"),
        names.arg = rep("", ncol(cluster_w_ip_count_pts)),
        args.legend = list(x = "topleft")
    )
    dev.off()

    #####################################################################################
    # For each cluster, plot: #####
    #   - Pairwise genetic distance distribution
    #   - Patient trace ordered by phylogeny
    #   - Heatmap of variants
    #####################################################################################
    out_group <- which.max(rowMeans(snp_dist))
    dna_var_reref <- t(apply(dna_aln, 1, function(x) x != dna_aln[out_group, ]))
    row.names(dna_var_reref) <- row.names(dna_aln)

    # Iterate over the clusters omitting the default cluster
    for (c in setdiff(unique(clusters), 1)) {
        # Omit clusters with only one member
        if (length(clusters[clusters == c]) == 1) {
            print(paste0("Cluster ", c, " has only one member!"))
            next()
        }

        # # Check that cluster contains all members of sub-clade
        # if (length(sub_trees[[c]]$tip.label) != sum(clusters == c)) {
        #     print(paste0("Cluster ", c, " does not contain all members of sub-clade!!"))
        #     # next; # CHECK: why is this commented out?
        # }

        # Get patients corresponding to cluster members
        cluster_members <- names(clusters)[clusters == c]
        cluster_pts <- seq2pt[as.character(cluster_members)]
        labels <- paste0(cluster_members, "(", cluster_pts, ")")

        # Get variants to plot
        var_pos <- which(colSums(dna_var_reref[cluster_members, ]) > 0 &
                             colSums(dna_var_reref[cluster_members, ]) < length(cluster_members))
        dna_cluster_bin <- t(dna_var_reref[cluster_members, var_pos])

        row.names(dna_cluster_bin) <- var_pos
        colnames(dna_cluster_bin) <- labels

        # # Determine which variants agree with subtree
        # is_mp <- apply(rbind(dna_aln[out_group, var_pos], dna_cluster), 2,
        #     FUN = function(x) {
        #         is.monophyletic(sub_trees[[c]], labels(dna_cluster)[x[2:length(x)] != x[1]])
        #     }
        # )
        colSide <- rep("white", nrow(dna_cluster_bin))
        # colSide[!is_mp] <- 'red'
        rowSide <- rep("white", ncol(dna_cluster_bin))
        rowSide[!(cluster_members %in% ip_seqs)] <- "#96a1a4"

        # Plot variant heatmap
        if (sum(rowSums(dna_cluster_bin) < (ncol(dna_cluster_bin) - 1)) > 1) {
            file <- paste0(
                "figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix,
                "_cluster", c, "_variant_map_nonMPred.pdf"
            )
            pdf(file)
            heatmap(
                apply(dna_cluster_bin[rowSums(dna_cluster_bin) < (ncol(dna_cluster_bin) - 1), ], 1, as.numeric),
                scale = "none", col = c("lightgray", "black"),
                distfun = function(x) dist(x, method = "manhattan"),
                RowSideColors = rowSide,
                ColSideColors = colSide[rowSums(dna_cluster_bin) < (ncol(dna_cluster_bin) - 1)],
                Colv = NA, cexCol = .3, labRow = colnames(dna_cluster_bin)
            )
            dev.off()
        }

        # Plot patient trace for the cluster
        plot_pt_trace(
            cluster_members, dates, ip_seqs, pt_trace, seq2pt, snp_dist, paste0(prefix, "_cluster", c)
        )

        # Plot floor trace if given
        if (!is.null(floor_trace) && length(unique(cluster_pts)) > 1) {
            plot_floor_trace(
                clusters[clusters == c], dates, ip_seqs, floor_trace,
                seq2pt, snp_dist, paste0(prefix, "_cluster", c)
            )
        }
    }
}

#' Plots the patient trace broken up by SNP threshold and labeled using designated clusters.
#'
#' @param dna_aln A DNA object with names corresponding to sequence IDs
#' @param pt_trace A data frame or matrix representing daily patient trace (rows) by patient ID (columns)
#' @param seq2pt A named vector mapping sequence IDs to patient IDs
#' @param snp_thresh Numeric threshold for breaking up trace heatmaps
#' @param ip_pt_seqs Vector of sequence IDs for intake positive patients
#' @param ip_seqs Vector of sequence IDs presumed to be imported
#' @param snp_dist A matrix of SNP distances between isolates
#' @param dates A named vector of isolate dates by sequence ID
#' @param prefix A string descriptor used for naming figure outputs
#' @param floor_trace Optional data for floor-level tracing
#'
#' @return A numeric vector indicating the cluster that each isolate belongs to.
#'
#' @details First, we cluster the isolates using a hard SNP distance threshold. Then, we plot
#' the patient trace heatmaps for each cluster. If a floor trace is provided, we also plot the
#' floor trace for each cluster.
#'
#' @importFrom stats hclust as.dist
#' @export
cluster_context_2 <- function(dna_aln, pt_trace, seq2pt, snp_thresh, ip_pt_seqs, ip_seqs,
                              snp_dist, dates, prefix, floor_trace = NULL) {
    # Get clusters at desired SNP threshold
    clusters <- get_tn_clusters_snp_thresh(dna_aln, snp_dist, snp_thresh)
    # Plot trace heatmaps
    for (c in setdiff(unique(clusters), 1)) {
        # Get patients corresponding to cluster members
        cluster_members <- names(clusters)[clusters == c]
        cluster_pts <- seq2pt[as.character(cluster_members)]
        # Plot trace
        plot_pt_trace(
            cluster_members, dates, ip_seqs, pt_trace, seq2pt, snp_dist,
            paste0(prefix, "_", snp_thresh, "_thresh_context_cluster_", c)
        )
        # Plot floor trace if given
        if (!is.null(floor_trace) && length(unique(cluster_pts)) > 1) {
            plot_floor_trace(
                clusters[clusters == c], dates, ip_seqs, floor_trace, seq2pt,
                snp_dist, paste0(prefix, "_", snp_thresh, "_thresh_context_cluster_", c)
            )
        }
    }

    # Return the clusters
    clusters
}

#' @importFrom stats dist hclust as.dendrogram
#' @importFrom ape keep.tip
#' @importFrom phytools nodeHeights
#' @importFrom dplyr mutate mutate_all
#' @importFrom ggtree ggtree geom_tiplab gheatmap geom_facet vexpand
#' @importFrom ggplot2 scale_fill_manual aes theme unit ggsave
plot_pt_trace <- function(sequence_ids, tree, dates, ip_seqs, pt_trace, seq2pt, snp_dist, prefix) {
    # Convert Excel-style numeric dates to R Date objects just once
    date_vec <- as.Date(dates, origin = "1899-12-30")

    # Build a data frame for all sequences
    # Each row has: sequence ID, patient ID, and the date
    df <- data.frame(
        id = sequence_ids,
        pt = seq2pt[sequence_ids],
        date = date_vec[sequence_ids]
    )

    # Identify which sequences/patients are present in pt_trace's columns
    df$in_trace <- df$pt %in% colnames(pt_trace)

    # Subset df to keep only those in pt_trace
    df_sub <- df[df$in_trace, ]

    # Build the trace_sub matrix by transposing only the necessary columns
    trace_sub <- t(pt_trace[, as.character(df_sub$pt), drop = FALSE])
    rownames(trace_sub) <- df_sub$id

    # Mark each isolate's date with a distinct value (2.5) in trace_sub
    # cbind() pairs the row name (isolate ID) and the date as a character
    trace_sub[cbind(df_sub$id, as.character(df_sub$date))] <- 2.5

    # Subset the tree provided to keep only the relevant sequences
    tree_sub <- keep.tip(tree, df_sub$id)

    # Row coloring: white for index patient isolates, lightblue otherwise
    row_colors <- ifelse(df_sub$id %in% ip_seqs, "white", "lightblue")

    # Column labeling: only label every 14th column i.e. every two weeks
    col_lab <- rownames(pt_trace)
    idx_labeled <- seq(1, length(col_lab), 14)
    col_lab[setdiff(seq_along(col_lab), idx_labeled)] <- ""

    # Row labels: ID and patient ID in parentheses
    labels <- paste0(df_sub$id, "(", df_sub$pt, ")")

    # Create file name, including approximate dendrogram height
    dend_height <- max(nodeHeights(tree_sub))
    file <- paste0(
        "figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix,
        "_trace_map_by_sequence_scale_", dend_height, "snps.pdf"
    )
    p <- gheatmap(
        ggtree(tree_sub),
        trace_sub |> as.data.frame() |> mutate_all(as.factor),
        color = NULL,
        width = 4,
        offset = 1,
        font.size = 2,
        custom_column_labels = col_lab,
        colnames_angle = 90
    ) +
        # adjust for long vertical column labels
        vexpand(.1, -1) +
        # color the tiles
        scale_fill_manual(
            values = c("white", "lightgray", "blue", "red", "yellow"),
            name = "Trace"
        ) +
        # add colors to the rows
        geom_facet(
            data = data.frame(id = df_sub$id, color = row_colors),
            aes(x = id, y = 0, fill = color),
            geom = "tile", inherit.aes = FALSE,
            width = 1, height = 0.5
        )
    ggsave(
        filename = file,
        plot = p,
        width = 50,
        height = 50,
        dpi = "retina",
        units = "cm"
    )
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr group_by slice_min ungroup rename left_join mutate arrange select
#' @keywords internal
plot_floor_trace <- function(clusters, dates, ip_seqs, floor_trace, seq2pt, snp_dist, prefix, cluster_class = NULL) {
    # Convert numeric Excel-style dates to R Date objects only once
    date_vec <- as.Date(dates, origin = "1899-12-30")

    # Prepare a data frame for all isolates
    df <- data.frame(
        id = names(clusters),
        cluster = as.character(clusters),
        pt = seq2pt[names(clusters)],
        date = date_vec[names(clusters)]
    )

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
        group_by(.data$cluster, .data$pt) |>
        slice_min(.data$date, with_ties = FALSE) |>
        ungroup() |>
        rename(earliest_id = .data$id)

    # Add earliest isolate info back to df
    df <- df |>
        left_join(df_earliest |> select(.data$cluster, .data$pt, .data$earliest_id), by = c("cluster", "pt"))

    # Identify isolate order: cluster_order, then by earliest isolate's date
    df_earliest <- df_earliest |>
        mutate(cluster_factor = factor(.data$cluster, levels = cluster_order)) |>
        arrange(.data$cluster_factor, .data$date)

    isolate_order <- df_earliest$earliest_id

    # Identify if each earliest isolate is in the user-supplied index set
    df_earliest <- df_earliest |> mutate(is_index = .data$earliest_id %in% ip_seqs)

    # Create named vector for rep_isolate: for each isolate, the earliest isolate in the same cluster & patient
    rep_isolate <- setNames(df$earliest_id, df$id)

    # Helper function to get max SNP distance within cluster or within cluster–patient
    get_max_snp <- function(iso_id, use_entire_cluster = FALSE) {
        clust <- df$cluster[df$id == iso_id]
        ptval <- df$pt[df$id == iso_id]
        sub_ids <- if (!use_entire_cluster) {
            # same cluster, same patient
            df$id[df$cluster == clust & df$pt == ptval]
        } else {
            # entire cluster
            df$id[df$cluster == clust]
        }
        max(snp_dist[sub_ids, sub_ids])
    }

    # For each earliest isolate, compute max distances
    max_dist <- sapply(isolate_order, get_max_snp, use_entire_cluster = FALSE)
    max_clust_dist <- sapply(isolate_order, get_max_snp, use_entire_cluster = TRUE)

    # Subset the floor_trace matrix
    sequence_ids <- df$id
    pt_in_trace <- seq2pt[sequence_ids] %in% colnames(floor_trace)
    trace_sub <- t(floor_trace[, as.character(seq2pt[sequence_ids[pt_in_trace]])])
    rownames(trace_sub) <- sequence_ids[pt_in_trace]

    # Mark the actual isolate day with 1.75
    these_ids <- sequence_ids[pt_in_trace]
    these_dates <- date_vec[these_ids]
    trace_sub[cbind(rep_isolate[these_ids], as.character(these_dates))] <- 1.75

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
        trace_date <- as.Date(colnames(trace_sub)[min(i_positive)])
        iso_date <- date_vec[id]
        (iso_date - trace_date < 7) && is_convert && !(id %in% ip_seqs)
    })

    # Index isolates
    is_index_isolate <- keep_rows %in% ip_seqs

    # Build annotation data
    patient_labels <- c("Secondary convert", "Index patient", "Convert patient")
    convert_class <- rep("Secondary convert", length(keep_rows))
    convert_class[is_index_isolate] <- "Index patient"
    convert_class[!is_index_isolate & is_convert_isolate] <- "Convert patient"

    annotation_row <- data.frame(
        Cluster          = df$cluster[match(keep_rows, df$id)],
        Cluster_class    = cluster_class[df$cluster[match(keep_rows, df$id)]],
        Convert          = convert_class,
        Intra_pt_dist    = max_dist[keep_rows] + 1,
        Intra_clust_dist = max_clust_dist[keep_rows] + 1
    )
    row.names(annotation_row) <- keep_rows

    # Annotation colors
    cluster_colors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
    bluescale <- colorRampPalette(c("white", "blue"))

    n_class <- length(unique(cluster_class))
    cluster_class_colors <- if (n_class > 0) {
        if (n_class > 3) {
            brewer.pal(min(n_class, 8), "Dark2")
        } else {
            colorRampPalette(c("lightgray", "darkgray"))(n_class)
        }
    } else {
        character(0)
    }

    ann_colors <- list(
        Cluster          = structure(cluster_colors[seq_along(unique(df$cluster))], names = sort(unique(df$cluster))),
        Cluster_class    = structure(cluster_class_colors, names = unique(cluster_class)),
        Convert          = structure(c("gray95", "gray", "red"), names = patient_labels),
        Intra_pt_dist    = bluescale(max(max_dist, 1)),
        Intra_clust_dist = bluescale(max(max_clust_dist, 1))
    )

    # Column labeling
    col_lab <- rownames(floor_trace)
    idx_keep_labels <- seq(1, length(col_lab), 14)
    col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

    # Row labels
    labels <- paste(keep_rows, "(", seq2pt[keep_rows], ")", sep = "")

    # Breaks and colors for pheatmap
    custom_breaks <- c(-1, 0, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6)
    custom_colors <- brewer.pal(7, "Pastel2")
    custom_colors <- c("white", custom_colors[1], "blue", "red", "yellow", custom_colors[c(2:5, 7)])

    # Gaps between clusters
    ctab <- table(annotation_row$Cluster)
    ctab <- ctab[as.character(cluster_order)]
    row_gaps <- cumsum(ctab)
    row_gaps <- row_gaps[row_gaps < nrow(annotation_row)]

    # Font size scaling
    fsize_row <- 2 * (1 - (nrow(trace_sub) / ncol(floor_trace)))

    # Output file
    file <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_floor_trace_map_w_clusters.pdf")
    pheatmap(
        trace_sub,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        breaks = custom_breaks,
        color = custom_colors,
        scale = "none",
        labels_col = col_lab,
        labels_row = labels,
        legend = TRUE,
        gaps_row = row_gaps,
        fontsize_row = fsize_row,
        fontsize_col = 2,
        fontsize = 2,
        annotation_row = annotation_row,
        annotation_colors = ann_colors,
        filename = file
    )
}

#' Produces summary statistics and diagnostic plots for transmission clusters.
#'
#' @param clusters A named vector of cluster assignments, where the names are
#' sequence IDs and the values are cluster IDs.
#' @param pt_trace A data frame or matrix of patient trace data, with rows
#' representing days and columns representing patient IDs.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_pt_seqs A vector of sequence IDs corresponding to intake positive
#' patients.
#' @param ip_seqs A vector of sequence IDs presumed to be imported.
#' @param snp_dist A matrix of SNP distances between isolates.
#' @param dates A named vector of isolate dates by sequence ID.
#' @param prefix A descriptor for naming figure outputs.
#' @param floor_trace (Optional) Data for floor-level tracing (rows: days,
#' columns: patients).
#' @param room_trace (Optional) Data for room-level tracing (rows: days,
#' columns: patients).
#'
#' @return A matrix whose rows correspond to cluster IDs (excluding singletons)
#' and whose columns contain various properties computed for each cluster,
#' such as the number of patients in the cluster, the total number of converts,
#' the cluster duration in days, and genetic distance statistics.
#'
#' @details Produce summary statistics and diagnostic plots for transmission clusters.
#' This function calculates summary properties (e.g., cluster size, number of
#' converts, genetic distance statistics) for each transmission cluster and
#' generates multiple diagnostic figures describing the distribution of these
#' properties.

#' @export
cluster_property_report <- function(clusters, pt_trace, seq2pt, ip_pt_seqs, ip_seqs, snp_dist,
                                    dates, prefix, floor_trace = NULL, room_trace = NULL) {
    # First, merge all singleton clusters (size == 1) into cluster ID 1 i.e. default
    # This step ensures we only consider meaningful clusters with more than 1
    # member as we generate the summary metrics.
    cluster_size <- table(clusters)
    singletons <- names(cluster_size)[cluster_size == 1]
    clusters[clusters %in% singletons] <- 1

    # Identify cluster IDs that remain after excluding the default singleton clusters
    cluster_ids <- sort(unique(clusters[clusters != 1]))

    # Calculate summary properties for each cluster using an auxiliary helper:
    # `cluster_properties`
    cluster_props <- t(
        sapply(cluster_ids, function(clust_id) {
            message("Processing cluster", clust_id)
            cluster_properties(
                names(clusters)[clusters == clust_id], pt_trace, seq2pt, ip_pt_seqs,
                ip_seqs, snp_dist, dates, floor_trace, room_trace
            )
        })
    )

    rownames(cluster_props) <- cluster_ids

    # Only keep clusters that have at least one convert.
    has_converts <- cluster_props[, "Number_of_converts"] > 0
    cluster_props <- cluster_props[has_converts, ]

    #---------------------------------------------------------------------------
    # Figure 1: Distribution of cluster size, number of converts, and duration.
    #---------------------------------------------------------------------------
    figure_file_1 <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_size_and_duration.pdf")
    pdf(figure_file_1)
    on.exit(dev.off(), add = TRUE) # Safeguard to close the PDF device on exit.
    par(mfrow = c(2, 2))

    # Count singletons for the cluster size histogram.
    singleton_count <- sum(clusters[!(names(clusters) %in% ip_pt_seqs)] == 1)

    # Histogram of cluster sizes, including singletons as size=1.
    hist(
        c(rep(1, singleton_count), cluster_props[, "Number_of_patients"]),
        seq(0, max(cluster_props[, "Number_of_patients"]), 1),
        main = "Cluster size",
        xlab = "Cluster size",
        ylab = "Number of clusters"
    )

    # Histogram of cluster durations.
    hist(
        cluster_props[, "Cluster_duration"],
        breaks = 20,
        main = "Duration of clusters",
        xlab = "Cluster duration (days)",
        ylab = "Number of clusters"
    )

    # Scatterplot of cluster size vs. number of converts.
    plot(
        jitter(cluster_props[, "Number_of_patients"], 0.5),
        jitter(cluster_props[, "Number_of_converts"], 0.5),
        main = "Number of patients vs. number of converts",
        xlab = "Cluster size",
        ylab = "Number of converts"
    )

    # Scatterplot of cluster size vs. cluster duration.
    plot(
        jitter(cluster_props[, "Number_of_patients"], 0.5),
        jitter(cluster_props[, "Cluster_duration"], 0.5),
        main = "Number of patients vs. duration of cluster",
        xlab = "Cluster size",
        ylab = "Cluster duration (days)"
    )

    #---------------------------------------------------------------------------
    # Figure 2: Genetic variability within clusters.
    #---------------------------------------------------------------------------
    figure_file_2 <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_cluster_genetic_var.pdf")
    pdf(figure_file_2)
    on.exit(dev.off(), add = TRUE) # Ensure all devices close properly.
    par(mfrow = c(2, 2))

    # Histogram of median pairwise genetic distance.
    hist(
        cluster_props[, "Median_genetic_distance"],
        breaks = 20,
        main = "Median pairwise genetic distance",
        xlab = "Median within-cluster distance",
        ylab = "Number of clusters"
    )

    # Histogram of maximum pairwise genetic distance.
    hist(
        cluster_props[, "Max_genetic_distance"],
        breaks = 20,
        main = "Maximum pairwise genetic distance",
        xlab = "Max within-cluster distance",
        ylab = "Number of clusters"
    )

    # Scatterplot of cluster size vs. median within-cluster genetic distance.
    plot(
        jitter(cluster_props[, "Number_of_patients"], 0.5),
        jitter(cluster_props[, "Median_genetic_distance"], 0.5),
        main = "Cluster size vs. median genetic distance",
        xlab = "Cluster size",
        ylab = "Median within-cluster distance"
    )

    # Scatterplot of cluster duration vs. median genetic distance.
    plot(
        jitter(cluster_props[, "Cluster_duration"], 0.5),
        jitter(cluster_props[, "Median_genetic_distance"], 0.5),
        main = "Duration of cluster vs. median genetic distance",
        xlab = "Cluster duration (days)",
        ylab = "Median within-cluster distance"
    )

    #---------------------------------------------------------------------------
    # Figure 3: Epidemiological summary stats.
    #---------------------------------------------------------------------------
    figure_file_3 <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_epi_summary.pdf")
    pdf(figure_file_3)
    on.exit(dev.off(), add = TRUE)

    n_converts_total <- sum(cluster_props[, "Number_of_converts"])

    # Fraction of clusters that contain at least one start index + at least one convert
    perc_index_clusters <- sum(
        cluster_props[, "Number_of_start_indexes"] > 0 & cluster_props[, "Number_of_converts"] > 0
    ) / sum(cluster_props[, "Number_of_converts"] > 0)

    # Fraction of clusters where all converts occur after the index
    perc_index_first <- sum(
        cluster_props[, "Number_of_converts"] == cluster_props[, "Number_of_converts_after_index"]
    ) / sum(cluster_props[, "Number_of_converts"] > 0)

    # Fraction of all converts that occur after the index
    perc_convert_after_index <- sum(cluster_props[, "Number_of_converts_after_index"]) / n_converts_total

    # Fraction of converts who have a known putative source
    perc_overlap <- sum(cluster_props[, "Number_of_converts_with_source"]) / n_converts_total

    barplot(
        c(
            perc_index_clusters,
            perc_index_first,
            perc_convert_after_index,
            perc_overlap
        ),
        names.arg = c(
            "%clusters with index",
            "% clusters with index first",
            "% conversions after index",
            "% conversion with overlap"
        ),
        ylim = c(0, 1)
    )

    #---------------------------------------------------------------------------
    # Figure 4: Explore cluster properties vs. cluster size or duration.
    #---------------------------------------------------------------------------
    figure_file_4 <- paste0("figures/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, "_epi_explanations.pdf")
    pdf(figure_file_4)
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 1))

    # Fraction of converts after index vs. cluster size
    plot(
        jitter(cluster_props[, "Number_of_patients"], 0.5),
        jitter(cluster_props[, "Number_of_converts_after_index"] / cluster_props[, "Number_of_converts"], 0.5),
        main = "Number of patients vs. known index",
        xlab = "Cluster size",
        ylab = "% converts after putative index"
    )
    # Fraction of converts with putative source vs. cluster size
    plot(
        jitter(cluster_props[, "Number_of_patients"], 0.5),
        jitter(cluster_props[, "Number_of_converts_with_source"] / cluster_props[, "Number_of_converts"], 0.5),
        main = "Number of patients vs. epi explanations",
        xlab = "Cluster size",
        ylab = "% converts with putative source"
    )

    # Return the final cluster property matrix to the user.
    cluster_props
}

#' Intra-Patient Genetic Variation Analysis
#'
#' @param clusters A named vector (by sequence IDs) with cluster IDs.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param ip_seqs A vector of sequence IDs corresponding to intake patient sequences (presumed imported).
#' @param ip_pt_seqs A vector of sequence IDs corresponding to intake positive patients.
#' @param snp_dist A matrix of SNP distances between isolates. Row and column names must be sequence IDs.
#' @param dates A named numeric (or Date) vector of isolate dates (named by sequence IDs).
#' @param prefix A descriptor string used to name the output figures.
#'
#' @return A two-column matrix with intra-cluster SNP distances and the corresponding patient class.
#'
#' @details This function computes intra-patient genetic diversity metrics within and between clusters.
#' For each patient with multiple sequences, it calculates the maximum pairwise SNP distance and the
#' corresponding time difference for both intra- and inter-cluster comparisons. Four plots are generated:
#' a scatter plot and three boxplots.
#'
#' @export
intra_pt_genetic_var_analysis <- function(clusters, seq2pt, ip_seqs, ip_pt_seqs, snp_dist, dates, prefix) {
    ## Remove sequences assigned to cluster "1"
    seq2pt <- seq2pt[clusters[names(seq2pt)] != 1]
    valid_seq_ids <- names(seq2pt)

    ## Identify patients with multiple sequences
    patients_multi_seq <- sort(names(which(table(seq2pt) > 1)))

    ## Helper function to compute maximum SNP distance and corresponding time difference (intra-cluster)
    compute_max_pair <- function(seqs) {
        combs <- t(combn(seqs, 2))
        dists <- snp_dist[combs[, 1], combs[, 2]]
        max_idx <- which.max(dists)
        c(
            distance = dists[max_idx],
            time_diff = abs(dates[combs[max_idx, 1]] - dates[combs[max_idx, 2]])
        )
    }

    ## Helper function for inter-cluster comparisons using two sets of sequences
    compute_max_pair_inter <- function(seqs1, seqs2) {
        grid_pairs <- expand.grid(seqs1, seqs2, stringsAsFactors = FALSE)
        dists <- snp_dist[grid_pairs[, 1], grid_pairs[, 2]]
        max_idx <- which.max(dists)
        c(
            distance = dists[max_idx],
            time_diff = abs(dates[grid_pairs[max_idx, 1]] - dates[grid_pairs[max_idx, 2]])
        )
    }

    ## Helper function to determine patient class for a set of sequences:
    ## returns TRUE (index) if any sequence is in ip_seqs, FALSE (convert) otherwise;
    ## returns NA if all sequences are convert.
    compute_patient_class <- function(seqs) {
        if (sum(seqs %in% ip_pt_seqs & !(seqs %in% ip_seqs)) < length(seqs)) {
            any(seqs %in% ip_seqs)
        } else {
            NA
        }
    }

    ## Intra-cluster metrics: list to store distance, time and class for each patient-cluster group
    intra_list <- list()
    cnt <- 0
    for (pt in patients_multi_seq) {
        pt_seqs <- valid_seq_ids[seq2pt == pt]
        pt_clusters <- unique(clusters[pt_seqs])
        for (cl in pt_clusters) {
            clust_seqs <- pt_seqs[clusters[pt_seqs] == cl]
            if (length(clust_seqs) > 1) {
                cnt <- cnt + 1
                res <- compute_max_pair(clust_seqs)
                intra_list[[cnt]] <- list(
                    dist = res["distance"],
                    time = res["time_diff"],
                    class = compute_patient_class(clust_seqs)
                )
            }
        }
    }

    intra_data <- if (cnt > 0) {
        do.call(rbind, lapply(intra_list, as.data.frame))
    } else {
        data.frame(dist = numeric(0), time = numeric(0), class = logical(0))
    }
    intra_data$dist <- as.numeric(as.character(intra_data$dist))
    intra_data$time <- as.numeric(as.character(intra_data$time))

    valid_intra <- !is.na(intra_data$dist) & !is.na(intra_data$time) & !is.na(intra_data$class)
    intra_data <- intra_data[valid_intra, ]

    ## Inter-cluster metrics: compute comparisons for patients in >1 cluster
    inter_list <- list()
    cnt_int <- 0
    for (pt in patients_multi_seq) {
        pt_seqs <- valid_seq_ids[seq2pt == pt]
        pt_clusters <- unique(clusters[pt_seqs])
        if (length(pt_clusters) > 1) {
            pairs <- t(combn(pt_clusters, 2))
            for (i in seq_len(nrow(pairs))) {
                seqs1 <- pt_seqs[clusters[pt_seqs] == pairs[i, 1]]
                seqs2 <- pt_seqs[clusters[pt_seqs] == pairs[i, 2]]
                cnt_int <- cnt_int + 1
                res <- compute_max_pair_inter(seqs1, seqs2)
                inter_list[[cnt_int]] <- list(
                    dist = res["distance"],
                    time = res["time_diff"],
                    class = any(c(seqs1, seqs2) %in% ip_pt_seqs)
                )
            }
        }
    }

    inter_data <- if (cnt_int > 0) {
        do.call(rbind, lapply(inter_list, as.data.frame))
    } else {
        data.frame(dist = numeric(0), time = numeric(0), class = logical(0))
    }
    inter_data$dist <- as.numeric(as.character(inter_data$dist))
    inter_data$time <- as.numeric(as.character(inter_data$time))

    valid_inter <- !is.na(inter_data$time)
    inter_data <- inter_data[valid_inter, ]

    ## Determine common axis limits for plotting
    xlim_max <- max(c(intra_data$dist, inter_data$dist), na.rm = TRUE) + 5
    ylim_max <- max(c(intra_data$time, inter_data$time), na.rm = TRUE) + 5

    ## Plot 1: Scatter plot of genetic distance versus time
    scatter_file <- file.path(
        "figures", paste0(format(Sys.time(), "%Y-%m-%d"), "_", prefix,
                          "_intra_pt_distance_scatter.pdf")
    )
    pdf(scatter_file)
    plot(intra_data$dist[intra_data$class],
        intra_data$time[intra_data$class],
        col = "blue", axes = FALSE, xlab = "", ylab = "",
        xlim = c(0, xlim_max), ylim = c(0, ylim_max)
    )
    par(new = TRUE)
    plot(intra_data$dist[!intra_data$class],
        intra_data$time[!intra_data$class],
        col = "blue", pch = 19, axes = FALSE, xlab = "", ylab = "",
        xlim = c(0, xlim_max), ylim = c(0, ylim_max)
    )
    par(new = TRUE)
    plot(inter_data$dist[inter_data$class],
        inter_data$time[inter_data$class],
        col = "red", axes = FALSE, xlab = "", ylab = "",
        xlim = c(0, xlim_max), ylim = c(0, ylim_max)
    )
    par(new = TRUE)
    plot(inter_data$dist[!inter_data$class],
        inter_data$time[!inter_data$class],
        col = "red", pch = 19,
        xlim = c(0, xlim_max), ylim = c(0, ylim_max),
        xlab = "Distance between pair of isolates",
        ylab = "Time between pair of isolates"
    )
    legend("top",
        legend = c(
            "Intra-cluster index", "Intra-cluster convert",
            "Inter-cluster index", "Inter-cluster convert"
        ),
        col = c("blue", "blue", "red", "red"),
        pch = c(1, 19, 1, 19)
    )
    dev.off()

    ## Plot 2: Boxplot of intra-cluster distances stratified by patient type
    boxplot_file <- file.path(
        "figures", paste0(format(Sys.time(), "%Y-%m-%d"), "_", prefix,
                          "_intra_pt_distance_boxplot.pdf")
    )
    pdf(boxplot_file)
    p_value <- wilcox.test(
        intra_data$dist[intra_data$class],
        intra_data$dist[!intra_data$class]
    )$p.value
    boxplot(intra_data$dist ~ as.numeric(intra_data$class),
        xlim = c(0.5, 2.5), names = c("Convert patients", "Index patients"),
        xlab = "", ylab = "SNV distance between patient isolates"
    )
    points(jitter(as.numeric(intra_data$class) + 1, factor = 1),
        intra_data$dist,
        pch = 19, cex = 0.75
    )
    text(1.5, max(intra_data$dist) - 5, paste0("Wilcox p = ", round(p_value, 4)))
    dev.off()

    ## Plot 3: Boxplot with distances controlled by a linear evolutionary rate model
    boxplot_ctl_file <- file.path(
        "figures", paste0(format(Sys.time(), "%Y-%m-%d"), "_", prefix,
                          "_intra_pt_distance_evoRateCTL_boxplot.pdf")
    )
    pdf(boxplot_ctl_file)
    valid_idx <- !intra_data$class & (intra_data$dist < 20)
    evo_rate <- glm(intra_data$dist[valid_idx] ~ intra_data$time[valid_idx])$coefficients[2]
    intra_ctl <- intra_data$dist - intra_data$time * evo_rate
    p_value_ctl <- wilcox.test(intra_ctl[intra_data$class], intra_ctl[!intra_data$class])$p.value
    boxplot(intra_ctl ~ as.numeric(intra_data$class),
        xlim = c(0.5, 2.5), names = c("Convert patients", "Index patients"),
        xlab = "", ylab = "SNV distance between patient isolates"
    )
    points(jitter(as.numeric(intra_data$class) + 1, factor = 1),
        intra_ctl,
        pch = 19, cex = 0.75
    )
    text(1.5, max(intra_ctl) - 5, paste0("Wilcox p = ", round(p_value_ctl, 4)))
    dev.off()

    ## Plot 4: Boxplot with distances controlled by a fixed evolutionary rate (5 mutations/genome/year)
    boxplot_fixed_file <- file.path(
        "figures", paste0(format(Sys.time(), "%Y-%m-%d"), "_", prefix,
                          "_intra_pt_distance_5mutPerYearCTL_boxplot.pdf")
    )
    pdf(boxplot_fixed_file)
    intra_ctl_fixed <- intra_data$dist - intra_data$time * 0.014
    p_value_fixed <- wilcox.test(intra_ctl_fixed[intra_data$class], intra_ctl_fixed[!intra_data$class])$p.value
    boxplot(intra_ctl_fixed ~ as.numeric(intra_data$class),
        xlim = c(0.5, 2.5), names = c("Convert patients", "Index patients"),
        xlab = "", ylab = "SNV distance between patient isolates"
    )
    points(jitter(as.numeric(intra_data$class) + 1, factor = 1),
        intra_ctl_fixed,
        pch = 19, cex = 0.75
    )
    text(1.5, max(intra_ctl_fixed) - 5, paste0("Wilcox p = ", round(p_value_fixed, 4)))
    dev.off()

    ## Return a matrix with intra-cluster distances and patient classes
    cbind(intra_data$dist, intra_data$class)
}
