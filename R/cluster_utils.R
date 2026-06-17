#' Create a lookup table for isolates along with their cluster assignments.
#'
#' @description
#' This function creates a lookup table for isolates. The columns are:
#' \itemize{
#'   \item cluster: the cluster assignment for the isolate.
#'   \item isolate_id: the isolate ID.
#'   \item patient_id: the patient ID.
#'   \item date: the date of the isolate.
#'   \item adm_pos: a logical value indicating if the isolate is admission positive.
#'   \item prev_surv: the previous surveillance date for the patient.
#'   \item prev_surv_neg: a logical value indicating whether the patient had a prior
#'     surveillance before this isolate and none of those prior surveillances were
#'     positive (i.e. a genuine previous-negative -> current-positive conversion).
#' }
#'
#' @param clusters A numeric vector of cluster assignments.
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param seq2pt A named vector mapping sequence IDs to patient IDs.
#' @param adm_seqs A vector of sequence IDs which correspond to admission positive patient sequences.
#' @param dates A named vector mapping sequence IDs to dates.
#' @param surv_df A data frame with surveillance information for isolates.
#'
#' @returns A data frame with columns: isolate_id, patient_id, date, cluster, adm_pos, prev_surv,
#'   prev_surv_neg.
#'
#' @importFrom stats na.omit
#' @export
get_isolate_lookup <- function(clusters, dna_aln, seq2pt, adm_seqs, dates, surv_df) {
    # core lookup table
    lookup <- data.frame(
        isolate_id = labels(dna_aln),
        patient_id = unname(seq2pt[labels(dna_aln)]),
        date = unname(dates[labels(dna_aln)]),
        cluster = unname(clusters[labels(dna_aln)]),
        adm_pos = labels(dna_aln) %in% adm_seqs
    ) |>
        na.omit()

    # nothing to annotate if no isolates survived the na.omit()
    if (nrow(lookup) == 0L) {
        lookup$prev_surv <- numeric(0)
        lookup$prev_surv_neg <- logical(0)
        return(lookup)
    }

    # pre-split surveillance dates and results by patient_id once
    surv_dates_by_pt <- with(surv_df, split(surv_date, patient_id))
    surv_results_by_pt <- with(surv_df, split(result, patient_id))

    # for each isolate compute:
    # - prev_surv: the latest surveillance date strictly before the isolate (else the
    #   isolate date when surveillance exists only on/after it, else NA).
    # - prev_surv_neg: whether the patient had a prior surveillance and none of the
    #   surveillances before this isolate were positive (result == 1). This marks a
    #   genuine previous-negative -> current-positive conversion; a prior positive
    #   means the patient was already colonized, so it is not a conversion. The date
    #   alone (prev_surv < date) cannot distinguish these, since a prior surveillance
    #   may itself have been positive.
    surv_info <- mapply(
        FUN = function(pid, iso_date) {
            key <- as.character(pid)
            surv_dates <- surv_dates_by_pt[[key]]
            # no surveillance at all for this patient
            if (is.null(surv_dates) || length(surv_dates) == 0L) {
                return(c(prev_surv = NA_real_, prev_surv_neg = NA))
            }
            prior <- surv_dates < iso_date
            if (any(prior)) {
                surv_results <- surv_results_by_pt[[key]]
                c(
                    # latest surveillance strictly before isolate
                    prev_surv = max(surv_dates[prior]),
                    # negative only if no prior surveillance was positive
                    prev_surv_neg = !any(surv_results[prior] == 1, na.rm = TRUE)
                )
            } else {
                # surveillance exists but all on/after isolate date (no prior screen)
                c(prev_surv = iso_date, prev_surv_neg = FALSE)
            }
        },
        pid = lookup$patient_id,
        iso_date = lookup$date
    )

    lookup$prev_surv <- unname(surv_info["prev_surv", ])
    lookup$prev_surv_neg <- as.logical(unname(surv_info["prev_surv_neg", ]))

    lookup
}

#' Get non-single patient clusters from a vector of cluster assignments.
#'
#' @description
#' This function gets non-single patient clusters from a vector of cluster assignments.
#'
#' @param isolate_lookup A lookup table for isolates and their clusters assignments which has
#'                       other relevant epidemiological information. For more information, see
#'                       [get_isolate_lookup()].
#'
#' @returns A numeric vector of cluster assignments with the non-single patient clusters.
#'
#' @export
get_non_single_patient_clusters <- function(isolate_lookup) {
    unique_clusters <- sort(unique(isolate_lookup$cluster))
    cluster_size <- vapply(
        unique_clusters,
        function(x) length(unique(isolate_lookup$patient_id[isolate_lookup$cluster == x])),
        integer(1)
    )
    unique_clusters[cluster_size > 1]
}

#' Remove singleton clusters from a vector of cluster assignments.
#'
#' @description
#' This function removes clusters with only one sequence from a vector of cluster assignments.
#'
#' @param clusters A numeric vector of cluster assignments.
#'
#' @returns A numeric vector of cluster assignments with the singleton clusters removed.
#'
#' @export
remove_singleton_clusters <- function(clusters) {
    # every cluster with only one sequence is a singleton
    singleton_clusters <- which(table(clusters) == 1)
    # remove the singleton clusters
    clusters[!(clusters %in% singleton_clusters)]
}

#' Remap cluster values: each unique value (including each special value) gets its own number
#' in order of appearance.
#'
#' @description
#' This function remaps cluster values so that each unique value (including each special value)
#' gets its own number in order of appearance.
#'
#' @param x A numeric vector of cluster assignments.
#' @param special_val A value to treat as special i.e. it will get its own number in order of appearance.
#'                    This is useful for values that are not part of the cluster assignment, such as
#'                    singleton clusters.
#'
#' @returns A numeric vector of remapped cluster assignments.
#'
#' @keywords internal
remap_cluster_values <- function(x, special_val = 0) {
    lookup <- new.env(hash = TRUE, parent = emptyenv())
    next_id <- 1
    out <- integer(length(x))
    for (i in seq_along(x)) {
        val <- x[i]
        # isTRUE used for NA handling
        if (isTRUE(val == special_val)) {
            out[i] <- next_id # every 'special' gets its own ID
            next_id <- next_id + 1
        } else {
            key <- paste0(val) # NA becomes the string "NA"
            if (exists(key, envir = lookup, inherits = FALSE)) {
                out[i] <- lookup[[key]]
            } else {
                lookup[[key]] <- next_id
                out[i] <- next_id
                next_id <- next_id + 1
            }
        }
    }
    names(out) <- names(x) # keep original names (if any)
    out
}

#' Force cluster assignments to be monophyletic with respect to a tree.
#'
#' @description
#' Checks each multi-isolate cluster for monophyly on `tree` and reconciles those that are not
#' monophyletic using the chosen strategy. The result is renumbered sequentially with
#' [remap_cluster_values()], and each reconciled cluster is reported via [message()] against its
#' final (remapped) cluster id, listing the isolates that make it up.
#'
#' @param clusters A named numeric vector of cluster assignments; the names are the isolate/tip
#'                 labels used to locate members on the tree.
#' @param tree A phylogenetic tree of class `phylo`.
#' @param monophyly_method How to make a non-monophyletic cluster monophyletic. `"expand"` grows the
#'                 cluster to the smallest clade of the tree that contains all of its members,
#'                 absorbing any intervening isolates. `"break_down"` instead splits the cluster into
#'                 the largest monophyletic clades it already contains, without pulling in foreign
#'                 isolates.
#' @param special_val An optional cluster value to leave untouched, e.g. `0` for unclustered isolates
#'                 that should not be treated as a single cluster. Default is `NULL`.
#'
#' @returns A numeric vector of cluster assignments that are monophyletic with respect to the tree,
#'          renumbered sequentially.
#'
#' @importFrom ape subtrees is.monophyletic
#' @keywords internal
enforce_monophyly <- function(
    clusters,
    tree,
    monophyly_method = c("expand", "break_down"),
    special_val = NULL
) {
    monophyly_method <- match.arg(monophyly_method)

    # enumerate every clade (subtree) of the tree once up front
    sub_trees <- subtrees(tree)
    sub_tree_tips <- lapply(sub_trees, function(st) st$tip.label)
    sub_tree_size <- vapply(sub_tree_tips, length, integer(1))

    # cluster ids to consider; a special value (e.g. unclustered isolates) is left untouched as it
    # represents many independent isolates rather than a single cluster
    cluster_ids <- setdiff(unique(clusters), special_val)

    # for each multi-isolate cluster that is not already monophyletic, reconcile it with the tree
    # using the chosen strategy. a single isolate is trivially monophyletic, so only clusters with
    # more than one member need checking. "break_down" assigns fresh ids beyond the current maximum.
    # we record the resulting isolate groups for each reconciled cluster so they can be reported
    # against the final (remapped) cluster ids once all reconciliation is done.
    next_id <- max(clusters) + 1
    reports <- list()
    for (cl in cluster_ids) {
        members <- names(clusters)[clusters == cl]
        # monophyly can only be assessed for members that are tips of the tree
        if (!all(members %in% tree$tip.label)) {
            warning(
                "Cluster ",
                cl,
                " contains isolates absent from the provided tree; skipping monophyly enforcement for it."
            )
            next
        }
        if (length(members) <= 1 || is.monophyletic(tree, members)) {
            next
        }
        if (monophyly_method == "expand") {
            # clades that contain every member of the cluster; the smallest is the monophyletic
            # group these isolates belong to. absorbing its tips into the cluster makes the cluster
            # monophyletic by construction (restricted to isolates we actually clustered, in case
            # the tree carries extra tips such as an outgroup).
            contains_all <- vapply(
                sub_tree_tips,
                function(tips) all(members %in% tips),
                logical(1)
            )
            smallest <- which(contains_all)[which.min(sub_tree_size[contains_all])]
            clade_isolates <- intersect(sub_tree_tips[[smallest]], names(clusters))
            clusters[clade_isolates] <- cl
            # record the cluster's membership before and after expansion so the change is reportable
            reports[[length(reports) + 1]] <- list(
                type = "expand",
                before = members,
                groups = list(clade_isolates)
            )
        } else {
            # "pure" clades: subtrees whose tips all belong to this cluster. these are the
            # monophyletic groups that can be carved out without pulling in foreign isolates.
            is_pure <- vapply(
                sub_tree_tips,
                function(tips) all(tips %in% members),
                logical(1)
            )
            # claim the largest pure clades first; because clades are either nested or disjoint, a
            # smaller clade nested inside an already-claimed one will have no members left to take,
            # so this keeps the cluster split into the fewest, largest monophyletic pieces.
            pure_idx <- which(is_pure)
            pure_idx <- pure_idx[order(sub_tree_size[pure_idx], decreasing = TRUE)]
            unassigned <- members
            # track which isolates land in each new cluster so they can be reported
            new_pieces <- list()
            for (pi in pure_idx) {
                tips <- sub_tree_tips[[pi]]
                if (all(tips %in% unassigned)) {
                    clusters[tips] <- next_id
                    new_pieces[[length(new_pieces) + 1]] <- tips
                    next_id <- next_id + 1
                    unassigned <- setdiff(unassigned, tips)
                }
            }
            # any members not contained in a pure multi-tip clade are monophyletic on their own
            for (m in unassigned) {
                clusters[m] <- next_id
                new_pieces[[length(new_pieces) + 1]] <- m
                next_id <- next_id + 1
            }
            reports[[length(reports) + 1]] <- list(
                type = "break_down",
                groups = new_pieces
            )
        }
    }

    # renumber clusters sequentially; this assigns the final cluster ids that callers see
    clusters <- remap_cluster_values(clusters, special_val)

    # report each reconciliation against the final (remapped) cluster ids. every recorded group of
    # isolates now shares a single final id, looked up from any member of the group.
    if (length(reports) > 0) {
        message(
            "Enforcing monophyly (",
            monophyly_method,
            ") reconciled ",
            length(reports),
            " non-monophyletic cluster(s):"
        )
        for (rep in reports) {
            if (rep$type == "expand") {
                grp <- rep$groups[[1]]
                message(
                    "  cluster ",
                    unname(clusters[grp[1]]),
                    " expanded from {",
                    paste(sort(rep$before), collapse = ", "),
                    "} to {",
                    paste(sort(grp), collapse = ", "),
                    "}"
                )
            } else {
                pieces <- vapply(
                    rep$groups,
                    function(grp) {
                        paste0(
                            unname(clusters[grp[1]]),
                            " {",
                            paste(sort(grp), collapse = ", "),
                            "}"
                        )
                    },
                    character(1)
                )
                message("  cluster broken into: ", paste(pieces, collapse = "; "))
            }
        }
    }

    clusters
}

#' Remove a node from a vector of cluster assignments.
#'
#' @description
#' This function removes a node from a vector of cluster assignments.
#'
#' @param clusters A numeric vector of cluster assignments.
#' @param node The node to remove.
#'
#' @returns A numeric vector of cluster assignments with the node removed.
#'
#' @keywords internal
remove_node_from_clusters <- function(clusters, node) {
    clusters[!(names(clusters) == node)]
}

#' Flatten cluster-patient categorization to a data frame
#'
#' @description
#' Converts the named list output of [cluster_patient_categorization()] into
#' a flat data frame with one row per cluster-patient pair.
#'
#' @param categorization A named list as returned by [cluster_patient_categorization()].
#'
#' @returns A data frame with columns `cluster`, `patient_id`, and `category`.
#'
#' @export
flatten_cluster_patient_categorization <- function(categorization) {
    rows <- lapply(names(categorization), function(cl) {
        cats <- categorization[[cl]]
        data.frame(
            cluster = cl,
            patient_id = names(cats),
            category = unname(cats),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

#' Flatten cluster-overlap categorization to a data frame
#'
#' @description
#' Converts the named vector output of [categorize_cluster_overlap()] into a
#' flat data frame with one row per cluster. The patient-level counterpart is
#' [flatten_cluster_patient_categorization()].
#'
#' @param categorization A named vector as returned by
#'   [categorize_cluster_overlap()], mapping each cluster ID to its overlap
#'   category.
#'
#' @returns A data frame with columns `cluster` and `category`.
#'
#' @export
flatten_cluster_overlap_categorization <- function(categorization) {
    data.frame(
        cluster = names(categorization),
        category = unname(categorization),
        stringsAsFactors = FALSE
    )
}
