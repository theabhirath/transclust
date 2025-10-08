#' Get SNP distance matrix from DNA alignment constructed using a model of DNA evolution.
#'
#' @description
#' Get SNP distance matrix from DNA alignment constructed using a model of DNA evolution.
#'
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param core Logical: if `TRUE`, return the core SNP distance matrix, otherwise
#'             return the full SNP distance matrix.
#'
#' @returns A numeric matrix representing the SNP distance between sequences.
#'
#' @importFrom ape dist.dna
#' @export
get_snp_dist_matrix <- function(dna_aln, core = TRUE) {
    dist.dna(dna_aln, model = "N", pairwise.deletion = !core, as.matrix = TRUE)
}

#' Get phylogenetic tree using neighbor-joining or maximum parsimony method.
#'
#' @description
#' This function constructs a neighbor-joining or maximum parsimony phylogenetic
#' tree from a DNA alignment object and a SNP distance matrix.
#'
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param snp_dist A numeric matrix representing the SNP distance between sequences.
#'                 See [`get_snp_dist_matrix`] for a useful function to generate this.
#' @param method A string indicating the method to use for tree construction. Options are
#'               "nj" (neighbor-joining) or "pars" (maximum parsimony).
#'
#' @returns An object of class `phylo` representing the phylogenetic tree.
#'
#' @importFrom phangorn NJ optim.parsimony as.phyDat midpoint
#' @importFrom ape root
#' @export
get_phylo_tree <- function(dna_aln, snp_dist, method = c("nj", "pars")) {
    # check if the method is valid
    method <- match.arg(method)

    # the out-group is the isolate with the maximum average SNP distance
    out_group <- which.max(rowMeans(snp_dist))

    # construct the neighbor-joining tree
    nj_tree <- NJ(snp_dist)
    # reroot the tree to the out-group
    nj_tree <- root(
        nj_tree,
        which(nj_tree$tip.label == row.names(snp_dist)[out_group]),
        resolve.root = TRUE
    )

    # if method is "pars", construct the maximum parsimony tree
    tree <- if (method == "pars") {
        optim.parsimony(midpoint(nj_tree, node.labels = "support"), as.phyDat(dna_aln), all = TRUE)
    } else {
        nj_tree
    }

    # return the tree
    tree
}

#' Remove singleton clusters from a vector of cluster assignments.
#'
#' @description
#' This function removes clusters with only one sequence from a vector of cluster assignments.
#'
#' @param clusters A numeric vector of cluster assignments.
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
remap_cluster_values <- function(x, special_val = 1) {
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

#' Remove a node from a vector of cluster assignments.
#'
#' @description
#' This function removes a node from a vector of cluster assignments.
#'
#' @param clusters A numeric vector of cluster assignments.
#' @param out_group The node to remove.
#'
#' @returns A numeric vector of cluster assignments with the node removed.
#'
#' @keywords internal
remove_node_from_clusters <- function(clusters, node) {
    clusters[!(names(clusters) == node)]
}
