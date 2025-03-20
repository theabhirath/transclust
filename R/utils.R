#' Get SNP distance matrix from DNA alignment constructed using a model of DNA evolution.
#'
#' @description
#' Get SNP distance matrix from DNA alignment constructed using a model of DNA evolution.
#'
#' @param dna_aln A DNA alignment object.
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
#' @param dna_aln A DNA alignment object.
#' @param snp_dist A SNP distance matrix.
#' @param method A string indicating the method to use for tree construction. Options are
#'               "nj" (neighbor-joining) or "pars" (maximum parsimony).
#'
#' @returns An object of class `phylo` representing the phylogenetic tree.
#'
#' @importFrom phangorn NJ optim.parsimony as.phyDat
#' @importFrom phytools reroot
#' @export
get_phylo_tree <- function(dna_aln, snp_dist, method = c("nj", "pars")) {
    # Check if the method is valid
    method <- match.arg(method)

    # If method is "nj" or "pars", we need to calculate the out-group: the isolate with the maximum average SNP distance
    out_group <- which.max(rowMeans(snp_dist))

    # For both "nj" and "pars", we need to construct the neighbor-joining tree first
    # and then reroot it to the out-group
    nj_tree <- NJ(snp_dist)
    nj_tree <- reroot(nj_tree, which(nj_tree$tip.label == row.names(snp_dist)[out_group]))

    # If method is "pars", return the maximum parsimony tree
    if (method == "pars") optim.parsimony(nj_tree, as.phyDat(dna_aln), all = TRUE) else nj_tree
}
