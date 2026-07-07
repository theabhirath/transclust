#' Compute a pairwise SNP distance matrix from a DNA alignment
#'
#' @description
#' Counts the number of nucleotide differences between every pair of sequences, via
#' [ape::dist.dna()] with model `"N"`.
#'
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param core Logical; if `TRUE`, sites with missing data in any sequence are
#'             dropped before counting, giving a core-genome distance. If `FALSE`, missing data
#'             are handled per pair of sequences (pairwise deletion).
#'
#' @returns A numeric matrix of pairwise SNP distances between sequences.
#'
#' @importFrom ape dist.dna
#' @export
get_snp_dist_matrix <- function(dna_aln, core = TRUE) {
    dist.dna(dna_aln, model = "N", pairwise.deletion = !core, as.matrix = TRUE)
}

#' Build a phylogenetic tree by neighbor-joining or maximum parsimony
#'
#' @description
#' Builds a tree from a DNA alignment and its SNP distance matrix. A neighbor-joining tree is
#' computed first and rooted on the most divergent isolate (the one with the largest mean SNP
#' distance); with `method = "pars"` this is then refined into a maximum-parsimony tree.
#'
#' @param dna_aln A DNA alignment object of class `DNAbin`.
#' @param snp_dist A numeric matrix of SNP distances between sequences. See [get_snp_dist_matrix()].
#' @param method Tree-construction method: `"nj"` (neighbor-joining) or `"pars"` (maximum
#'               parsimony).
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

    # Return the tree
    tree
}
