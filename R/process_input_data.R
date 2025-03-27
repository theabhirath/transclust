#' Read in Sequence Alignment
#'
#' This function reads in a sequence alignment file and processes the sequence names.
#'
#' @param aln_file A character string specifying the path to the alignment file in FASTA format.
#' @return A DNA object with cleaned sequence names.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads the DNA sequences from the specified FASTA file.
#'   \item Cleans up the sequence names by removing specific substrings and patterns.
#'   \item Replaces certain sequence names with predefined values.
#'   \item Prints the cleaned sequence names.
#'   \item Assigns the cleaned sequence names to the DNA object.
#' }
#' @keywords internal
#' @export
read_in_seq_aln <- function(aln_file) {
    dna_aln <- ape::read.dna(aln_file, format = "fasta")

    # Clean up sequence names
    dna_labels <- labels(dna_aln)
    dna_labels <- gsub("Rush_KPC_|Sample_|001_final_|_R1_|_R_|__|_$", "", dna_labels)
    dna_labels <- gsub("_S.*$", "", dna_labels)
    dna_labels[grepl("KPNIH1", dna_labels)] <- "KPNIH1"

    # # CHECK: Do we really need to replace these labels?
    # replacements <- c("43714" = "1", "43715" = "2", "43716" = "3", "43717" = "4",
    #                   "43718" = "8", "43719" = "14", "43720" = "18", "43721" = "19",
    #                   "43722" = "29", "43723" = "30")
    # dna_labels <- vapply(dna_labels, function(label) {
    #     if (label %in% names(replacements)) replacements[[label]] else label
    # }, character(1))

    # print(dna_labels)
    # # CHECK: Why the dimension check?
    # if (is.null(dim(dna))) {
    #     names(dna) <- dna_labels
    # } else {
    #     row.names(dna) <- dna_labels
    # }

    row.names(dna_aln) <- dna_labels
    dna_aln
}
