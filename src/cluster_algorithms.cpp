#include <unordered_map>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeSharedMatrix(CharacterMatrix dna_aln, int out_group) {
    int n_isolates = dna_aln.nrow();
    int n_cols = dna_aln.ncol();
    out_group = out_group - 1; // Convert to 0-based index

    // Pre-convert the dna_aln matrix into a vector of strings (one per isolate)
    std::vector<std::string> aln(n_isolates);
    for (int i = 0; i < n_isolates; i++) {
        std::string row_str = "";
        for (int j = 0; j < n_cols; j++) {
            // Assume each element is a single-character string.
            std::string cell = Rcpp::as<std::string>(dna_aln(i, j));
            row_str.push_back(cell[0]);
        }
        aln[i] = row_str;
    }

    // Precompute validity of each column for the out_group row
    std::vector<bool> out_valid(n_cols);
    for (int k = 0; k < n_cols; k++) {
        char base_out = aln[out_group][k];
        out_valid[k] = (base_out != '-' && base_out != 'n');
    }

    // Matrix to store intermediate shared counts
    NumericMatrix temp(n_isolates, n_isolates);
    int max_shared = 0;

    // Compute shared counts for each pair of isolates (only once per pair)
    for (int i = 0; i < n_isolates; i++) {
        for (int j = i + 1; j < n_isolates; j++) {
            int count = 0;
            for (int k = 0; k < n_cols; k++) {
                char base_i = aln[i][k];
                char base_out = aln[out_group][k];
                char base_j = aln[j][k];
                // Check all conditions at once.
                if (base_i != base_out && base_i == base_j &&
                    base_i != '-' && base_i != 'n' && out_valid[k]) {
                    count++;
                }
            }
            temp(i, j) = count;
            temp(j, i) = count;
            if (count > max_shared) max_shared = count;
        }
    }

    // Create the final shared variant matrix:
    // Each off-diagonal element is max_shared - shared_count, and diagonal is Inf.
    NumericMatrix shared_mat(n_isolates, n_isolates);
    for (int i = 0; i < n_isolates; i++) {
        for (int j = 0; j < n_isolates; j++) {
            if (i == j) {
                shared_mat(i, j) = R_PosInf;
            } else {
                shared_mat(i, j) = max_shared - temp(i, j);
            }
        }
    }

    return shared_mat;
}

// [[Rcpp::export]]
IntegerVector computeDefiningVariants(CharacterMatrix dna_aln, CharacterVector isolate_names, List subtrees) {
    int n_isolates = dna_aln.nrow();
    int n_cols = dna_aln.ncol();
    int n_subtrees = subtrees.size();

    // Pre-convert the dna_aln matrix into a vector of strings (one per isolate)
    std::vector<std::string> aln(n_isolates);
    for (int i = 0; i < n_isolates; i++) {
        std::string row_str = "";
        for (int j = 0; j < n_cols; j++) {
            // Assume each element is a single-character string
            std::string cell = Rcpp::as<std::string>(dna_aln(i, j));
            row_str.push_back(cell[0]);
        }
        aln[i] = row_str;
    }

    // Precompute mapping from isolate names to indices (0-based)
    std::unordered_map<std::string, int> isolate_index_map;
    for (int i = 0; i < n_isolates; i++) {
        isolate_index_map[Rcpp::as<std::string>(isolate_names[i])] = i;
    }

    // Vector to store the number of defining variants for each subtree
    IntegerVector defining_variants(n_subtrees);

    // Iterate over each subtree
    for (int s = 0; s < n_subtrees; s++) {
        // Current subtree
        List subtree = subtrees[s];
        // Get all tip labels from the subtree
        CharacterVector subtree_tip_labels = subtree["tip.label"];
        int n_subtree_tips = subtree_tip_labels.size();

        // Compute subtree indices and mark them in a boolean vector for quick lookup
        std::vector<int> subtree_indices;
        std::vector<bool> in_subtree(n_isolates, false);
        for (int t = 0; t < n_subtree_tips; t++) {
            std::string tip_label = Rcpp::as<std::string>(subtree_tip_labels[t]);
            int index = isolate_index_map[tip_label]; // 0-based index
            subtree_indices.push_back(index);
            in_subtree[index] = true;
        }

        // Get the indices of isolates not in the subtree
        std::vector<int> outside_indices;
        for (int i = 0; i < n_isolates; i++) {
            if (!in_subtree[i]) outside_indices.push_back(i);
        }

        // // Create a vector of indices for all isolates
        // IntegerVector all_indices = seq(0, n_isolates - 1);
        // // Find indices of tip labels in isolate_names
        // IntegerVector subtree_indices = match(subtree_tip_labels, isolate_names);
        // // Now get the indices of the isolates not in the subtree
        // IntegerVector outside_indices = setdiff(all_indices, subtree_indices);

        // Initialize the defining variant count for this subtree
        int dvcount = 0;
        // Iterate over each column in the alignment
        for (int j = 0; j < n_cols; j++) {
            // Check that all subtree isolates have the same base at this column
            char first_base = aln[subtree_indices[0]][j];
            if (first_base == 'n') continue; // early continue if ambiguous base

            bool same_in_subtree = true;
            for (int k : subtree_indices) {
                if (aln[k][j] != first_base) {
                    same_in_subtree = false;
                    break;
                }
            }
            // // Get all bases corresponding to the subtree isolates
            // std::vector<char> subtree_bases;
            // for (int index : subtree_indices) {
            //     char base = aln[index - 1][j]; // -1 for 0-based index
            //     subtree_bases.push_back(base);
            // }

            // // Check that all bases are the same, and then check that the base is not 'n'
            // char first_base = subtree_bases[0];
            // if (first_base == 'n') continue; // early continue if ambiguous base
            // bool same_in_subtree = true;
            // for (int k : subtree_indices) {
            //     if (aln[k - 1][j] != first_base) { // - 1 for 0-based index
            //         same_in_subtree = false;
            //         break;
            //     }
            // }
            if (!same_in_subtree) continue; // early continue if bases differ

            // Use a fixed array to track which bases are seen - we assume the bases are standard ASCII characters
            bool seen[128] = { false };
            int unique_count = 0;
            for (int k : outside_indices) {
                char base = aln[k][j];
                if (base == 'n') continue; // ignore 'n'
                int ascii_code = static_cast<int>(base);
                if (!seen[ascii_code]) {
                    seen[ascii_code] = true;
                    unique_count++;
                }
            }

            // Condition: outside must have more than one unique base (ignoring 'n')
            // and none of those bases should equal first_base.
            if (unique_count > 1 && !seen[static_cast<int>(first_base)]) {
                dvcount++;
            }
 
            // // The outside should have at least more than one unique base ignoring 'n'
            // std::vector<char> outside_bases;
            // for (int index : outside_indices) {
            //     char base = aln[index - 1][j]; // -1 for 0-based index
            //     // Ignore 'n'
            //     if (base != 'n') {
            //         outside_bases.push_back(base);
            //     }
            // }
            // std::set<char> unique_outside_bases(outside_bases.begin(), outside_bases.end());
            // bool outside_different = unique_outside_bases.size() > 1;

            // // Final check: none of the outside bases should match the subtree base
            // bool outside_match = false;
            // for (char base : unique_outside_bases) {
            //     if (base == first_base) {
            //         outside_match = true;
            //         break;
            //     }
            // }

            // // If all conditions are met, increment the defining variant count
            // if (outside_different && !outside_match) {
            //     dvcount++;
            // }
        }
        defining_variants[s] = dvcount;
    }
    return defining_variants;
}
