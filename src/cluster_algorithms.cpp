#include <unordered_map>
#include <Rcpp.h>
using namespace Rcpp;

// unsigned 128-bit integer
typedef unsigned __int128 uint128_t;

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

    // Define a getter function for indices from the map
    auto get_index = [&](const std::string &tip) -> int {
        auto it = isolate_index_map.find(tip);
        try {
            if (it == isolate_index_map.end()) throw std::out_of_range("Index not found");
            return it->second;
        } catch (const std::out_of_range &e) {
            Rcpp::stop("Tip label '" + tip + "' not found in isolate_names.");
        }
    };

    // Vector to store the number of defining variants for each subtree
    IntegerVector defining_variants(n_subtrees);

    // Iterate over each subtree
    for (int s = 0; s < n_subtrees; s++) {
        // Current subtree
        List subtree = subtrees[s];
        // Get all tip labels from the subtree
        CharacterVector subtree_tip_labels = subtree["tip.label"];
        int n_subtree_tips = subtree_tip_labels.size();

        // Compute subtree indices
        std::vector<int> subtree_indices;
        std::vector<bool> in_subtree(n_isolates, false);
        subtree_indices.reserve(n_subtree_tips);
        for (int t = 0; t < n_subtree_tips; t++) {
            std::string tip_label = Rcpp::as<std::string>(subtree_tip_labels[t]);
            int index = get_index(tip_label);
            subtree_indices.push_back(index);
            in_subtree[index] = true;
        }
        // Get the indices of isolates not in the subtree
        std::vector<int> outside_indices;
        outside_indices.reserve(n_isolates - subtree_indices.size());
        for (int i = 0; i < n_isolates; i++) {
            if (!in_subtree[i]) outside_indices.push_back(i);
        }

        // Initialize the defining variant count for this subtree
        int dvcount = 0;
        // Iterate over each column in the alignment
        for (int j = 0; j < n_cols; j++) {
            // Check that all subtree isolates have the same base at this column
            char first_base = aln[subtree_indices[0]][j]; 
            if (first_base == 'n') continue; // early continue if ambiguous base
            bool same_in_subtree = true;
            for (int idx : subtree_indices) {
                if (aln[idx][j] != first_base) {
                    same_in_subtree = false;
                    break;
                }
            }
            if (!same_in_subtree) continue; // early continue if bases differ

            // This is an alternative to avoid having to use a boolean set – bitwise operations are faster
            // We use a 128-bit mask to track which bases are seen and assume the bases are standard ASCII characters
            uint128_t seen = 0;
            for (int idx: outside_indices) {
                char base = aln[idx][j];
                if (base == 'n') continue; // ignore 'n'
                int ascii_code = static_cast<int>(base);
                seen |= ((uint128_t)1 << ascii_code); // mark the base as seen
            }
            // Count unique bases by splitting the 128-bit integer into two 64-bit halves
            // and tallying total bits set
            uint64_t low = (uint64_t) seen;
            uint64_t high = (uint64_t) (seen >> 64);
            int unique_count = __builtin_popcountll(low) + __builtin_popcountll(high);

            // Condition: outside must have at least one unique base (ignoring 'n')
            // and the bit corresponding to the first base must not be set
            if (unique_count >= 1 && !((seen >> static_cast<int>(first_base)) & 1)) {
                dvcount++;
            }
        }
        // Store the count of defining variants for this subtree
        defining_variants[s] = dvcount;
    }

    // Return the vector of defining variant counts
    return defining_variants;
}
