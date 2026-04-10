# transclust

The project will involve the implementation and evaluation of algorithms
for detecting transmission clusters in healthcare settings using genomic
data.

> \[!WARNING\] This project is still in the early stages of development
> and the API may change without notice.

## Installation

`transclust` is not available on CRAN yet, but you can install it from
GitHub using `devtools`:

``` r
devtools::install_github("theabhirath/transclust")
```

## Usage

### Epidemiology data

We also provide epidemiological trace data for the isolates in the
dataset. This is used in conjunction with the genomic data to identify
clusters of transmission. The data is provided in an Rdata file, which
can be loaded into R using the
[`load()`](https://rdrr.io/r/base/load.html) function. It has multiple
data objects, including:

- `trace_mat`: A matrix of trace data. Rows are days, columns are
  patients.
- `floor_mat`: A matrix tracking patient locations on the floor of the
  hospital. Rows are days, columns are patients.
- `room_mat`: A matrix tracking the room assignments of patients. Rows
  are days, columns are patients.
- `dates`: A vector of admission dates for each patient.
- `dna_pt_labels`: A vector mapping the sequence IDs for isolates to
  patient IDs.
- `ip_seqs`: A vector of sequence IDs which correspond to intake patient
  sequences presumed to be imported.
- `ip_seqs_3days`: A vector of sequence IDs which correspond to all
  intake patient sequences that tested positive within 3 days of
  admission.

### Example dataset for using the tool

The example dataset is a sequence alignment file that is a smaller
subset of the data generated only from ST258 isolates. This is provided
in the `inst/extdata` directory of the package as `example.FASTA`. The
epidemiological data matching these isolates is provided in the example
dataset as `example.RData`. Dates have been shifted relative to the
start date to protect patient privacy, but the structure of the data is
the same as the complete dataset. The example dataset is used in the
package documentation and vignettes to demonstrate how to use the
package.

This example dataset is not too small, but small enough to still be
manageable for testing and demonstration purposes. It is a subset of the
complete dataset, which is much larger and contains more sequence types.
The example dataset is used to demonstrate the functionality of the
package without overwhelming users with too much data. Locally, when
tested on an M1 MacBook Air, the example dataset only takes about a
minute to run the clustering algorithm and generate the phylogenetic
tree, which is a reasonable amount of time for a demonstration.
