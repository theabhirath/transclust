# hospitraceR

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/theabhirath/hospitraceR/graph/badge.svg)](https://app.codecov.io/gh/theabhirath/hospitraceR)
[![R-CMD-check](https://github.com/theabhirath/hospitraceR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theabhirath/hospitraceR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The project will involve the implementation and evaluation of algorithms for detecting transmission clusters in healthcare settings using genomic data.

> [!WARNING]
> This project is still in the early stages of development and the API may change without notice.

## Installation

`hospitraceR` is not available on CRAN yet, but you can install it from GitHub using `devtools`:

```r
devtools::install_github("theabhirath/hospitraceR")
```

## Usage

### Epidemiology data

We also provide epidemiological trace data for the isolates in the dataset. This is used in conjunction with the genomic data to identify clusters of transmission. The data is provided in an Rdata file, which can be loaded into R using the `load()` function. It has multiple data objects, including:

- `facility_trace`

### Example dataset for using the tool

The example dataset is a sequence alignment file that is a small subset of the data generated only from ST258 isolates. This is provided in the `inst/extdata` directory of the package as `example.FASTA`. The epidemiological data matching these isolates is provided in the example dataset as `example.RData`. Dates have been shifted relative to the start date to protect patient privacy, but the structure of the data is the same as the complete dataset. The example dataset is used in the package documentation and vignettes to demonstrate how to use the package.
