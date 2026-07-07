# hospitraceR

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/theabhirath/hospitraceR/graph/badge.svg)](https://app.codecov.io/gh/theabhirath/hospitraceR)
[![R-CMD-check](https://github.com/theabhirath/hospitraceR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theabhirath/hospitraceR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**hospitraceR** detects transmission clusters of healthcare-associated pathogens by combining
whole-genome sequencing data with patient location traces, surveillance cultures, and admission status.

Beyond clustering, the package can:

- compare two cluster assignments, both structurally and epidemiologically;
- categorize each patient's transmission role and attempt to "explain" transmission scenarios
  in clusters using spatiotemporal overlap information;
- summarize per-cluster genetic distances and acquisition timing;

and more. See the [introductory article](https://theabhirath.github.io/hospitraceR/articles/hospitraceR.html)
for a guided introduction.

Plotting lives in the companion package
[hospitraceRVisualize](https://github.com/theabhirath/hospitraceRVisualize).

## Installation

`hospitraceR` is not on CRAN yet. Install the development version from GitHub with `devtools`:

```r
devtools::install_github("theabhirath/hospitraceR")
```
