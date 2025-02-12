# Software Requirements Specification

## Introduction

This is an R package that aims to provide users with functions to analyze and visualize the transmission clusters of pathogens. It aims to help researchers and epidemiologists understand the spread and evolution of infectious diseases by providing tools for cluster detection, visualization, and statistical analysis. The package includes functions for generating summary plots, performing permutation tests, and analyzing intra-patient genetic variation, among other features.

## Overall Description

**Input**: The package takes as input an alignment of whole sequence genomes of pathogens isolated from patients in a long-term acute care hospital (LTACH). The exact method for generating this alignment is described in the datasets.md file. The alignment is in FASTA format and contains the sequences of the parts of the genome that differ from the reference genome.

**Output**: The package generates a variety of plots and tables that summarize the transmission clusters, including phylogenetic trees, dendograms and heatmaps. It also provides statistical analysis of the clusters, such as permutation tests to determine the significance of the clusters.

### Requirements

The input data format **must** be a FASTA file with the sequences being the parts of the genome that differ from the reference genome. For the software, this package requires the following R packages: `ape`, `phytools`, `phangorn`, `ggplot2`, `heatmap3`, `xlsx`, `igraph`, `RColorBrewer`.