# Dataset description

## Pre-processing

The sequence alignment that the package takes is generated from the pipeline described in the paper by Hawken et al. (2022)[^1]:

> All whole genomes first undergo QC with FastQC and Trimmomatic, and then raw reads are mapped to sequence-type (ST) specific reference genomes using BWA-MEM algorithm from BWA. Variant calling is performed with samtools. Reference genomes with matching sequence types (STs) to study strains are chosen to maximize core-genome  orthologous regions for SNV identification and decrease erroneous variant calls. Consensus alignments generated from read mapping to ST specific reference genomes undergo detection and masking of putative recombinant regions using gubbins. Lastly, reference genome annotations are generated using Prokka version 1.14.5 and variant annotations are predicted using snpEff version 4.3t.

## Epidemiology data

We also provide epidemiological trace data for the isolates in the dataset. This is used in conjunction with the genomic data to identify clusters of transmission. The data is provided in an Rdata file, which can be loaded into R using the `load()` function. It has multiple data objects, including:

- `floor_mat`: A matrix tracking patient locations on the floor of the hospital. Rows are days, columns are patients.
- `room_mat`: A matrix tracking the room assignments of patients. Rows are days, columns are patients.
- `dates`: A vector of admission dates for each patient.
- `dna_pt_labels`: A vector mapping the sequence IDs for isolates to patient IDs.
- `ip_seqs`: A vector of sequence IDs which correspond to intake patient sequences presumed to be imported.
- `ip_seqs_3days`: A vector of sequence IDs which correspond to all intake patient sequences that tested positive within 3 days of admission.

## Example dataset for using the tool

The example dataset is a sequence alignment file that is a smaller subset of the data generated only from ST258 isolates. This is provided in the `inst/extdata` directory of the package as `example.FASTA`. The epidemiological data matching these isolates is provided in the example dataset as `example.Rdata`. Dates have been shifted relative to the start date to protect patient privacy, but the structure of the data is the same as the complete dataset. The example dataset is used in the package documentation and vignettes to demonstrate how to use the package.

This example dataset is not too small, but small enough to still be manageable for testing and demonstration purposes. It is a subset of the complete dataset, which is much larger and contains more sequence types. The example dataset is used to demonstrate the functionality of the package without overwhelming users with too much data. Locally, when tested on an M1 MacBook Air, the example dataset only takes about a minute to run the clustering algorithm and generate the phylogenetic tree, which is a reasonable amount of time for a demonstration.

## Extended dataset (Real dataset for answering a biological question using the tool.)

The extended dataset is a larger subset of the complete dataset, on which the `extended-vignette` was run. This contains additional FASTA sequences with ST258 isolates and is not provided for the public. The epidemiological data matching these isolates is the same as the example dataset. This is a good dataset to run the tool to see how it performs on a larger dataset because the results still look similar to our example but the analysis is more robust. The extended dataset allows for more comprehensive testing and validation of the tool's capabilities.

The extended dataset allows us to verify that the statistical analyses and visualization functions yield results that are consistent across slightly different types of data. For example, specifically wrt. mutation signatures within clusters, we expect to see similar results across the example and extended datasets at the population level but divergences at the cluster level. This is because all the sequences are ST258 and should average out to similar mutation signatures, but the clusters should still differ significantly enough from the population average to be detectable. This is a good test of the tool's ability to detect clusters and their mutation signatures, as well as its robustness to different datasets.

## Complete dataset

The dataset is comprised of alignment files that have been processed from 462 whole genome isolates of *Klebsiella pneumoniae* from patients at a long-term acute care hospital (LTACH) to only retain the parts of each genome where they differ from the reference. Each FASTA file comprises sequences from one sequence-type of *Klebsiella pneumoniae*. The pre-processing steps to generate these are described earlier.

## Reference

[^1]: Hawken, Shawn E., et al. “Threshold-Free Genomic Cluster Detection to Track Transmission Pathways in Health-Care Settings: A Genomic Epidemiology Analysis.” The Lancet Microbe, vol. 3, no. 9, Sept. 2022, pp. e652–62. DOI.org (Crossref), https://doi.org/10.1016/S2666-5247(22)00115-X.
