# Code spatial reconstruction enterocytes

This is a collection of analysis scripts that were used in our manuscript "Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis".

## Prerequisites

* The raw and interemediary input data files are hosted on Zenodo: [DOI 10.5281/zenodo.1320734](https://doi.org/10.5281/zenodo.1320734)
* R 3.5
* Matlab R2015_b

## Scripts

* [01_prepare_scrnaseq_data](01_prepare_scrnaseq_data.md) - R: Preparation of the scRNAseq data
* [02_Seurat_processing](02_Seurat_processing.md) - R: Processing of the scRNAseq data
* [03_matlab_zonation_reconstruction](matlab_03_zonation_reconstruction.m) - Matlab: spatial reconstruction
* [04_scRNAseq_plots](04_figures_R.md) - R: tSNE plots
* [05_annotation_files](05_annotation_files.md) - R: Ensembl annotation
* [06_zonation_plots](06_zonation_plots.md) - R: Zonation plots
* [07_rnaseq_germfree_antibiotics](07_rnaseq_germfree_antibiotics.md) - R: Germfree and antibiotic analysis
* [08_pseudotime](08_pseudotime.md) - R: pseudotime analysis



