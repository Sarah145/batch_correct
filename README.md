## Batch correction pipeline

This repository contains scripts to run a Nextflow pipeline to compare different batch correction methods for single-cell RNA-seq data. This is mostly just a clone of the [BatchBench](https://github.com/cellgeni/batchbench) pipeline from the CellGen IT team at Sanger but I couldn't get that to run so made some edits and added one or two extra things.

The input files for this pipeline must be .Rds files of the uncorrected data as a SingleCellExperiment object (all batches in one object) with batch labels stored in the `batch_key` ('Batch' by default) column and cell type labels stored in the `celltype_key` ('cell_type1' by default) column.



The pipeline will run 7 different batch correction methods on the data:

- Scanorama
- BBKNN
- Seurat 3
- Combat
- Harmony
- limma
- MNNCorrect

For each method, 5 different evaluation metrics are returned:

- Batch entropy (from [BatchBench](https://www.biorxiv.org/content/10.1101/2020.05.22.111211v2)) - measure of how well batches are aligned after correction - related to the probability that for each cell, its *k* nearest neighbors come from a different batch - value reported is average entropy scaled between 0-1 - high batch entropy = well-mixed batches, low batch entropy = poorly-mixed batches.
- Cell type entropy (from [BatchBench](https://www.biorxiv.org/content/10.1101/2020.05.22.111211v2)) - same as batch entropy but using cell type labels instead - high cell type entropy = mixing of cell types (not good), low cell type entropy = cell types are not mixing (good).
- Batch ASW (from [scIB](https://github.com/theislab/scib)) - average silhouette width of batches - scaled between -1-1 - high batch ASW = dense, well-separated batches (bad), low batch ASW = well mixed batches (good).
- Cell type ASW (from [scIB](https://github.com/theislab/scib)) - same as batch ASW but for cell type labels - high cell type ASW = good, low cell type ASW = bad.
- Recovery of marker genes - this idea was taken from the BatchBench paper but couldn't find code for it so wrote my own - not sure if it's right. For methods that correct the expression matrix (Scanorama, Seurat3, Combat, limma, MNNCorrect), found marker genes for each cell type (by batch and in the merged dataset), before and after batch correction, then compared the list of total marker genes identified before batch correction to the list of total marker genes identified after batch correction and calculated the Jaccard similarity index of the two lists. High Jaccard index = gene expression was not distorted too much by batch correction, most markers genes could still be identified (good), low Jaccard index = batch correction highly distorted the gene expression values so not as many marker genes could be recovered (bad). Jaccard index = 1 - all marker genes recovered, Jaccard index = 0 - no marker genes recovered.



To run pipeline:

- Need to have Nextflow and Singularity installed.
- Clone this repo and `cd` into it.
- Pull Singularity image - `singularity pull shub://Sarah145/batch_correct`.
- Edit the `nextflow.config` script with location of data, batch key, cell type key, etc. *Note: profile section of the nextflow.config script in this repo is configured to run on cluster with slurm.*
- Edit `dataset_list.txt` file with name of files - one file on each line, no file extension.
- Run pipeline with `nextflow run main.nf -profile singularity -with-trace trace.txt -with-dag flowchart.png`.
- Compile html report of run by running `./compile_report.R <sample_name>`.

**Overview of pipeline**

<img src="https://github.com/Sarah145/batch_correct/blob/master/imgs/flowchart.png?raw=true">

