#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to rds input file' 
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "RNA",
    type = 'character',
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-c", "--corrected_assay"),
    action = "store",
    default = "corrected",
    type = 'character',
    help = 'Corrected counts assay name'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Minimum number of cells for a gene to be expressed in.'
  ),
  make_option(
    c("-o", "--output_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(limma))

#args 
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 
batch_key <- opt$batch_key

# read input file
dataset <- readRDS(opt$input_object)

# cell batch label vector
batch_vector <- as.character(dataset[[batch_key]])

# Run limma
assay(dataset, corrected_assay) <- removeBatchEffect(x = assay(dataset, assay_name), batch = batch_vector)

# save corrected object
saveRDS(dataset, file = opt$output_object) 
print("Limma worked!")