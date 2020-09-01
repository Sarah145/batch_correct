#!/usr/bin/env Rscript

#Convert SCE to h5ad object

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
    c("-o", "--output_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad output file'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

#suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SingleCellExperiment)) 
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(loomR))
suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(stringr))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 

# read input object
sce <- readRDS(opt$input_object)

if(class(sce) == "Seurat"){
  s <- as.SingleCellExperiment(sce, assay = corrected_assay)
  assay(s, assay_name) <- sce@assays[[assay_name]][rownames(sce@assays[[corrected_assay]]), colnames(sce@assays[[corrected_assay]])]
  assayNames(s) <- c(corrected_assay, assay_name)
  sce <- s
}


reducedDimNames(sce) <- str_extract(reducedDimNames(sce), '[^X_]+')

# convert sce2h5ad and save
if(corrected_assay %in% assayNames(sce)){
  main_layer_assay <- corrected_assay
} else{
  main_layer_assay <- assay_name
}

sceasy:::sce2anndata(obj = sce, outFile = opt$output_object, main_layer = main_layer_assay, transfer_layers = assay_name)
print("SCE successfully converted to H5ad object.")
