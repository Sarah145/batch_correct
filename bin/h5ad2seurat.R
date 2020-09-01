#!/usr/bin/env Rscript

#Convert h5ad to seurat object

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad input file' 
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
    help = 'Path to rds output file'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

#suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SingleCellExperiment)) 
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(stringr))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 

conv_h5ad2seurat <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  #read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  #Extract Counts   
  #set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  seurat_counts <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  
  # extract uncorrected counts layer
  uncorrected_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$layers['uncorrected'], index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$layers['uncorrected']$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  seurat_uncorrected_counts <- t(as.matrix(reticulate::py_to_r(uncorrected_counts)))
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Create seurat object
  seurat_h5ad <- CreateSeuratObject(counts = seurat_counts, assay = corrected_assay, meta.data = meta_data)
  seurat_h5ad[[assay_name]] <- CreateAssayObject(counts = seurat_uncorrected_counts)
  return(seurat_h5ad)
}

#read in data
h5ad_path <- opt$input_object

#convert with seurat function
h5ad2seurat <- tryCatch({
  ReadH5AD(h5ad_file, assay = "corrected")
}, error=function(e){conv_h5ad2seurat(h5ad_path)})

#write out rds
saveRDS(h5ad2seurat, file = opt$output_object)
