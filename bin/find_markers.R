#!/usr/bin/env Rscript

#Find markers and compute jaccard similarity

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
    help = 'Name of batch column'
  ),
  make_option(
    c("-t", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Name of cell type column'
  ),
  make_option(
    c("-o", "--output_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to output file'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

#suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SingleCellExperiment)) 
suppressPackageStartupMessages(library(Seurat))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
batch_key <- opt$batch_key
celltype_key <- opt$celltype_key

# read input object
dataset <- readRDS(opt$input_object)

# convert to Seurat
if(class(dataset) == "SingleCellExperiment"){
  dataset1 <- as.Seurat(dataset, counts = assay_name, assay = corrected_assay, data = corrected_assay)
  dataset1[[assay_name]] <- CreateAssayObject(counts = assay(dataset, assay_name))
  dataset <- dataset1
}
  
# Find markers function 
find_markers_seurat <- function(dataset, batch){
  Idents(dataset) <- dataset[[celltype_key]]
  markers_list <- list() #list of dataframes with marker gene info
  for (i in 1:length(levels(dataset))){
    tryCatch({
      print(paste0("** Cell type ", i, "--", levels(dataset)[i]))
      #if computing over the whole dataset, 
      if(batch == F){
        markers_list[["whole_dataset"]][[levels(dataset)[i]]] <- FindMarkers(object = dataset, slot = "data", ident.1 = levels(dataset)[i], ident.2 = NULL, min.pct = 0.5, logfc.threshold = 2)
      }
      if (batch == T){
        markers_list[[levels(dataset)[i]]] <- FindMarkers(object = dataset,  slot = "data", ident.1 = levels(dataset)[i], ident.2 = NULL, min.pct = 0.5,  logfc.threshold = 2)
      }
      # continue if error
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(markers_list)
}

#. Find_markers_by_batch
find_markers_by_batch <- function(dataset){
  table_batches <- table(dataset[[batch_key]])[table(dataset[[batch_key]]) > 0]
  batch_list <- lapply(1:length(table_batches), function(x) dataset[, dataset[[batch_key]]==names(table_batches)[x]])
  names(batch_list) <- names(table_batches)
  by_batch_markers_list <- list()
  #loop through the batches
  for(i in 1:length(batch_list)){
    print(paste0("* Batch number ", i, ": ", names(batch_list)[i]))
    #execute function
    markers_list <- find_markers_seurat(batch_list[[i]], batch = T)
    by_batch_markers_list[[paste0("Batch_", names(batch_list)[i])]] <- markers_list
  }
  return(by_batch_markers_list)
}

# get markers for raw/uncorrected dataset
DefaultAssay(dataset) <- assay_name
raw_merged_markers <- find_markers_seurat(dataset, batch = F)
raw_merged_markers_names <- c()
for(i in 1:length(raw_merged_markers$whole_dataset)){
  raw_merged_markers_names <- c(raw_merged_markers_names, rownames(raw_merged_markers$whole_dataset[[i]]))
}
raw_merged_markers_names <- unique(raw_merged_markers_names)
raw_by_batch_markers <- find_markers_by_batch(dataset)
raw_by_batch_markers_names <- c()
for(i in 1:length(names(raw_by_batch_markers))){
  for(j in 1:length(raw_by_batch_markers[[names(raw_by_batch_markers)[i]]])){
    b <- raw_by_batch_markers[[names(raw_by_batch_markers)[i]]]
    raw_by_batch_markers_names <- c(raw_by_batch_markers_names, rownames(b[[j]]))
  }
}
raw_by_batch_markers_names <- unique(raw_by_batch_markers_names)
raw_total_markers <- unique(raw_merged_markers_names, raw_by_batch_markers_names)


# get markers for corrected dataset
DefaultAssay(dataset) <- corrected_assay
cor_merged_markers <- find_markers_seurat(dataset, batch = F)
cor_merged_markers_names <- c()
for(i in 1:length(cor_merged_markers$whole_dataset)){
  cor_merged_markers_names <- c(cor_merged_markers_names, rownames(cor_merged_markers$whole_dataset[[i]]))
}
cor_merged_markers_names <- unique(cor_merged_markers_names)
cor_by_batch_markers <- find_markers_by_batch(dataset)
cor_by_batch_markers_names <- c()
for(i in 1:length(names(cor_by_batch_markers))){
  for(j in 1:length(cor_by_batch_markers[[names(cor_by_batch_markers)[i]]])){
    b <- cor_by_batch_markers[[names(cor_by_batch_markers)[i]]]
    tryCatch({
      cor_by_batch_markers_names <- c(cor_by_batch_markers_names, rownames(b[[j]]))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
cor_by_batch_markers_names <- unique(cor_by_batch_markers_names)
cor_total_markers <- unique(cor_merged_markers_names, cor_by_batch_markers_names)

# Calculate Jaccard simmilarity index
jacc <- (length(intersect(raw_total_markers, cor_total_markers)))/length(unique(c(raw_total_markers, cor_total_markers)))
write.table(jacc, file = opt$output_file, col.names = 'Jacc_index', row.names = F, quote = F)
