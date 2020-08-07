#function to build sce
build_sce <- function(counts_mat, meta_data, assay_name, emb_list){
  sce <- SingleCellExperiment(assays = list(corrected = counts_mat),
                              colData = meta_data,
                              rowData = list('feature_names' = rownames(counts_mat)),
                              reducedDims = emb_list)
  sce
}



h5ad2sce <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  #Read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  #Extract Counts   
  #set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  counts_mat <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Extract Embeddings
  emb_names <- reticulate::py_to_r(h5ad_file$obsm_keys())
  emb_list <- lapply(as.list(emb_names), function(x) reticulate::py_to_r(h5ad_file$obsm[x]))
  names(emb_list) <- emb_names
  
  #Build SCE
  sce <- build_sce(counts_mat = counts_mat, meta_data = meta_data, emb_list = emb_list)
  sce
}


sc <- h5ad2sce('data/3822_merge.h5ad')


#save converted object
saveRDS(sc, file = '3822_merge_sce.Rds')
print("h5ad successfully converted to SCE object!")