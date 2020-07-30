library(reticulate)
library(Seurat)
library(SingleCellExperiment)
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
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Create seurat object
  seurat_h5ad <- CreateSeuratObject(counts = seurat_counts, meta.data = meta_data)
}

#convert with seurat funciton
h5ad2seurat <- tryCatch({
  ReadH5AD('~/Documents/Project/data/3822_merge.h5ad', assay = "corrected")
}, error=function(e){conv_h5ad2seurat('~/Documents/Project/data/3822_merge.h5ad')})

#save seurat object
saveRDS(h5ad2seurat, file = 'data/3822_merge.Rds')
print("h5ad successfully converted to Seurat object!")


dataset <- readRDS('data/3822_merge.Rds')
assay_name <- 'test'
batch_key <- 'sample'

# convert to seurat object
dataset <- as.Seurat(dataset, counts = assay_name, data = assay_name)
batch_vector <- as.character(dataset$sample)
batch_names  <- names(table(batch_vector))
N_batches <- length(batch_names)

# split seurat object into batches
batch_list <- lapply(1:N_batches, function(x) {abc <- dataset[, dataset[[batch_key]] == batch_names[x]]})
# Normalize and find HVG
for (i in 1:N_batches) {
  #batch_list[[i]] <- NormalizeData(object = batch_list[[i]], 
   #                                verbose = FALSE)
  batch_list[[i]] <- FindVariableFeatures(object = batch_list[[i]], 
                                          selection.method = 'dispersion', 
                                          nfeatures = 10000, 
                                          verbose = F)
}

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors 
if(any(sapply(batch_list, ncol)) < 200) {
  k_filter <- (min(sapply(batch_list, ncol)))
}else{k_filter =200}
# Find integration anchors
n_anchors <- 30
anchors <- FindIntegrationAnchors(object.list = batch_list, anchor.features = 10000, dims = 1:n_anchors, k.filter = k_filter)
# Integrate subsets
corrected_assay <- 'batch_corrected'
integrated <- IntegrateData( new.assay.name = corrected_assay, anchorset = anchors, dims = 1:n_anchors)
# save seurat3 corrected object
saveRDS(integrated, opt$output_object)




# set during IntegrateData
DefaultAssay(integrated) <- "batch_corrected"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
DimPlot(integrated, reduction = "umap", group.by = "sample")
DimPlot(integrated, reduction = "umap", group.by = "cell_type1", label = TRUE, 
              repel = TRUE)
seu <- as.data.frame(Embeddings(integrated[['umap']]))
seu$sample <- integrated[['sample']]$sample
seu$cell_type <- integrated[['cell_type1']]$cell_type1
seu$sample <- factor(seu$sample, levels = c('3822d0', '3822TRT', '3822R'))
ggplot(seu, aes(x = UMAP_1, y = UMAP_2, col=sample)) +
  #geom_polygon(data = hull, aes(col = NULL), alpha = 0, lty = 2, col='black', size=0.1, show.legend = F) +
  geom_point(size = 0.4, alpha = 0.7, show.legend = T) +
  scale_color_manual(values = paletteer_d('unikn::pal_signal')) +
  theme_void(base_size = 22) 
ggplot(seu, aes(x = UMAP_1, y = UMAP_2, col=cell_type)) +
 # geom_polygon(data = hull, aes(col = NULL), alpha = 0, lty = 2, col='black', size=0.1, show.legend = F) +
  geom_point(size = 0.4, alpha=0.8, show.legend = T) +
  scale_color_manual(values = cols) +
  theme_void()
