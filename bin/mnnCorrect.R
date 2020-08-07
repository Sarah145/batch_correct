library(batchelor)

# args
batch_key <- 'sample'
assay_name <- 'corrected'
k_num <- 30
n_pcs <- 25
corrected_assay <- 'batch_corrected'
svd_dim <-  2

# input file
dataset <- readRDS('3822_merge_sce.Rds')
# batch cell label vector
batch_vector <- as.character(dataset[[batch_key]])
batch_names <- unique(batch_vector)
N_batches <- length(batch_names)
# split object by batches
batch_list <- lapply(1:N_batches, function(x) {assay(dataset[, batch_vector == batch_names[x]], assay_name)})

# run mnnCorrect
corrected <- do.call('mnnCorrect', c(batch_list, c(k=k_num, sigma=sigma, svd.dim=svd_dim)))
# append the corrected expression matrix to object preserving column order
assay(dataset, corrected_assay) <- assay(corrected, "corrected")[, colnames(assay(dataset, assay_name))]
saveRDS(dataset, file = 'data/3822_corrected_mnn.Rds') 
print("mnnCorrect worked!")
