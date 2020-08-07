suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(limma))

# args 
assay_name <- 'corrected'
corrected_assay <- 'batch_corrected' 
# read input file
dataset <- readRDS('3822_merge_sce.Rds')
# cell batch label vector
batch_vector <- as.character(dataset$sample)
# Run limma
assay(dataset, corrected_assay) <- removeBatchEffect(x = assay(dataset, assay_name), batch = batch_vector)
# save corrected object
saveRDS(dataset, '3822_corrected_limma.Rds')
print("Limma worked!")