suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(sva))

#args 
assay_name <- 'corrected'
corrected_assay <- 'batch_corrected' 
batch_key <- 'sample'
# input file
dataset <- readRDS('3822_merge_sce.Rds')
batch_vector <- as.character(dataset[[batch_key]])
# run ComBat
mod_data <- as.data.frame(t(as.matrix(assay(dataset,assay_name))))
mod0 = model.matrix(~ 1, data = mod_data)

assay(dataset, corrected_assay) <- ComBat(
  dat = t(mod_data),
  batch = as.character(dataset$sample),
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)
# save corrected object
saveRDS(dataset, file = 'data/3822_corrected_combat.Rds') 
print("ComBat worked!")
