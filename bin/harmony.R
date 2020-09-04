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
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Batch key in cell metadata'
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "RNA",
    type = 'character',
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-e", "--corrected_emb"),
    action = "store",
    default = "corrected_emb",
    type = 'character',
    help = 'Name of the fastMNN corrected low-dimensional embedding.'
  ),
  make_option(
    c("-n", "--n_pcs"),
    action = "store",
    default = 25,
    type = 'integer',
    help = 'Number of PCs to perform dimensionality reduction'
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

suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(harmony))
suppressPackageStartupMessages(require(Seurat))


# args
assay_name <- opt$assay_name
n_pcs <- opt$n_pcs
batch_key <- opt$batch_key
corrected_emb <- opt$corrected_emb

# read input object
dataset <- readRDS(opt$input_object)
# run PCA
dataset <- runPCA(dataset, exprs_values = assay_name, ncomponents = n_pcs)
seu_dataset <- as.Seurat(dataset, counts = assay_name, data = assay_name, assay = assay_name)

pca <- as.matrix(Embeddings(seu_dataset[['PCA']]))

# cell batch label vector
batch_vector <- as.character(dataset[[batch_key]])
# run Harmony
harmony_embs <- HarmonyMatrix(pca, batch_vector, theta=4, do_pca = F)
seu_dataset[[corrected_emb]] <- CreateDimReducObject(embeddings = harmony_embs, key = "corrected_", assay = DefaultAssay(seu_dataset))
dataset <- as.SingleCellExperiment(seu_dataset, assay = assay_name)
assayNames(dataset) <- c(assay_name, paste0(assay_name, '1'))
reducedDimNames(dataset) <- c('X_pca', 'X_umap', 'PCA', corrected_emb)

# save object with corrected embedding
saveRDS(dataset, opt$output_object)
print("Harmony worked!")
