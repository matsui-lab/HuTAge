
#######################################################################
# Preprocessing of GTEx data and creation of GTEx gene expression data
#######################################################################

# GTEx tpm --> normalized ------------------------------------------------------

library(data.table)
library(foreach)
library(doParallel)

# Loading tissue information
sample_tissue <- readRDS("preprocessing/intermediate_file/Global/sample_tissue.rds")

rawdat <- fread("data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", data.table = F)
gname <- paste(rawdat[, 1], rawdat[, 2], sep = "|")
gexpr <- rawdat[, -(1:2)]

unique_tissue <- unique(sample_tissue$SMTS)

# Set the number of cores to use
numCores <- 31  # Use one less to avoid overloading the system
registerDoParallel(cores = numCores)

# Parallel processing using foreach loop
foreach(i = seq_along(unique_tissue)) %dopar% {
  select_sampid <- sample_tissue$SAMPID[sample_tissue$SMTS == unique_tissue[i]]
  bulk.matrix <- as.matrix(gexpr[, colnames(gexpr) %in% select_sampid])
  rownames(bulk.matrix) <- gname
  # bulk.matrix <- normalizeQuantiles(bulk.matrix)
  bulk.matrix <- bulk.matrix[!duplicated(rownames(bulk.matrix)), ]
  saveRDS(bulk.matrix, file = paste0("preprocessing/intermediate_file/CellComposition/normalized/", unique_tissue[i], ".rds"))
}

# Stop the parallel backend
stopImplicitCluster()

# normalized --> normalized_protein coding -------------------------------------

library(viridis)
library(hrbrthemes)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(fgsea)
library(data.table)


deconv_file <- list.files("CellComposition/", full.names = TRUE)
tissue_list <- gsub("_gtex_decomp_all.rda", "", basename(deconv_file[grep("_gtex_decomp_all.rda", deconv_file)]))

sample_tissue <- readRDS("preprocessing/intermediate_file/Global/GTEx_v8_SampleAttributes.rds")

hgnc <- readRDS("preprocessing/intermediate_file/Global/hgnc.rds")

gmt <- readRDS("preprocessing/intermediate_file/Global/tsp_gmt.rds")

for(ti in seq_along(tissue_list)){
  select_tissue <- tissue_list[ti]
  markergene_list <- gmt[grep(paste0("^", select_tissue, "-"), names(gmt))]
  
  # select_tissue
  bulk.matrix <- readRDS(paste0("preprocessing/intermediate_file/normalized/", select_tissue, ".rds"))
  gnm <- sapply(strsplit(rownames(bulk.matrix), "\\|"), function(x) x[2])
  select_idx <- gnm %in% hgnc$symbol[hgnc$locus_group == "protein-coding gene"]
  bulk.matrix <- bulk.matrix[select_idx, ]
  rownames(bulk.matrix) <- gnm[select_idx]
  saveRDS(bulk.matrix, file = paste0("CellComposition/gexpr/", select_tissue, ".rds")) # = normalized_proteincoding
}


#########################
# Deconvolution by bisque
#########################

# Prepare single-cell gene expression from Tabula Sapiens

infiles <- list.files("data/tabulasp/h5ad_unzipped")
infiles <- infiles[grep("h5seurat", infiles)][-1]
tsp_tisname <- gsub("TS_|.h5seurat|.h5ad", "", basename(infiles))

library(doParallel)
# Set up the cluster (specify the number of cores to use)
cl <- makeCluster(length(infiles) * 2)
registerDoParallel(cl)

# Parallel processing in foreach loop
foreach(i = seq_along(infiles)) %dopar% {
  library(Seurat)
  library(data.table)
  library(rhdf5)
  library(SeuratDisk)
  library(Biobase)
  sc <- LoadH5Seurat(
    file = file.path("data/tabulasp/h5ad_unzipped", infiles[i]),
    assays = "RNA",
    reductions = NULL,
    graphs = NULL,
    images = NULL
  )
  
  # Get expression matrix
  mat <- GetAssayData(sc, slot = "counts", assay = "RNA")
  # Save the total expression counts per cell as nCount_RNA
  sc$nCount_RNA <- colSums(mat)
  # Count the number of genes with expression greater than 0 per cell and save it as nFeature_RNA
  sc$nFeature_RNA <- colSums(mat > 0)
  sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
  sc_norm <- NormalizeData(sc, verbose = FALSE, normalization.method = "RC", scale.factor = 1e6) # Normalization
  sc_norm <- FindVariableFeatures(sc_norm, verbose = FALSE, nfeatures = 3000) # Dynamic genes
  var.gene <- VariableFeatures(sc_norm, selection.method = "vst", assay = "RNA")
  all.genes <- rownames(sc_norm)
  sc_norm <- ScaleData(sc_norm, features = all.genes, verbose = FALSE) # Scaling (3 minutes)
  sc_norm <- RunPCA(sc_norm, features = VariableFeatures(object = sc_norm), verbose = FALSE) # PCA (1 minute)
  sc_norm <- RunUMAP(sc_norm, reduction = "pca", dims = 1:30, verbose = FALSE) # UMAP (1 minute)
  
  Idents(sc_norm) <- "cell_ontology_class"
  
  # Exclude cells with less than 3 cells (in line with the method in Bisuque's paper)
  cellcount <- table(sc_norm@meta.data$cell_ontology_class)
  many <- names(cellcount[cellcount >= 3])
  sc_norm2 <- subset(sc_norm, subset = cell_ontology_class %in% many)
  
  library(fgsea)
  
  gmt <- gmtPathways("data/Tabula_Sapiens.txt")
  geneset <- gmt
  # geneset <- gmt[grep("Muscle-", names(gmt))]
  select.gene <- unique(unlist(geneset))
  
  sc.counts.matrix.all <- as.matrix(sc_norm2@assays$RNA@counts)
  
  meta <- sc_norm2@meta.data
  individual.labels <- meta$donor
  cell.type.labels <- meta$cell_ontology_class
  sample.ids <- rownames(meta)
  
  sc.pheno <- data.frame(check.names = F, check.rows = F,
                         stringsAsFactors = F,
                         row.names = sample.ids,
                         SubjectName = individual.labels,
                         cellType = cell.type.labels)
  
  sc.meta <- data.frame(labelDescription = c("SubjectName",
                                             "cellType"),
                        row.names = c("SubjectName",
                                      "cellType"))
  
  sc.pdata <- new("AnnotatedDataFrame",
                  data = sc.pheno,
                  varMetadata = sc.meta)
  
  sc.eset.all <- Biobase::ExpressionSet(assayData = sc.counts.matrix.all,
                                        phenoData = sc.pdata)
  
  saveRDS(sc.eset.all, paste0("preprocessing/intermediate_file/CellComposition/sc_eset_all/", tsp_tisname[i], "_sc.eset.all.rda"))
  return(NULL)
}

# Close the cluster
stopCluster(cl)

# Prepare bulk gene expression GTEx

library(data.table)
# Loading tissue information
sample_tissue <- readRDS("preprocessing/intermediate_file/Global/GTEx_v8_SampleAttributes.rds")

rawdat <- fread("data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", data.table = F)
gname <- rawdat[, 2]
gexpr <- rawdat[, -(1:2)]

tisli <- list(GTEx = sort(unique(sample_tissue$SMTS)), TSP = sort(tsp_tisname))
tisli

tisTb <- readRDS("preprocessing/intermediate_file/Global/tisTb.rds")

data_count <- table(sample_tissue$SMAFRZE, sample_tissue$SMTS)
include_tissue <- colnames(data_count)[data_count["RNASEQ", ] > 20]

tisTb <- tisTb[tisTb[, 1] %in% include_tissue, ]


# Estimate cell composition

library(doParallel)
cl <- makeCluster(nrow(tisTb) * 2)  # Use 4 cores
registerDoParallel(cl)
clusterExport(cl, c("tisTb"))

# Define the function to be parallelized
process_tissue <- function(i) {
  library(BisqueRNA)
  library(Biobase)
  library(data.table)
  
  # Loading tissue information
  sample_tissue <- fread("data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", data.table = F)
  
  rawdat <- fread("data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", data.table = F)
  gname <- rawdat[, 2]
  gexpr <- rawdat[, -(1:2)]
  
  select_sampid <- sample_tissue$SAMPID[sample_tissue$SMTS == tisTb[i, 1]]
  bulk.matrix <- as.matrix(gexpr[, colnames(gexpr) %in% select_sampid])
  rownames(bulk.matrix) <- gname
  bulk.matrix <- bulk.matrix[!duplicated(rownames(bulk.matrix)), ]
  
  sc.eset.all <- readRDS(file.path("preprocessing/intermediate_file/CellComposition/sc_eset_all", paste(tisTb[i, 2], "sc.eset.all.rda", sep = "_")))
  
  if(length(levels(sc.eset.marker$SubjectName)) < 2){
    return(NULL)
  }
  
  test_sample_idx <- 1:ncol(bulk.matrix)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix[, test_sample_idx]) 
  
  res_all <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset.all, markers = NULL, use.overlap = FALSE)
  
  saveRDS(res_all, file = file.path("CellComposition/deconv/", paste(tisTb[i, 1], "gtex_decomp_all.rda", sep = "_")))
}

# Perform parallel processing using foreach loop
foreach(i = 1:nrow(tisTb)) %dopar% {
  process_tissue(i)
}

# Clean up the cluster
stopCluster(cl)

