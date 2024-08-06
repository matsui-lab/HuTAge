
################
# Preprocessing
################

library(dplyr)
library(data.table)
library(SummarizedExperiment)

# GTEx metadata
sample_tissue <- readRDS("preprocessing/intermediate_file/Global/sample_tissue.rds")

# HUGO gene symbol
hgnc <- readRDS("preprocessing/intermediate_file/Global/hgnc.rds")

# GTEx raw data
rawdat <- fread("data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", data.table = F)
rawdat <- rawdat %>%
  filter(Description %in% hgnc$symbol) %>%
  distinct(Description, .keep_all = TRUE)

gene <- rawdat[, 2]
counts <- rawdat[, -c(1, 2)]
counts <- as.matrix(counts)
rownames(counts) <- gene

idx <- match(colnames(counts), sample_tissue$SAMPID)
colData <- sample_tissue[idx,]

# SummarizedExperiment object
se <- SummarizedExperiment(assays = list(count = counts), colData = colData)
gene <- rownames(assay(se))
# Separate by tissue
tiss <- unique(colData$SMTS)
se_tiss <- list()
for(i in seq_along(tiss)){
  se1 <- se[, se$SMTS == tiss[i]]
  # se1 <- se1[, se1$GROUP %in% c("old", "young")]
  rownames(se1) <- gene
  se_tiss[[i]] <- se1
}
names(se_tiss) <- tiss

tisTb <- readRDS("preprocessing/intermediate_file/Global/tisTb.rds")
tisTb <- as.data.frame((tisTb))
se_tiss <- se_tiss[names(se_tiss) %in% tisTb$GTEx]

saveRDS(se_tiss, "preprocessing/intermediate_file/CellCellInteraction/se_tiss.rds")


############################
# BulkSignalR & CellChat
############################


library(dplyr)
library(data.table)
library(SummarizedExperiment)
library(BulkSignalR)
library(CellChat)
library(doParallel)
library(foreach)

# BULKSIGNALR: Estimating LR correlations and statistical significance (q-values) by group (young, middle, old) -------

# Determine significant LR by group ---> input to CellChat
# Calculate single sample LR activity values based on LR correlations and p-values estimated by generation

# Input data
se_tiss <- readRDS("preprocessing/intermediate_file/CellCellInteraction/se_tiss.rds")

# Parallel model settings
n.proc <- 17
cl <- makeCluster(n.proc)
registerDoParallel(cl)

foreach(i = seq_along(se_tiss), .packages = c("BulkSignalR", "SummarizedExperiment")) %dopar% {
  library(BulkSignalR)
  library(SummarizedExperiment)
  
  se <- se_tiss[[i]]
  age = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
  
  bsrdm2_age <- list()
  bsrinf_age <- list()
  scoresLR_age <- list()
  
  for(j in seq_along(ageclass)){
    
    # STEP1: Extract one age group
    se_age <- se[, se$AGE == age[j]]
    # saveRDS(i, paste0("preprocessing/supplement/logfile/", i, "_", j, ".rds"))
    
    # STEP2: BulkSignalR Estimate significant triples (L-R-pathway) (correlation values and q-values are calculated by group)
    set.seed(123)
    bsrdm <- prepareDataset(counts = assay(se_age))
    bsrdm2 <- learnParameters(bsrdm, plot.folder = "preprocessing/supplement/figures/", filename = "dm")
    bsrinf <- initialInference(bsrdm2, min.cor = 0.25)
    
    # STEP3: BulkSignalR Remove redundancy
    bsrinf.redBP <- reduceToBestPathway(bsrinf)
    
    # STEP4: BulksignalR Calculate activity of L-R interactions (calculated for each single sample)
    bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres = 0.001)
    scoresLR <- scoreLRGeneSignatures(bsrdm2, bsrsig.redBP, name.by.pathway = FALSE)
    
    # Save
    bsrdm2_age[[j]] <- bsrdm2
    bsrinf_age[[j]] <- bsrinf
    scoresLR_age[[j]] <- scoresLR
    names(bsrdm2_age)[j] <- ageclass[j]
    names(bsrinf_age)[j] <- ageclass[j]
    names(scoresLR_age)[j] <- ageclass[j]
  }
  saveRDS(bsrdm2_age, paste0("preprocessing/intermediate_file/CellCellInteraction/bulksigobj_age/bsrdm2_age_", names(se_tiss)[i], ".rds"))
  saveRDS(bsrinf_age, paste0("preprocessing/intermediate_file/CellCellInteraction/bulksigobj_age/bsrinf_age_", names(se_tiss)[i], ".rds"))
  saveRDS(scoresLR_age, paste0("preprocessing/intermediate_file/CellCellInteraction/bulksigobj_age/scoresLR_age_", names(se_tiss)[i], ".rds"))
}
stopCluster(cl)

# CELLCHAT: Input significant LR from each group into CellChat ----------------------------------

# STEP1: Prepare single cell data: Convert single cell data from Tabula Sapiens to Seurat
tisTb <- as.data.frame(readRDS("preprocessing/intermediate_file/Global/tisTb.rds"))
dataDir <- "data/tabulasp"
infiles <- list.files(file.path(dataDir, "h5ad_unzipped"))
infiles <- infiles[grep("h5seurat", infiles)][-1]
tsp_tisname <- gsub("TS_|.h5seurat|.h5ad", "", basename(infiles))
keep <- match(tisTb[, "TSP"], tsp_tisname)
infiles <- infiles[keep]

# Parallel model settings
n.proc <- 17
cl <- makeCluster(n.proc)
registerDoParallel(cl)

# Parallel processing in foreach loop
foreach(i = seq_along(infiles)) %dopar% {
  library(Seurat)
  library(data.table)
  library(rhdf5)
  library(SeuratDisk)
  library(Biobase)
  library(CellChat)
  library(BulkSignalR)
  
  tsp_tiss <- gsub("TS_", "", infiles[i])
  tsp_tiss <- gsub(".h5seurat", "", tsp_tiss)
  gtex_tiss = tisTb$GTEx[tisTb$TSP == tsp_tiss]
  
  # STEP2: single cell data: loading and preprocessing
  sc <- LoadH5Seurat(
    file = file.path(dataDir, "h5ad_unzipped", infiles[i]),
    assays = "RNA",
    reductions = NULL,
    graphs = NULL,
    images = NULL
  )
  
  # Get expression matrix
  mat <- GetAssayData(sc, slot = "counts", assay = "RNA")
  sc$nCount_RNA <- colSums(mat)
  sc$nFeature_RNA <- colSums(mat > 0)
  sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
  sc_norm <- NormalizeData(sc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000) # Normalization
  cellcount <- table(sc_norm@meta.data$cell_ontology_class)
  many <- names(cellcount[cellcount >= 3])
  sc_norm2 <- subset(sc_norm, subset = cell_ontology_class %in% many)
  Idents(sc_norm2) <- "cell_ontology_class"
  
  # Average expression by cell type
  avg <- AverageExpression(sc_norm2, assay = 'RNA', slot = 'data')
  meta <- sc_norm2@meta.data
  
  # CellChat object by tissue
  options(stringsAsFactors = FALSE)
  
  # STEP3: Call significant L-R groups estimated by BulkSignalR
  file_bsrinf = list.files("preprocessing/intermediate_file/CellCellInteraction/bulksigobj_age", pattern = "bsrinf", full.names = TRUE)
  bsrinf_age = readRDS(file_bsrinf[grep(paste0(gtex_tiss, ".rds"), file_bsrinf)])
  lrinter_age = lapply(bsrinf_age, function(x) reduceToBestPathway(x))
  lrinter_age = lapply(bsrinf_age, function(x) LRinter(x))
  lrinter_age = lapply(lrinter_age, function(x) mutate(x, LR = paste0("{", x$L, "} / {", x$R, "}")))
  lr_sig = lapply(lrinter_age, function(x) filter(x, qval < 0.001))
  lr_sig = lapply(lr_sig, function(x) distinct(x, LR, .keep_all = TRUE))
  
  # STEP4: Replace CellChatDB with significant L-R estimated by BulkSignalR
  # User-defined CellChatDB
  CellChatDB <- CellChatDB.human
  
  for(j in seq_along(lr_sig)){
    interaction_input_user <- data.frame(interaction_name = lr_sig[[j]]$LR, ligand = lr_sig[[j]]$L, receptor = lr_sig[[j]]$R, row.names = lr_sig[[j]]$LR)
    geneInfo <- data.frame(Symbol = unique(c(lr_sig[[j]]$L, lr_sig[[j]]$R)))
    db.user <- interaction_input_user
    db.new <- updateCellChatDB(db = db.user, gene_info = geneInfo)
    
    # STEP5: Run CellChat: input single cell data and LR
    cellchat <- createCellChat(object = sc_norm2, group.by = "cell_ontology_class")
    cellchat@idents <- droplevels(cellchat@idents)  
    cellchat@DB <- db.new
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat) # take time
    cellchat <- aggregateNet(cellchat)
    cellchat@LR$LRsig$annotation <- "aging"
    cellchat@LR$LRsig$evidence <- "empty"
    df.net <- subsetCommunication(cellchat)
    saveRDS(cellchat, file = file.path("preprocessing/intermediate_file/CellCellInteraction/cellchatobj", paste(gtex_tiss, names(lr_sig)[j], "cellchat.rds", sep = "_")))
  }
  
  return(NULL)
}
stopCluster(cl)

# CELLCHAT: Obtain LR Probability for app table,
#           Obtain Weight of LR Probability for each Cell->Cell for app heatmap-------

file = list.files("preprocessing/intermediate_file/CellCellInteraction/cellchatobj", pattern = ".rds", full.names = TRUE)
tis = basename(file)
tis = sapply(strsplit(tis, "_"), function(x) x[1])

age = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
tis_list = unique(tis)

cellchat_df_tis <- list()
lr_prob_df_tis <- list()
for(i in seq_along(tis_list)){
  trg_tis = tis_list[i]
  
  # STEP1: Call CellChat object
  cellchat_age <- list()
  lr_prob_age <- list()
  for(j in seq_along(age)){
    trg_file = file[grep(paste0(trg_tis, "_"), file)]
    cellchat = readRDS(trg_file[j])
    
    # STEP2: Obtain LR Probability = strength of cell-cell interaction
    df.net <- subsetCommunication(cellchat)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net$ligand_receptor <- paste(df.net$ligand, df.net$receptor, sep = "_")
    df.net$key <- paste(df.net$source_target, df.net$ligand_receptor, sep = "/")
    lr_prob_age[[j]] <- df.net
    names(lr_prob_age)[j] <- age[j]
    
    # STEP3: Summarize LR Probability by Cell-Cell group and obtain Weight = strength of cell-cell interaction
    df <- df.net %>%
      group_by(source_target) %>%
      summarise(sum_prob = sum(prob, na.rm = TRUE))
    colnames(df) <- c("source_target", age[j]) 
    cellchat_age[[j]] <- df
    names(cellchat_age)[j] <- age[j]
  }
  
  # STEP4: Combine each age data into one tissue data frame
  
  lr_prob_age <- lapply(lr_prob_age, function(x) select(x, key, prob))
  lr_prob_age <- purrr::map2(lr_prob_age, names(lr_prob_age), function(x, y) {
    colnames(x)[2] <- y
    return(x)
  })
  lr_prob_age_bind <- purrr::reduce(lr_prob_age, full_join, by = "key")
  lr_prob_df_tis[[i]] <- lr_prob_age_bind
  names(lr_prob_df_tis)[i] <- tis_list[i]
  
  cellchat_age_bind <- purrr::reduce(cellchat_age, full_join, by = "source_target")
  cellchat_df_tis[[i]] <- cellchat_age_bind
  names(cellchat_df_tis)[i] <- tis_list[i]
  
}

cellchat_df_tis <- lapply(cellchat_df_tis, function(x) separate(x, source_target, into = c("sender", "receiver"), sep = "_", remove = FALSE))
lr_prob_df_tis <- lapply(lr_prob_df_tis, function(x) separate(x, key, into = c("cellcell", "lr"), sep = "/", remove = FALSE))
lr_prob_df_tis <- lapply(lr_prob_df_tis, function(x) separate(x, cellcell, into = c("sender", "receiver"), sep = "_", remove = FALSE))

# Save
saveRDS(cellchat_df_tis, "CellCellInteraction/cellchat_cellcell_heat_tis.rds")
saveRDS(lr_prob_df_tis, "CellCellInteraction/cellchat_lr_table_tis.rds")
