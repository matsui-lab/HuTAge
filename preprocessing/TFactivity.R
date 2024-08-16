
##########################################
# GTEx order logistic regression analysis
##########################################

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(data.table)
library(dorothea)
library(dplyr)
library(tidyr)
library(tibble)
library(ordinal)

### import dataset ###

sample_tissue <- readRDS("preprocessing/intermediate_file/Global/sample_tissue.rds")

#protein coding gene list（19247 genes）
hgnc_data <- readRDS("preprocessing/intermediate_file/Global/hgnc.rds")
genelist <- hgnc_data$symbol

rawdat <- fread("data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",data.table=F)
gname <- rawdat[,2]
gexpr <- rawdat[,-(1:2)]

idx <- match(colnames(gexpr),sample_tissue$SAMPID)
sample_tissue <- sample_tissue[idx,]
# grp_id <- unique(sample_tissue$SMTSD)
# grp_detail <- sample_tissue$SMTSD
grp_abbr <- sample_tissue$SMTS

tb <- sort(table(grp_abbr))
rmname <- names(tb)[tb < 10]

subgname <- unique(gname[gname %in% genelist])
subexpr <- rawdat[rawdat$Description %in% subgname,]
subexpr <- subexpr[!duplicated(subexpr$Description), ]

rownames(subexpr) <- subexpr$Description
subexpr <- subexpr[,-c(1,2)]

tissue_id <- unique(sample_tissue$SMTS)


# order logis ------------------------------------------------------------------

age_dep_tiss <- vector("list",length(tissue_id))
for(i in seq_along(tissue_id)){
  # select tissue
  tissue <- tissue_id[i]
  cidx <- match(sample_tissue$SAMPID[sample_tissue$SMTS %in% tissue],colnames(subexpr))

  subexpr_tissue <- subexpr[,cidx]
  subexpr_tissue <- subexpr_tissue[apply(subexpr_tissue,1,function(x)mean(x) > 10),]

  # matching sample id to extract age
  idx <- match(colnames(subexpr_tissue),sample_tissue$SAMPID)
  age <- sample_tissue$AGE[idx]
  age <- factor(age,levels = c("20-29","30-39","40-49","50-59","60-69","70-79"),ordered = TRUE)
  sex <- sample_tissue$SEX[idx]
  sex <- factor(sex,levels = c(1,2))
  rin <- sample_tissue$SMRIN[idx]
  model_list <- vector("list",nrow(subexpr_tissue))
  for(j in 1:nrow(subexpr_tissue)){
    #dat <- data.frame(x = unlist(subexpr_tissue[j,]),y=age,row.names = NULL)
    dat <- data.frame(expr = unlist(subexpr_tissue[j,]),sex = sex,rin = rin, age=age,row.names = NULL)
    if(length(unique(sex)) < 2){
      model <- clm(age ~ expr + rin,data = dat)
    }else{
      model <- clm(age ~ expr + sex + rin,data = dat)
    }
    model_list[[j]] <- model
    cat(j,"\n")
  }
  coefmat <- sapply(model_list,function(result)coef(result))
  colnames(coefmat) <- rownames(subexpr_tissue)
  
  ## if only expr pvalue
  pvals <- sapply(model_list,function(result)summary(result)$coefficients["expr", "Pr(>|z|)"])
  qvals <- p.adjust(pvals,method="fdr")
  
  ## if all pvalue
  #pvals <- t(sapply(model_list,function(result)summary(result)$coefficients[c("expr"), "Pr(>|z|)"]))
  #qvals <- apply(pvals,2,function(p)p.adjust(p,method="fdr"))
  #rownames(pvals) <- rownames(qvals) <- rownames(subexpr_tissue)
  tbl <- data.frame(gene=colnames(coefmat),t(coefmat),pvalue=pvals,qvalue=qvals,signif_pval = pvals <= .5,signif_qval = qvals <= .05,check.names = F)
  age_dep_tiss[[i]] <- tbl
  cat(i,"\n")
}
names(age_dep_tiss) <- tissue_id

saveRDS(age_dep_tiss, "TFactivity/age_dep_glob_tiss.rds")


#######################################
# TF activity inference (differential)
#######################################

library(dorothea)
library(tidyr)
library(tibble)
library(ordinal)
library(decouple)
library(dplyr)

# Estimation of transcription factor activity

# Prepare the TF-regulon network
net <- read.table("data/omnipathnet_downloaded_on_20240401.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
net <- as_tibble(net)

# Prepare the DEGs
age_dep_tiss <- readRDS("TFactivity/age_dep_glob_tiss.rds")

# Format DEGs for input
res_act <- list()
for(i in seq_along(age_dep_tiss)){
  input <- age_dep_tiss[[i]]
  input <- as.matrix(input[,"expr", drop = FALSE])
  
  # Run ulm
  contrast_acts <- run_ulm(mat=input, net=net, .source='source', .target='target',
                           .mor='mor', minsize = 5)
  res_act[[i]] <- contrast_acts
  names(res_act)[i] <- names(age_dep_tiss)[i]
}

saveRDS(res_act, "TFactivity/tfact_grpcomp.rds")
saveRDS(res_act, "preprocessing/intermediate_file/TFacticity/tfact_grpcomp.rds")


########################################
# TF activity inference (single sample)
########################################

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(data.table)
library(dorothea)
library(tidyr)
library(tibble)
library(ordinal)
library(decoupleR)
library(plotly)
library(dplyr)
library(flashClust)


# Estimation of transcription factor activity
net <- read.table("data/omnipathnet_downloaded_on_20240401.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
net <- as_tibble(net)

# Run ulm
sample_acts <- run_ulm(mat = as.matrix(subexpr[, 1:2]), net = net, .source = 'source', .target = 'target',
                       .mor = 'mor', minsize = 5)

saveRDS(sample_acts, "preprocessing/intermediate_file/TFacticity/tfact_allsample.rds")


########################
# Heatmap visualization
########################

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(data.table)
library(dorothea)
library(tidyr)
library(tibble)
library(ordinal)
library(decoupleR)
library(plotly)
library(dplyr)
library(flashClust)

# Age dependent TF changes across tissues ----------------------------------------

diftf_tis <- readRDS("preprocessing/intermediate_file/TFacticity/tfact_grpcomp.rds")
head(diftf_tis$Muscle)
dim(diftf_tis$Muscle)

first_df <- data.frame(source = unique(unlist(lapply(diftf_tis, function(x) x$source))))
for(i in seq_along(diftf_tis)){
  diftf <- diftf_tis[[i]]
  diftf <- dplyr::select(diftf, source, score)
  colnames(diftf)[2] <- names(diftf_tis)[i]
  first_df <- left_join(first_df,diftf, by="source")
}

saveRDS(first_df, "TFactivity/agedep_tfact_acrosstis.rds")


####################################
# Single-cell TF activity inference
####################################

# Estimation of TF activity ----------------------------------------------------
# Matching GTEx tissues with TSP tissues
tisTb <- readRDS("preprocessing/intermediate_file/Global/tisTb.rds")
tisTb <- as.data.frame((tisTb))

infiles <- list.files("data/tabulasp/h5ad_unzipped/")
infiles <- infiles[grep("h5seurat",infiles)][-1]
# infiles <- infiles[-14] # Exclude Salivary Gland
tsp_tisname <- gsub("TS_|.h5seurat|.h5ad","",basename(infiles))
keep <- match(tisTb$TSP,tsp_tisname)
infiles <- infiles[keep]
tsp_tisname <- tsp_tisname[keep]
gtex_tisname <- tisTb$GTEx

library(doParallel)
# Set up cluster (specify number of cores to use)
cl <- makeCluster(length(infiles))
registerDoParallel(cl)

# Parallelize the processing with foreach loop
foreach(i = seq_along(infiles)) %dopar% {
  library(Seurat)
  library(data.table)
  library(rhdf5)
  library(SeuratDisk)
  library(Biobase)
  
  # Create Seurat Object
  sc <- LoadH5Seurat(
    file = file.path("data/tabulasp/h5ad_unzipped",infiles[i]),
    assays = "RNA",
    reductions = NULL,
    graphs = NULL,
    images = NULL
  )
  mat1 <- GetAssayData(sc, slot = "counts", assay = "RNA")
  sc$nCount_RNA <- colSums(mat1)
  sc$nFeature_RNA <- colSums(mat1 > 0)
  sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
  sc_norm <- NormalizeData(sc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000) # Normalization
  cellcount <- table(sc_norm@meta.data$cell_ontology_class)
  many <- names(cellcount[cellcount >= 3])
  data <- subset(sc_norm, subset = cell_ontology_class %in% many)
  Idents(data) <- "cell_ontology_class"
  
  
  # Estimation of TF activity
  library(decoupleR)
  library(OmnipathR)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
  
  net <- read.table("data/omnipathnet_downloaded_on_20240401.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  net <- as_tibble(net)
  mat <- as.matrix(data@assays$RNA@data)
  acts <- run_ulm(mat = mat, net = net, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
  
  # Extract ulm and store it in tfsulm in pbmc
  data[['tfsulm']] <- acts %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                       values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  DefaultAssay(object = data) <- "tfsulm"
  
  # Scale the data
  data <- ScaleData(data)
  data@assays$tfsulm@data <- data@assays$tfsulm@scale.data
  saveRDS(data, file = file.path("preprocessing/intermediate_file/TFacticity/sctfact", paste(tsp_tisname[i], "sctfact.rds", sep = "_")))
  saveRDS(acts, file = file.path("preprocessing/intermediate_file/TFacticity/sctfact", paste(tsp_tisname[i], "sctfact_ref.rds", sep = "_")))
  
  return(NULL)
}
stopCluster(cl)

# Remove unnecessary information from the Seurat object, keeping only UMAP data---------

# Load Seurat objects
file1 = list.files("preprocessing/intermediate_file/TFacticity/sctfact/", pattern = "sctfact.rds", full.names = T)
file2 = list.files("preprocessing/intermediate_file/TFacticity/sctfact/", pattern = "sctfact_ref.rds", full.names = T)
tis = basename(file1)
tis = gsub("_sctfact.rds", "", tis)

for(i in seq_along(file1)){
  data = readRDS(file1[i])
  # Load acts data
  acts = readRDS(file2[i])
  umap_data = FetchData(data, vars = c("umap_1", "umap_2", "ident"))
  
  data[['tfsulm']] <- acts %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  Seurat::DefaultAssay(object = data) <- "tfsulm"
  
  # Scale the data
  data <- Seurat::ScaleData(data, do.center = F)
  data@assays$tfsulm@data <- data@assays$tfsulm@scale.data
  
  scaledscore <- as.data.frame(data@assays$tfsulm@data)
  scaledscore <- tibble::rownames_to_column(scaledscore, var = "source")
  
  # Convert to long format using pivot_longer
  scaledscore <- scaledscore %>%
    pivot_longer(
      cols = -source,  # Long format for all columns except the source column
      names_to = "condition",  # New column name
      values_to = "scaledscore"       # New column name for values
    )
  
  acts_scaled <- left_join(acts, scaledscore, by = c("source", "condition"))
  
  umap_data$condition = rownames(umap_data)
  acts_scaled = left_join(acts_scaled, umap_data, by = "condition")
  acts_scaled$statistic = NULL
  acts_scaled$condition = NULL
  acts_scaled$score = NULL
  
  saveRDS(acts_scaled, paste0("/preprocessing/intermediate_file/TFacticity/sctfact/", tis[i], "_sctfact_data.rds"))
}

# Compress file sizes -----------------------------------------------------------

file_list <- list.files("preprocessing/intermediate_file/TFacticity/sctfact", pattern = "data.rds", full.names = TRUE)

for(i in seq_along(file_list)){
  
  data = readRDS(file_list[i])
  
  # Wide format data frame with values from p_value column
  data_wide_pvalue <- data %>%
    select(-scaledscore) %>%
    pivot_wider(
      names_from = source,
      values_from = p_value
    )
  
  # Wide format data frame with values from scaledscore column
  data_wide_scaledscore <- data %>%
    select(-p_value) %>%
    pivot_wider(
      names_from = source,
      values_from = scaledscore
    )
  
  # Combine two data frames into a list
  data_list <- list(
    p_value = data_wide_pvalue,
    scaledscore = data_wide_scaledscore
  )
  
  filename = gsub("data","data_wide",basename(file_list[i]))
  
  saveRDS(data_list, paste0("TFactivity/sctfact/",filename), compress = "xz")
  
}


# speed up
library(fst)
file = list.files("ShinyApp/TFactivity/sctfact", pattern = "sctfact_data_wide.rds", full.names = TRUE)
for (i in seq_along(file)) {
  data <- readRDS(file[i])
  p_value <- data$p_value
  scaledscore <- data$scaledscore
  fst_file <- gsub("rds", "fst", file[i])
  p_value_fst <- gsub("data", "p_value", fst_file)
  scaledscore_fst <- gsub("data", "scaledscore", fst_file)
  write_fst(p_value, p_value_fst)
  write_fst(scaledscore, scaledscore_fst)
}

# compress
library(float)
file = list.files("ShinyApp/TFactivity/sctfact",full.names = TRUE)
tissue = gsub("_sctfact_p_value_wide.fst","",basename(file))
tissue = gsub("_sctfact_scaledscore_wide.fst","",tissue)
tissue = unique(tissue)

for(i in seq_along(tissue)){
  target_tissue = tissue[i]
  target_file = file[grep(target_tissue,file)]
  pvalue_path = target_file[grep("p_value",target_file)]
  pvalue = read_fst(pvalue_path)
  scaledscore_path = target_file[grep("scaledscore",target_file)]
  scaledscore = read_fst(scaledscore_path)
  scaledscore[] <- lapply(scaledscore, function(x) {
    if(is.numeric(x)) {
      as.numeric(fl(x))  
    } else {
      x  
    }
  })
  write_fst(pvalue, pvalue_path, compress = 100)
  write_fst(scaledscore,scaledscore_path, compress = 100)
}
