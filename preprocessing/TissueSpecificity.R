
###########
# function
###########

getMart2 <- function (x, y, z) {
  # Set URL corresponding to dataset name and version
  version.url <- switch(y,
                        "79" = "https://mar2015.archive.ensembl.org",
                        "80" = "https://may2015.archive.ensembl.org",
                        "81" = "https://july2015.archive.ensembl.org",
                        "82" = "https://sep2015.archive.ensembl.org",
                        "83" = "https://dec2015.archive.ensembl.org",
                        "84" = "https://mar2016.archive.ensembl.org",
                        "85" = "https://july2016.archive.ensembl.org",
                        "86" = "https://oct2016.archive.ensembl.org",
                        "87" = "https://dec2016.archive.ensembl.org",
                        "88" = "https://mar2017.archive.ensembl.org",
                        "89" = "https://may2017.archive.ensembl.org",
                        "90" = "https://aug2017.archive.ensembl.org",
                        "108" = "https://oct2022.archive.ensembl.org"
  )
  
  # Set dataset name corresponding to species
  dataset <- switch(x,
                    "mouse" = "mmusculus_gene_ensembl",
                    "human" = "hsapiens_gene_ensembl",
                    "macaque" = "mmulatta_gene_ensembl",
                    stop("species not supported. choose human, macaque, or mouse.")
  )
  
  mart <- biomaRt::getBM(
    mart = biomaRt::useMart(biomart = "ensembl", dataset = dataset, host = version.url),
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), useCache = FALSE
  )
  
  # Format and merge annotation information
  mart$ensembl_gene_id <- as.character(mart$ensembl_gene_id)
  row.names(mart) <- mart$ensembl_gene_id
  mart <- mart[, c(2:ncol(mart))]
  common_id <- intersect(row.names(mart), row.names(z))
  # common_id <- intersect(mart$external_gene_name, row.names(z))
  z <- z[match(common_id, rownames(z)),]
  # mart <- mart[match(common_id, mart$external_gene_name),]
  mart <- mart[match(common_id, rownames(mart)),]
  # mart <- mart[row.names(mart) %in% row.names(z), ]
  # mart <- mart[order(row.names(mart)), ]
  output <- cbind(mart, z)
  return(output)
}


tau_quantileNormalize <- function (x) 
{
  x.cols <- names(x)
  x.rows <- row.names(x)
  #x[x == 0] <- NA
  x_m <- as.matrix(x)
  #x <- round(preprocessCore::normalize.quantiles(x_m), digits = 3)
  #x <- round(quantile_normalize(x_m), digits = 3)
  x <- round(normalizeQuantiles(x_m),digits = 3)
  x[is.na(x)] <- 0
  x <- data.frame(x)
  names(x)[c(seq_along(x.cols))] <- x.cols
  row.names(x)[c(seq_along(x.rows))] <- x.rows
  return(x)
}


###########################
# Calculation of tau score
###########################

library(tispec)
library(data.table)
library(foreach)
library(doParallel)
library(preprocessCore)

# Use GTEx unnormalized data

gtex_file <- list.files("preprocessing/intermediate_file/TissueSpecificity/unnormalized/", full.names = TRUE)
tissue_list <- gsub(".rds", "", basename(gtex_file[grep(".rds", gtex_file)]))

sample_tissue <- readRDS("preprocessing/intermediate_file/Global/sample_tissue.rds")

gmt <- readRDS("preprocessing/intermediate_file/Global/tsp_gmt.rds")

tis_mean_mat <- vector("list", length(gtex_file))
names(tis_mean_mat) <- gsub(".rds", "", basename(gtex_file))
for(i in seq_along(gtex_file)){
  gexpr <- readRDS(gtex_file[i])
  if(ncol(gexpr) == 0){next}
  age <- sample_tissue$AGE[match(colnames(gexpr), sample_tissue$SAMPID)]
  # Create a dataframe
  df <- data.frame(t(gexpr), age = age, check.names = FALSE)
  df <- split(df, f = df$age)
  mean_mat <- sapply(df, function(x) colMeans(x[, !colnames(x) %in% "age"]))
  colnames(mean_mat) <- names(df)
  tis_mean_mat[[i]] <- mean_mat
  cat(i, "\n")
}

sapply(tis_mean_mat, colnames)
mat1 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, "20-29"]))
mat2 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, "30-39"]))
mat3 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, "40-49"]))
mat4 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, "50-59"]))
mat5 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, "60-69"]))
mat6 <- do.call(cbind, lapply(tis_mean_mat, function(x) x[, colnames(x) %in% "70-79"]))
matlist <- list(mat1, mat2, mat3, mat4, mat5, mat6)
age_level <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")

names(matlist) <- age_level

for(i in seq_along(age_level)){
  saveRDS(matlist[[i]],
          file = paste0("preprocessing/intermediate_file/TissueSpecificity/tissue_mean_matrix/tissue_mean_matrix_", age_level[i], ".rds"))
}

library(limma)
library(edgeR)
library(biomaRt)
library(dplyr)

rdsFiles <- list.files("preprocessing/intermediate_file/TissueSpecificity/tissue_mean_matrix", full.names = TRUE)


age_list <- list()
for(i in seq_along(matlist)){
  mat <- matlist[[i]]
  mat[mat < 1] <- 0
  zero_idx <- rowSums(mat) == 0
  mat <- mat[!zero_idx,]
  gname <- sapply(strsplit(rownames(mat), "\\|"), function(x) x[1])
  gname <- sapply(strsplit(gname, "\\."), function(x) x[1])
  rownames(mat) <- gname
  # filt <- filterByExpr(mat)
  # mat <- mat[rownames(mat) %in% hgnc$symbol,]
  log2Exp <- log2Tran(mat)
  log2Exp <- as.data.frame(log2Exp)
  mat_norm <- data.frame(tau_quantileNormalize(log2Exp), check.names = F)
  
  tauExp <- calcTau(mat_norm)
  tauAnno <- getMart2(x = 'human', y = "108", z = tauExp)
  
  common_geneid <- intersect(rownames(tauExp), rownames(tauAnno))
  tauExp <- tauExp[match(common_geneid, rownames(tauExp)),]
  tauAnno <- tauAnno[match(common_geneid, rownames(tauAnno)),]
  
  mat_norm <- mat_norm[rownames(mat_norm) %in% rownames(tauAnno),]
  
  result <- list(tauExp = tauExp, tauAnno = tauAnno, qnExp = mat_norm)
  age_list[[i]] <- result
  names(age_list)[i] <- names(matlist)[i]

}

# Format the results
aged_score <- vector("list", length(age_list))
names(aged_score) <- names(age_list)

for(i in seq_along(age_list)){
  li <- list() 
  result <- age_list[[i]]
  tauExp <- result$tauExp
  tauAnno <- result$tauAnno
  qnExp <- result$qnExp
  
  score <- sapply(colnames(qnExp), function(x) getTissue(x, qnExp, tauAnno)$score)
  rownames(score) <- rownames(qnExp)
  
  geneid <- cbind(rownames(tauAnno), tauAnno[, 1:2])
  rownames(geneid) <- NULL
  li$geneid <- geneid
  li$score <- score
  li$tauExp <- tauExp
  li$tauAnno <- tauAnno
  aged_score[[i]] <- li
}

common_id <- Reduce(intersect, lapply(aged_score, function(x) rownames(x$tauExp)))
tau_age <- Reduce(cbind, lapply(aged_score, function(x) x$tauExp[match(common_id, rownames(x$tauExp)), 1, drop = F]))
score <- Reduce(cbind, lapply(aged_score, function(x) x$tauExp[match(common_id, rownames(x$score)), 1, drop = F]))
colnames(tau_age) <- colnames(score) <- names(aged_score)
gene_id <- Reduce(rbind, lapply(aged_score, function(x) x$geneid))
gene_id <- unique(gene_id)

common_id <- intersect(rownames(tau_age), intersect(rownames(score), gene_id[, 1]))
tau_age <- tau_age[match(common_id, rownames(tau_age)),]
score <- score[match(common_id, rownames(score)),]
gene_id <- gene_id[match(common_id, gene_id[, 1]),]

tb <- table(gene_id$external_gene_name)
multiple_gene <- names(tb[tb > 1])
idx <- !(gene_id$external_gene_name %in% multiple_gene)
gene_id <- gene_id[idx,]

tau_age <- tau_age[match(gene_id[, 1], rownames(tau_age)),]
score <- score[match(gene_id[, 1], rownames(score)),]
rownames(tau_age) <- rownames(score) <- gene_id$external_gene_name

output <- list(tau_age = tau_age, score = score, gene_id = gene_id)

saveRDS(output, "TissueSpecificity/tau_age.rds")

tauAnno_by_age <- vector("list", length(aged_score))
names(tauAnno_by_age) <- names(aged_score)
for(i in seq_along(aged_score)){
  tauAnno_by_age[[i]] <- aged_score[[i]]$tauAnno
}

saveRDS(tauAnno_by_age, "TissueSpecificity/tauAnno_by_age.rds")