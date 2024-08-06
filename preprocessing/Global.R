# hgnc data --------------------------------------------------------------------
# downloaded on 20240401
library(fread)
hgnc <- fread("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt", data.table = F)
saveRDS("preprocessing/intermediate_file/Global/hgnc.rds")

# GTEx metadata ----------------------------------------------------------------
sample_age <- fread("data/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",data.table = F)
sample_tissue <- fread("data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",data.table = F)

# Modify metadata
subj <- sample_tissue$SAMPID
subj <- sapply(strsplit(subj, "-"), function(x) paste(x[1:2], collapse = "-"))
sample_tissue$SUBJID <- subj 

sample_age <- sample_age %>%
  mutate(GROUP = case_when(
    AGE %in% c("60-69", "70-79") ~ "old",
    AGE %in% c("40-49", "50-59") ~ "middle",
    AGE %in% c("20-29", "30-39") ~ "young")) %>%
  mutate(GROUP = as.factor(GROUP))

sample_tissue <- left_join(sample_tissue, sample_age, by = "SUBJID")

saveRDS(sample_tissue,"preprocessing/intermediate_file/Global/sample_tissue.rds")


# Consensus tissue -------------------------------------------------------------
tisTb <- matrix(c("Adipose Tissue","Fat",
                  "Adrenal Gland",NA,
                  "Bladder","Bladder",
                  "Blood","Blood",
                  "Blood Vessel","Vasculature",
                  "Bone Marrow","Bone_Marrow",
                  "Brain",NA,
                  "Breast",NA,
                  "Cervix Uteri",NA,
                  "Colon","Large_Intestine",
                  "Esophagus",NA,
                  "Fallopian Tube",NA,
                  "Heart","Heart",
                  "Kidney","Kidney",
                  "Liver","Liver",
                  "Lung","Lung",
                  "Muscle","Muscle",
                  "Nerve",NA,
                  "Ovary",NA,
                  "Pancreas","Pancreas",
                  "Pituitary",NA,
                  "Prostate","Prostate",
                  "Salivary Gland","Salivary_Gland",
                  "Skin","Skin",
                  "Small Intestine","Small_Intestine",
                  "Spleen","Spleen",
                  "Stomach",NA,
                  "Testis",NA,
                  "Thyroid",NA,
                  "Uterus","Uterus",
                  "Vagina",NA
),byrow = TRUE,ncol=2,dimnames = list(NULL,c("GTEx","TSP")))
tisTb <- na.omit(tisTb)

saveRDS(tisTb, file = "preprocessing/intermediate_file/Global/tisTb.rds")


# Create gmt file from Tabula Sapiens ------------------------------------------

tisTb <- readRDS("preprocessing/intermediate_file/Global/tisTb.rds")

deconv_file <- list.files("CellComposition", full.names = TRUE, pattern = "rda")
sample_age <- readRDS("preprocessing/intermediate_file/Global/GTEx_v8_SubjectPhenotypes.rds")
tissue_list <- gsub("_gtex_decomp_all.rda", "", basename(deconv_file[grep("_gtex_decomp_all.rda", deconv_file)]))

# Get Tabula_Sapiens.txt from Enrichr library
# "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Tabula_Sapiens"
gmt <- gmtPathways(file.path(datadir, "data/Tabula_Sapiens.txt"))

gmt_tisnm <- sapply(strsplit(names(gmt), "-"), function(x) x[1])
gmt_cellnm <- sapply(strsplit(names(gmt), "-"), function(x) x[2])
common_nm <- intersect(gmt_tisnm, tisTb[, 2])
gmt <- gmt[gmt_tisnm %in% common_nm]
gmt_tisnm <- sapply(strsplit(names(gmt), "-"), function(x) x[1])
gmt_tisnm2 <- tisTb[match(gmt_tisnm, tisTb[, 2]), 1]
gmt_cellnm <- sapply(strsplit(names(gmt), "-"), function(x) x[2])
gmt_newnm <- paste(gmt_tisnm2, gmt_cellnm, sep = "-")
names(gmt) <- gmt_newnm

saveRDS(gmt, file = "CellComposition/tsp_gmt.rds")
saveRDS(gmt, file = "preprocessing/intermediate_file/Global/tsp_gmt.rds")
