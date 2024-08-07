# HuTAge

HuTAge is a web-based application designed to address the complexities of aging by providing a comprehensive resource that integrates cross-tissue age-related information. Aging is a multifaceted process involving interorgan and intercellular interactions, and understanding these dynamics requires robust data resources. HuTAge fills a critical gap by combining data from the Genotype-Tissue Expression (GTEx) project and Tabula Sapiens to deliver tissue- and cell-specific aging molecular information in humans.

-   Tissue Specificity Module: Investigates age-dependent changes in tissue specificity of genes, using heatmaps and tau scores to identify tissue-specific genes.
-   Cell Type Composition Module: Analyzes age-dependent changes in cell type proportions within tissues, with visualizations of significant differences between age groups.
-   Transcription Factor Module: Examines the age dependency of transcription factor activities across tissues and visualizes these changes using heatmaps and bar plots.
-   Cell-Cell Interaction Module: Supports the examination of age-dependent changes in cell-cell interaction strength, visualizing interaction probabilities across age groups.

HuTAge is developed using the RStudio R Shiny package and provides an intuitive interface for exploring and visualizing age-related molecular changes across human tissues and cells. The application can be accessed via https, and the source code is provided in this GitHub repository.

In this GitHub repository, executing app.R, located under the ShinyApp directory, requires input data. This input data is generated by running the R script located in the preprocessing directory. The initial input data sources, including GTEx and Tabula Sapiens, used in preprocessing were downloaded using the script below:

-   GTEx

``` bash
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
```

-   Tabula Sapiens

``` bash
wget https://figshare.com/ndownloader/articles/14267219/versions/5
```
