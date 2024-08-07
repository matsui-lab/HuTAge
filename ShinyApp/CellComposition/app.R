library(shiny)
library(ggplot2)
library(bslib)
library(gridlayout)
library(plotly)
library(DT)
library(shiny)
library(plotly)
library(DT)
library(parallel)
library(flashClust)
library(tidyverse)

#############
no_cores <- 4

deconv_file <- list.files("deconv",full.names = TRUE)
deconv_file <- deconv_file[grep("decomp", deconv_file)]
tissue_list <- gsub("_gtex_decomp_all.rda", "", basename(deconv_file[grep("_gtex_decomp_all.rda", deconv_file)]))
all_deconv_results <- mclapply(deconv_file, readRDS, mc.cores = no_cores)
names(all_deconv_results) <- tissue_list

sample_age <- readRDS("GTEx_v8_SubjectPhenotypes.rds")
sample_tissue <- readRDS("GTEx_v8_SampleAttributes.rds")

#gtex_expr_file <-list.files(file.path(datadir,"expr/normalized_proteincoding"),full.names = TRUE)
gtex_expr_file <- list.files("gexpr",full.names = TRUE)
gtex_expr_file <- gtex_expr_file[gsub(".rds", "", basename(gtex_expr_file)) %in% tissue_list]
gexpr_all <- mclapply(gtex_expr_file, readRDS, mc.cores = no_cores)
names(gexpr_all) <- gsub(".rds", "", basename(gtex_expr_file))

celltype_gmt <- readRDS("tsp_gmt.rds")
celltype_gmt_names_mat <- t(do.call(cbind, strsplit(names(celltype_gmt), "-")))
celltype_gmt_names_mat <- data.frame(1:nrow(celltype_gmt_names_mat), celltype_gmt_names_mat, check.names = F)
colnames(celltype_gmt_names_mat) <- c("id", "tissue", "celltype")

ui <- page_navbar(
  title = "Cell Type Aging",
  selected = "Cell Type Proportion",
  collapsible = TRUE,
  theme = bslib::bs_theme(preset = "spacelab"),
  sidebar = sidebar(
    title = "HuTAge",
    id = "tabset-default-id",
    actionButton(
      inputId = "gotispec",
      label = "Tissue Specificity",
      onclick = "window.open('https://igcore.cloud/GerOmics/HuTAge/TissueSpecificity/', '_self')"
    ),
    actionButton(
      inputId = "gocelltype",
      label = "Cell Type Composition",
      onclick = "window.open('https://igcore.cloud/GerOmics/HuTAge/CellComposition/', '_self')"
    ),
    actionButton(
      inputId = "gotop_tf",
      label = "Transcriptional Factor",
      onclick = "window.open('https://igcore.cloud/GerOmics/HuTAge/TFactivity/', '_self')"
    ),
    actionButton(
      inputId = "gocellcell",
      label = "Cell-Cell Interaction",
      onclick = "window.open('https://igcore.cloud/GerOmics/HuTAge/CellCellInteraction/', '_self')"
    )
  ),
  nav_panel(
    title = "Cell Type Proportion",
    grid_container(
      layout = c(
        "area1 area2"
      ),
      col_sizes = c(
        "0.6fr",
        "1.4fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        card_header(strong("Tissue")),
        card_body(
          selectInput(
            inputId = "select_tissue",
            label = "Select Tissue",
            choices = tissue_list
          ),
          card(
            full_screen = TRUE,
            card_header(strong("t-test")),
            card_body(
              selectInput(
                inputId = "comparison",
                label = "Comparison condition",
                choices = list(
                  "Young > Middle" = "ym_greater",
                  "Young < Middle" = "ym_less",
                  "Young > Elderly" = "ye_greater",
                  "Young < Elderly" = "ye_less",
                  "Middle > Elderly " = "me_greater",
                  "Middle < Elderly" = "me_less"
                ),
                selected = "ye_less"
              ),
              numericInput(
                inputId = "significance",
                label = "Significance Level",
                value = 0.05,
                min = 0.001,
                max = 1,
                step = 0.001
              )
            )
          )
        )
      ),
      grid_card(
        area = "area2",
        card_header(strong("Cell Type Proportion")),
        card_body(
          card(
            card_body(
              plotlyOutput(outputId = "cp_boxplot")
            )
          ),
          card(
            card_body(
              DTOutput(
                outputId = "cell_proportion_table",
                width = "100%"
              )
            )
          )
        )
      )
    )
  ),
  nav_panel(
    title = "Cell Marker Gene",
    grid_container(
      layout = c(
        "area1 area2",
        "area3 area2"
      ),
      row_sizes = c(
        "0.7fr",
        "1.3fr"
      ),
      col_sizes = c(
        "0.66fr",
        "1.34fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        card_header(strong("Tissue of Interest")),
        card_body(
          selectInput(
            inputId = "selectTissue_marker",
            label = "Select Tissue",
            choices = unique(celltype_gmt_names_mat$tissue)
          )
        )
      ),
      grid_card(
        area = "area2",
        card_header(strong("Gene Expression")),
        card_body(
          card(
            card_body(
              plotlyOutput(outputId = "markergene_lineplot")
            )
          ),
          card(
            card_body(
              plotlyOutput(outputId = "gene_boxplot")
            )
          )
        )
      ),
      grid_card(
        area = "area3",
        full_screen = TRUE,
        card_header(strong("Marker Gene Filter")),
        card_body(
          card(
            full_screen = TRUE,
            card_header(strong("Cell Type")),
            card_body(
              selectInput(
                inputId = "selectCelltype_marker",
                label = "Select Cell Type",
                choices = NULL
              )
            )
          ),
          card(
            full_screen = TRUE,
            card_header(strong("Marker Genes")),
            card_body(
              selectizeInput("selected_gene",
                             "Select a Gene",
                             choices = NULL, 
                             options = list('placeholder' = 'Please select a gene'))
            )
          )
        )
      )
    )
  ),
  nav_spacer(),
  nav_panel(
    title = "User Guide",
    tags$iframe(style = "height:100%; width:100%; border:none;", src = "20240729_Hutage_userguide.pdf")
  )
)

# Skip UI and data loading settings

server <- function(input, output, session) {
  
  # Load and process data reactively
  processedData <- reactive({
    select_tissue <- input$select_tissue
    #deconv_result <- readRDS(deconv_file[match(select_tissue, tissue_list)])
    deconv_result <- all_deconv_results[[select_tissue]]
    
    age <- sample_age$AGE[match(sapply(strsplit(colnames(deconv_result$bulk.props), "-"), 
                                       function(x) paste(x[1:2], collapse = "-")), 
                                sample_age$SUBJID)]
    
    young <- age %in% c("20-29", "30-39")
    middle <- age %in% c("40-49", "50-59")
    elderly <- age %in% c("60-69", "70-79")
    
    young_prop <- deconv_result$bulk.props[, young]
    middle_prop <- deconv_result$bulk.props[, middle]
    elderly_prop <- deconv_result$bulk.props[, elderly]
    
    # 
    # young_vs_elderly <- sapply(1:nrow(young_prop),function(ii)t.test(young_prop[ii,],elderly_prop[ii,])$p.value)
    # names(young_vs_elderly) <- rownames(young_prop)
    # target_cell <- names(which(young_vs_elderly <= input$significance))
    
    # Perform t-test based on input$comparison
    if (input$comparison == "ym_greater") {
      test_results <- sapply(1:nrow(young_prop), function(ii) t.test(young_prop[ii,], middle_prop[ii,], alternative = "greater")$p.value)
    } else if (input$comparison == "ym_less") {
      test_results <- sapply(1:nrow(young_prop), function(ii) t.test(young_prop[ii,], middle_prop[ii,], alternative = "less")$p.value)
    } else if (input$comparison == "ye_greater") {
      test_results <- sapply(1:nrow(young_prop), function(ii) t.test(young_prop[ii,], elderly_prop[ii,], alternative = "greater")$p.value)
    } else if (input$comparison == "ye_less") {
      test_results <- sapply(1:nrow(young_prop), function(ii) t.test(young_prop[ii,], elderly_prop[ii,], alternative = "less")$p.value)
    } else if (input$comparison == "me_greater") {
      test_results <- sapply(1:nrow(middle_prop), function(ii) t.test(middle_prop[ii,], elderly_prop[ii,], alternative = "greater")$p.value)
    } else if (input$comparison == "me_less") {
      test_results <- sapply(1:nrow(middle_prop), function(ii) t.test(middle_prop[ii,], elderly_prop[ii,], alternative = "less")$p.value)
    }
    
    names(test_results) <- rownames(young_prop)
    target_cell <- names(which(test_results <= input$significance))
    
    
    dfli <- vector("list", nrow(young_prop))
    for(j in seq_along(target_cell)){
      df_young <- data.frame(prop = young_prop[target_cell[j], ], age = "young")
      df_middle <- data.frame(prop = middle_prop[target_cell[j], ], age = "middle")
      df_elderly <- data.frame(prop = elderly_prop[target_cell[j], ], age = "elderly")
      df <- rbind(df_young, df_middle, df_elderly)
      df$age <- factor(df$age, levels = c("elderly", "middle", "young"), ordered = TRUE)
      df$celltype <- target_cell[j]
      dfli[[j]] <- df
    }
    df_temp <- do.call(rbind, dfli)
    #df_temp$celltype_age <- paste(df_temp$celltype, df_temp$age, sep = "-")
    
    list(df = df_temp, 
         mean_prop = cbind(young = rowMeans(young_prop[target_cell, ]), 
                           middle = rowMeans(middle_prop[target_cell, ]), 
                           elderly = rowMeans(elderly_prop[target_cell, ])))
  })
  
  # Render box plot
  output$cp_boxplot <- renderPlotly({
    data <- processedData()$df
    # data$age <- factor(data$age, levels = c("young", "middle", "elderly"))
    #data <- data %>% arrange(desc(age))
    plot_ly(data = data, x = ~celltype, y = ~prop, color = ~age, type = 'box',
            boxpoints = 'all', jitter = 0.3, pointpos = 0, marker = list(opacity = 0.3)) %>%
      layout(boxmode = "group", yaxis = list(title = 'Proportion'), xaxis = list(title = 'Cell Type', tickangle = -45))
    
  })
  
  # Render mean proportion table
  output$cell_proportion_table <- renderDT({
    mean_prop <- processedData()$mean_prop
    datatable(round(mean_prop, 2), selection = 'single', options = list(scrollX = TRUE))
    
  })
  
  ###############cell type marker gene##############
  observeEvent(input$selectTissue_marker, {
    # Filter celltype choices based on selected tissue
    selected_tissue <- input$selectTissue_marker
    filtered_celltypes <- celltype_gmt_names_mat %>%
      filter(tissue == selected_tissue) %>% `$`('celltype')
    
    # Update choices for selectCelltype_marker
    updateSelectInput(session, "selectCelltype_marker",
                      label = "Select Cell Type",
                      choices = unique(filtered_celltypes))
  })
  
  observe({
    
    req(input$selectTissue_marker, input$selectCelltype_marker)
    
    selected_tissue <- input$selectTissue_marker
    selected_celltype <- input$selectCelltype_marker
    
    gmt_idx <- celltype_gmt_names_mat$id[
      celltype_gmt_names_mat$tissue == selected_tissue &
        celltype_gmt_names_mat$celltype == selected_celltype
    ]
    cellmarker_gene <- unlist(celltype_gmt[gmt_idx])
    
    # Get expression data corresponding to selected tissue
    expr_data <- gexpr_all[[selected_tissue]]
    
    common_genes <- intersect(rownames(expr_data), cellmarker_gene)
    
    # Update the list of gene names
    updateSelectInput(session, "selected_gene", choices = common_genes, selected = NULL)
  })
  
  
  get_expr_data <- function(selected_tissue, selected_celltype, selected_gene = NULL) {
    gmt_idx <- celltype_gmt_names_mat$id[celltype_gmt_names_mat$tissue == selected_tissue & celltype_gmt_names_mat$celltype == selected_celltype]
    tissue_idx <- which(selected_tissue == names(gexpr_all))
    
    if(length(gmt_idx) > 0 & length(tissue_idx) > 0) {
      filter_gene <- celltype_gmt[[gmt_idx]]
      expr_data <- gexpr_all[[tissue_idx]]
      
      # If selected_gene is specified, filter by that gene
      if(!is.null(selected_gene) && selected_gene %in% rownames(expr_data)) {
        expr_data <- expr_data[rownames(expr_data) %in% selected_gene, , drop = FALSE]
      } else {
        expr_data <- expr_data[rownames(expr_data) %in% filter_gene, ]
      }
      
      # Create subsets by age group
      age <- sample_age$AGE[match(sapply(strsplit(colnames(expr_data), "-"), function(x) paste(x[1:2], collapse = "-")), sample_age$SUBJID)]
      age_groups <- sort(unique(age))
      age_expr <- lapply(age_groups, function(ag) expr_data[, age == ag, drop = FALSE])
      names(age_expr) <- age_groups
      
      # Calculate mean by age group
      age_expr_mean <- sapply(age_expr, rowMeans, simplify = "array")
      age_expr_mean <- data.frame(age_expr_mean)
      
      # If vector, convert to data frame
      if(is.vector(age_expr_mean)) {
        age_expr_mean <- data.frame(mean = age_expr_mean)
      }
      age_expr_mean <- tibble::rownames_to_column(age_expr_mean, var = "Gene")
      
      # Return list
      return(list(expr_data = expr_data, 
                  age_expr = age_expr, 
                  age_expr_mean = age_expr_mean, 
                  age_groups = age_groups,
                  age = age))
    } else {
      return(NULL)
    }
  }
  
  
  output$markergene_lineplot <- renderPlotly({
    
    selected_tissue <- input$selectTissue_marker
    selected_celltype <- input$selectCelltype_marker
    
    result <- get_expr_data(selected_tissue, selected_celltype)
    
    if(is.null(result)) return(NULL)
    
    expr_data_mean <- result$age_expr_mean
    colnames(expr_data_mean) <- gsub("X", "", colnames(expr_data_mean))
    colnames(expr_data_mean) <- gsub("\\.", "-", colnames(expr_data_mean))
    expr_data_mean_long <- expr_data_mean %>%
      pivot_longer(cols = -Gene, names_to = "Age", values_to = "Expression")
    
    # Render line plot by gene
    
    expr_data_mean_long %>%
      plot_ly(x = ~Age, y = ~Expression, type = 'scatter', mode = 'lines+markers',
              color = ~Gene, text = ~Gene, hoverinfo = 'text+x+y') %>%
      layout(title = paste(selected_celltype, ": Marker Genes"),
             xaxis = list(title = "Age"),
             yaxis = list(title = "Expression Level"))
    
  })
  
  output$gene_boxplot <- renderPlotly({
    
    req(input$selectTissue_marker, input$selectCelltype_marker, input$selected_gene)
    result <- get_expr_data(input$selectTissue_marker, input$selectCelltype_marker, input$selected_gene)
    
    if(is.null(result)) return(NULL)
    
    # Prepare expression data by age group
    expr_data <- result$expr_data
    age_groups <- result$age_groups
    age <- result$age
    
    plot_data <- data.frame(Age = age, Expression = as.vector(unlist(expr_data)))
    #print(plot_data)
    
    # Render box plot using plotly
    p <- plot_ly(data = plot_data, y = ~Expression, color = ~Age, type = 'box',
                 boxpoints = 'all', jitter = 0.3, pointpos = 0, marker = list(opacity = 0.5)) %>%
      layout(yaxis = list(title = 'Expression Level'),
             xaxis = list(title = 'Age'),
             title = paste("Gene: ", input$selected_gene))
    p
    
  })
  
  
}

shinyApp(ui, server)
