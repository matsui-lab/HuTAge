library(shiny)
library(ggplot2)
library(bslib)
library(gridlayout)
library(plotly)
library(DT)
library(data.table)
library(flashClust)
library(shinyHeatmaply)
library(dplyr)
library(tibble)
library(tidyr)
library(ggdendro)
library(heatmaply)
library(shinyWidgets)
library(plotly)
library(RVAideMemoire)
library(metap)
#library(uwot)


tf_agedep <- readRDS("agedep_tfact_acrosstis.rds")
net <- readRDS("omnipathnet.rds")
tfact_grpcomp <- readRDS("tfact_grpcomp.rds")
age_dep_glob_tiss <- readRDS("age_dep_glob_tiss.rds")

file = list.files("sctfact", pattern = "sctfact_data_wide.rds", full.names = TRUE)

library(tidyr)
library(dplyr)

###################################
# Comment out after development #
# file = file[c(1, 8)]
###################################

sctf_data <- lapply(file, function(x) readRDS(x))
names(sctf_data) <- gsub("_", " ", gsub("_sctfact_data_wide.rds", "", basename(file)))

ui <- page_navbar(
  title = "Transcriptional Factor",
  selected = "Tissue",
  collapsible = TRUE,
  theme = bslib::bs_theme(preset = "spacelab"),
  sidebar = sidebar(
    title = "HuTAge",
    id = "tabset-default-id",
    actionButton(
      inputId = "gotop_tf",
      label = "Transcriptional Factor",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/TFact/', '_self')"
    ),
    actionButton(
      inputId = "gotispec",
      label = "Tissue Specificity",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/TisSpec/', '_self')"
    ),
    actionButton(
      inputId = "gocelltype",
      label = "Cell Type Composition",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/CellComp/', '_self')"
    ),
    actionButton(
      inputId = "gocellcell",
      label = "Cell-Cell Interaction",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/CCI/', '_self')"
    )
  ),
  nav_panel(
    title = "Tissue",
    grid_container(
      layout = c(
        "area1 area2"
      ),
      col_sizes = c(
        "0.3fr",
        "1.7fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        card_header(strong("Tissue of interest")),
        card_body(
          selectInput(
            inputId = "select_tissue",
            label = "Select a tissue",
            choices = names(tfact_grpcomp)
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
      ),
      grid_card(
        area = "area2",
        card_body(
          grid_container(
            layout = c(
              "area1 area1 area1",
              "area2 area3 area4"
            ),
            row_sizes = c(
              "1.2fr",
              "0.8fr"
            ),
            col_sizes = c(
              "0.7fr",
              "0.75fr",
              "0.55fr"
            ),
            gap_size = "10px",
            grid_card(
              area = "area1",
              card_header(strong("Aging-Dependent TF Activity")),
              card_body(
                plotlyOutput(outputId = "tf_tis_main_heatmap")
              )
            ),
            grid_card(
              area = "area2",
              card_header(
                textOutput(outputId = "estimation_table_title")
              ),
              card_body(
                DTOutput(outputId = "estimation_table", width = "100%")
              )
            ),
            grid_card(
              area = "area3",
              full_screen = TRUE,
              card_header(
                textOutput(outputId = "tf_tissue_barplot_title")
              ),
              card_body(
                plotlyOutput(outputId = "tf_tissue_barplot")
              )
            ),
            grid_card(
              area = "area4",
              full_screen = TRUE,
              card_header(
                textOutput(outputId = "tf_target_volcano_title")
              ),
              card_body(plotlyOutput(outputId = "tf_target_volcano"))
            )
          )
        )
      )
    )
  ),
  nav_panel(
    title = "Cell",
    grid_container(
      layout = c(
        "area1 area2"
      ),
      col_sizes = c(
        "0.4fr",
        "1.6fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        card_header(strong("Tissue of Interest")),
        card_body(
          tags$div(
            textOutput(outputId = "selected_tissue"),
            style = "font-weight: bold; border: 1px solid #ccc; padding: 5px; margin-bottom: 20px; box-sizing: border-box; border-radius: 4px; max-width: 300px;"
            # textOutput(
            #   outputId = "selected_tissue"
          ),
          radioButtons(
            inputId = "direction",
            label = "Activity Status",
            choices = c("Activation" = "up", "Inactivation" = "down"),
            selected = "up",
            width = "100%"
          ),
          numericInput(
            inputId = "pval",
            label = "Significance Level",
            value = 0.05,
            min = 0.00001,
            max = 1,
            step = 0.000001
          )
        )
      ),
      grid_card(
        area = "area2",
        card_body(
          grid_container(
            layout = c(
              "area1 area1",
              "area2 area3" 
            ),
            row_sizes = c(
              "1.2fr",
              "0.8fr"
            ),
            col_sizes = c(
              "1.0fr",
              "1.0fr"
            ),
            gap_size = "10px",
            grid_card(
              area = "area1",
              card_header(strong("Single-Cell Transcription Factor Activity")),
              card_body(
                plotlyOutput(outputId = "tf_cell_main_umap")
              )
            ),
            grid_card(
              area = "area2",
              card_header(strong("Mean Transcription Factor Activity by Cell Type")),
              card_body(
                DTOutput(outputId = "estimation_table2", width = "100%")
              )
            ),
            grid_card(
              card_header(strong("Cell Type Enrichment")),
              area = "area3",
              card_body(plotlyOutput(outputId = "tf_cell_main_barplot"))
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

server <- function(input, output) {
  
  # -----------------------
  # Tissue
  # -----------------------
  
  filtered_tf_df <- reactiveVal() 
  
  # Reactive expression of tf_agedep filtered based on the selected tissue
  filtered_tf_agedep <- reactive({
    
    req(input$select_tissue) # Ensure tissue is selected
    
    # Retrieve data frame related to the selected tissue
    selected_df <- tfact_grpcomp[[input$select_tissue]]
    
    # Filter rows where p_value is less than or equal to input$significance
    filtered_df <- selected_df %>%
      filter(p_value <= input$significance)
    filtered_df <- filtered_df[,-1]
    colnames(filtered_df) <- c("TF", "Activity", "P-value")
    
    filtered_tf_df(filtered_df)  # Update reactiveVal
    
    # Retrieve only the filtered TF data from tf_agedep
    clean_tf_agedep <- tf_agedep[filtered_df$TF, , drop = FALSE]
    
    clean_tf_agedep <- as.matrix(clean_tf_agedep)
    clean_tf_agedep[is.na(clean_tf_agedep)] <- 0
    clean_tf_agedep[is.nan(clean_tf_agedep)] <- 0
    clean_tf_agedep[is.infinite(clean_tf_agedep)] <- 0
    
    return(clean_tf_agedep)
  })
  
  output$tf_tis_main_heatmap <- renderPlotly({
    req(filtered_tf_agedep()) # Ensure data is available
    
    # Start displaying progress bar
    withProgress(message = 'Calculating heatmap...', value = 0, {
      # Set initial value of progress bar
      setProgress(value = 0.1)
      
      # Retrieve data from filtered_tf_agedep
      mat <- as.matrix(filtered_tf_agedep())
      # Update progress bar
      setProgress(value = 0.5, message = "Preparing heatmap...")
      
      cols <- rep("white", ncol(mat))  # First set all columns to white (or default)
      cols[which(colnames(mat) == input$select_tissue)] <- "brown"  # Set the column corresponding to the selected tissue to brown
      cols <- data.frame(Selected = cols)
      rownames(cols) <- colnames(mat)
      
      p <- heatmaply(
        mat,
        ColSideColors = cols,
        dendrogram = "both",
        hclust_method = "ward",
        scale = "none",
        show_dendrogram = c(FALSE, TRUE),
        showticklabels = c(TRUE, FALSE),
        labRow = NULL,
        labCol = colnames(mat),
        xlab = "Tissues",
        ylab = "Transcription Factors",
        scale_fill_gradient_fun = ggplot2::scale_fill_viridis_b(limits = c(-3, 3))
      )
      
      # Update progress bar to final value
      setProgress(value = 1, message = "Heatmap ready!")
      
      p
    })
  })
  
  # Display data frame filtered based on selected tissue, including only rows where p_value is less than or equal to the specified significance level
  output$estimation_table <- renderDT({
    # Ensure tissue is selected
    req(input$select_tissue)
    
    filtered_tf_df <- filtered_tf_df()
    
    dt <- datatable(filtered_tf_df, options = list(pageLength = 5, scrollX = TRUE), 
                    rownames = FALSE, selection = list(mode = 'single', selected = 1)) %>%
      formatRound(columns = "Activity", digits = 2) %>%
      formatRound(columns = "P-value", digits = 3)  # Correct this
    
    dt
  })
  
  # Render bar plot based on the selected TF
  output$tf_tissue_barplot <- renderPlotly({
    # Do nothing if no row is selected
    req(input$estimation_table_rows_selected)
    
    filtered_df <- filtered_tf_df()  # Retrieve filtered data frame from reactiveVal
    
    selected_rows <- input$estimation_table_rows_selected
    
    # Retrieve the TF corresponding to the selected row
    selected_tf <- filtered_df$TF[selected_rows]
    
    # Retrieve activity data of the selected TF from tf_agedep
    activity_data <- tf_agedep[match(selected_tf, rownames(tf_agedep)), , drop = FALSE]
    
    # Retrieve tissue names
    tissue_names <- colnames(tf_agedep)
    
    # Prepare data for plotting
    plot_data <- data.frame(
      Tissue = tissue_names,
      Activity = as.numeric(activity_data),
      Color = ifelse(tissue_names == input$select_tissue, 'rgb(255,127,14)', 'rgb(158,202,225)')
    )
    
    plot_data <- plot_data %>%
      arrange(desc(Activity))
    
    plot_data$Tissue <- factor(plot_data$Tissue, levels = plot_data$Tissue, ordered = TRUE)
    
    # Render bar plot using Plotly
    p <- plot_ly(data = plot_data, x = ~Tissue, y = ~Activity, type = 'bar',
                 marker = list(color = plot_data$Color, line = list(color = 'rgb(8,48,107)', width = 1.5))) %>%
      layout(title = paste("TF: ", selected_tf),
             xaxis = list(title = "Tissue", tickangle = -45),
             yaxis = list(title = "Activity"))
    return(p)
  })
  
  output$tf_target_volcano <- renderPlotly({
    req(input$estimation_table_rows_selected)
    
    filtered_df <- filtered_tf_df()  # Retrieve filtered data frame from reactiveVal
    selected_rows <- input$estimation_table_rows_selected
    selected_tf <- filtered_df$TF[selected_rows]  # Retrieve the TF corresponding to the selected row
    
    coef_mat <- age_dep_glob_tiss[[input$select_tissue]]
    target_genes <- net$target[net$source == selected_tf]  # Target genes associated with the selected TF
    
    # Prepare data for volcano plot
    plot_data <- coef_mat[rownames(coef_mat) %in% target_genes,]
    plot_data$logp <- -log10(plot_data$pvalue)  # Transform p-values to -log10
    
    # Calculate transparency based on the absolute value of expr
    max_expr <- max(abs(plot_data$expr))  # Retrieve the maximum absolute value of expr
    plot_data$alpha <- abs(plot_data$expr) / max_expr  # Calculate transparency
    
    # Assign colors
    plot_data$color <- ifelse(plot_data$pvalue < 0.05, 
                              ifelse(plot_data$expr < 0, 
                                     paste0('rgba(100, 100, 255, ', plot_data$alpha, ')'), 
                                     paste0('rgba(255, 100, 100, ', plot_data$alpha, ')')), 
                              'rgba(200, 200, 200, 0.5)')  # Always gray with same transparency if p-value is above 0.05
    
    # Create plot using Plotly
    p <- plot_ly(data = plot_data, x = ~expr, y = ~logp, text = rownames(plot_data),
                 mode = 'markers',
                 marker = list(size = 10, color = ~color, line = list(color = ~color, width = 1)),
                 hoverinfo = 'text+x+y') %>%
      layout(xaxis = list(title = "Age effect on gene expression"),
             yaxis = list(title = "-log10(p-value)"),
             plot_bgcolor = "white")
    return(p)
  })
  
  # Reactive value to store the selected TF name
  selected_tf_name <- reactiveVal(NULL)
  
  # Update TF name based on the selected row in the estimation_table
  observeEvent(input$estimation_table_rows_selected, {
    filtered_df <- filtered_tf_df()  # Retrieve filtered data frame
    selected_rows <- input$estimation_table_rows_selected
    # Retrieve the TF name corresponding to the selected row
    if (length(selected_rows) > 0 && nrow(filtered_df) >= selected_rows) {
      selected_tf_name(filtered_df$TF[selected_rows])
    } else {
      selected_tf_name(NULL)  # Set to NULL if no selection
    }
  })
  
  # Title based on the selected tissue name
  output$estimation_table_title <- renderText({
    req(input$select_tissue)
    paste("TF Activity Table of ", input$select_tissue)
  })
  
  output$tf_tissue_barplot_title <- renderText({
    req(input$select_tissue)
    paste("TF activity Across Tissue")
  })
  
  # Text output to display the selected TF name
  output$tf_target_volcano_title <- renderText({
    req(selected_tf_name())  # Ensure TF name is set
    paste(selected_tf_name(), "Target DEGs")
  })
  
  # --------------------
  # Cell
  # --------------------
  
  tf_sub <- reactive({
    # Change filtering conditions based on user input
    if(input$tf_gene_input_checkbox && !is.null(input$tf_gene_input_gene) && input$tf_gene_input_gene != "") {
      net[net$target %in% input$tf_gene_input_gene, "source", drop = TRUE]
    } else {
      diftf <- tfact_grpcomp[[input$tf_tis_input_tis]]  # Clarify the use of the tissue variable, either input or another reactive value
      diftf[diftf$p_value < 0.2, "source", drop = TRUE]
    }
  })
  
  # Selected tissue name
  output$selected_tissue <- renderText({
    req(input$select_tissue)
    print(input$select_tissue)
  })
  
  filtered_sctf_df <- reactive({
    data <- sctf_data[[input$select_tissue]]
    selected_rows <- input$estimation_table_rows_selected
    filtered_df <- filtered_tf_df()
    selected_tf <- filtered_df$TF[selected_rows]
    # data <- filter(data, source == selected_tf)
    data$scaledscore <- select(data$scaledscore, umap_1, umap_2, ident, scaledscore = selected_tf)
    data$p_value <- select(data$p_value, umap_1, umap_2, ident, p_value = selected_tf)
    data <- left_join(data$scaledscore, data$p_value, by = c("umap_1", "umap_2", "ident"))
    return(data)
  })
  
  output$tf_cell_main_umap <- renderPlotly({
    req(c(input$select_tissue, input$estimation_table_rows_selected)) # Ensure data is available
    
    # Start displaying progress bar
    withProgress(message = 'Downloading data...', value = 0, {
      # Set initial value of progress bar
      setProgress(value = 0.1)
      
      data <- filtered_sctf_df()
      
      setProgress(value = 0.5, message = "Preparing umap...")
      
      # filter
      if (input$direction == "up") {
        data <- data %>%
          mutate(scaledscore = ifelse(scaledscore > 0 & p_value < input$pval, scaledscore, 0))
      } else {
        data <- data %>%
          mutate(scaledscore = ifelse(scaledscore < 0 & p_value < input$pval, scaledscore, 0))
      }
      
      # Generate UMAP plot
      p <- plot_ly(data,
                   x = ~umap_1,
                   y = ~umap_2, 
                   text = ~ident,
                   type = 'scatter', 
                   mode = 'markers',
                   showlegend = FALSE,
                   marker = list(
                     size = 3,
                     color = ~scaledscore,
                     colorscale = list(c(0, "rgb(0,0,255)"),  # Blue - minimum
                                       c(0.49, "rgb(240,240,240)"), # Very light gray, almost white
                                       c(0.51, "rgb(240,240,240)"), # Keep around the mid-value
                                       c(1, "rgb(255,0,0)")), 
                     cmax = 2,
                     cmin = -2,
                     colorbar = list(
                       yanchor = "middle",
                       y = 0.5,
                       ypad = 20,
                       thickness = 20,
                       len = 0.5
                     )
                   )
      ) %>%
        layout(title = paste("TF: ", selected_tf_name()), 
               xaxis = list(title = 'UMAP 1', zeroline = FALSE), 
               yaxis = list(title = 'UMAP 2', zeroline = FALSE))
      
      setProgress(value = 1, message = "UMAP ready!")
      
      return(p)
      
    })
  })
  
  output$tf_cell_main_barplot <- renderPlotly({
    data <- filtered_sctf_df()
    
    if (input$direction == "up") {
      counts <- data %>%
        group_by(ident) %>%
        summarise(
          positive_count = sum(scaledscore > 0 & p_value < input$pval),
          negative_count = n() - positive_count
        )
    } else {
      counts <- data %>%
        group_by(ident) %>%
        summarise(
          positive_count = sum(scaledscore < 0 & p_value < input$pval),
          negative_count = n() - positive_count
        )
    }
    
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$ident
    counts$ident <- NULL
    
    # Fisher multiple testing
    fisher_res <- fisher.multcomp(as.matrix(counts), p.method = "BH")
    tri_mat = fisher_res$p.value
    
    melt_mat = melt(tri_mat)
    cell = unique(c(melt_mat$Var1, melt_mat$Var2))
    
    # Combine p-values with corresponding cell type for each cell
    fisher_df <- data.frame(cell = cell, log.pval = numeric(length(cell)))
    for(i in seq_along(cell)){
      trgcell <- cell[i]
      melt_mat$Var1 <- as.character(melt_mat$Var1)
      melt_mat$Var2 <- as.character(melt_mat$Var2)
      filtmat <- melt_mat[melt_mat$Var1 == trgcell | melt_mat$Var2 == trgcell,]
      filtmat <- drop_na(filtmat)
      temp <- sumlog(filtmat$value, log.p = TRUE)
      fisher_df$log.pval[i] <- temp$p
    }
    
    fisher_df <- fisher_df %>%
      mutate(log10.pval = -log.pval / log(10)) %>%
      arrange(desc(log10.pval))
    
    # Create bar plot using Plotly
    plot_ly(fisher_df,
            x = ~reorder(cell, -log10.pval),
            y = ~log10.pval,
            type = 'bar',
            marker = list(color = 'rgb(158,202,225)', line = list(color = 'rgb(8,48,107)', width = 1.5))) %>%
      layout(
        title = paste("TF: ", selected_tf_name()),
        yaxis = list(title = "-log10(p-value)"),  # Change y-axis label
        xaxis = list(title = "", tickangle = -45)  # Adjust x-axis label angle
      )
  })
  
  # Display data frame filtered based on selected tissue, including only rows where p_value is less than or equal to the specified significance level
  output$estimation_table2 <- renderDT({
    # Ensure tissue is selected
    req(input$select_tissue)
    
    data <- filtered_sctf_df()
    data <- group_by(data, ident)
    data <- summarise(data, mean_score = mean(scaledscore))
    colnames(data)[1] <- "cell_type"
    
    dt <- datatable(data, 
                    caption = htmltools::tags$caption(style='caption-side: top; 
                                                      text-align: center;
                                                      color: black;
                                                      font-size: 100% ;',
                                                      paste("TF: ", selected_tf_name())),
                    options = list(pageLength = 5),
                    rownames = FALSE,
                    selection = list(mode = "single")) %>%
      formatRound(columns = "mean_score", digits = 3)
    dt
  })
  

}

shinyApp(ui, server)
