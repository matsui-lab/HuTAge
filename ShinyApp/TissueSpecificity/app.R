library(shiny)
library(ggplot2)
library(bslib)
library(gridlayout)
library(DT)
library(plotly)
library(shinythemes)
library(plotly)
library(flashClust)

tispecObj <- readRDS("tau_age.rds")

tau_age <- tispecObj$tau_age

score_mat <- tispecObj$score

gene_id <- tispecObj$gene_id

gene_symbols <- sort(gene_id$external_gene_name)

tauAnno_by_age <- readRDS("tauAnno_by_age.rds")
for(i in seq_along(tauAnno_by_age)){
  tmp <- tauAnno_by_age[[i]]
  colnames(tmp) <- gsub("\\."," ",colnames(tmp))
  tauAnno_by_age[[i]] <- tmp
}

tauAnno_by_age <- lapply(tauAnno_by_age, function(x) x[x$external_gene_name %in% rownames(score_mat), ])
tissueList_tau <- names(tauAnno_by_age$`20-29`)[-c(1:3)]

ui <- page_navbar(
  title = "Tissue Specificity",
  selected = "By Gene",
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
    title = "By Gene",
    #    selected = "By Gene",,,,,,,,,,,
    grid_container(
      layout = c(
        "area1 area2",
        "area3 area2"
      ),
      row_sizes = c(
        "0.75fr",
        "1.25fr"
      ),
      col_sizes = c(
        "0.5fr",
        "1.5fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        full_screen = TRUE,
        card_header(strong("Age of Interest")),
        card_body(
          selectInput(
            inputId = "age_select",
            label = "Select Age",
            choices = c(
              "20-29",
              "30-39",
              "40-49",
              "50-59",
              "60-69",
              "70-79"
            ),
            selected = "20-29"
          ),
          "",
          strong("Tau Score Threshold Range"),
          sliderInput(
            inputId = "tau_score_range",
            label = "",
            min = 0,
            max = 1,
            value = c(
              0.8,
              1
            ),
            width = "100%",
            step = 0.05,
            ticks = TRUE,
            animate = FALSE
          )
        )
      ),
      grid_card(
        area = "area2",
        card_header(strong("Tau Scores")),
        card_body(
          card(
            card_body(
              plotlyOutput(outputId = "tau_age", height = "450px")
            )
          ),
          card(
            card_body(
              plotlyOutput(outputId = "tau_linePlot", height = "150px")
            )
          )
        )
      ),
      grid_card(
        area = "area3",
        full_screen = TRUE,
        card_header(strong("Tau Score Table")),
        card_body(
          strong("Select a Gene"),
          DTOutput(
            outputId = "tauTable",
            width = "100%",
            height = "auto"
          )
        )
      )
    )
  ),
  nav_panel(
    title = "By Tissue",
    grid_container(
      layout = c(
        "area1 area2"
      ),
      col_sizes = c(
        "0.5fr",
        "1.5fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        full_screen = TRUE,
        card_header(strong("Tissue of Interest")),
        card_body(
          selectInput(
            inputId = "selectTissue",
            label = "Select Tissue",
            choices = sort(tissueList_tau)
          ),
          numericInput(
            inputId = "tau_thres_finder",
            label = "Tissue Fraction Score",
            value = 0.8,
            max = 1,
            min = 0,
            step = 0.01
          )
        )
      ),
      grid_card(
        area = "area2",
        card_header(textOutput(outputId = "tau_age_tissue_title")),
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
            "0.7fr",
            "1.3fr"
          ),
          gap_size = "10px",
          grid_card(
            area = "area1",
            card_body(
              plotlyOutput(
                outputId = "tau_age_tissue",
                height = "600px"
              )
            )
          ),
          grid_card(
            area = "area2",
            card_body(
              DTOutput(outputId = "combinedTable", height = "200px")
            )
          ),
          grid_card(
            area = "area3",
            card_body(
              plotlyOutput(outputId = "tissue_fraction_barplot")
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
  
  # -------------------------
  # By Gene 
  # -------------------------
  
  selectData <- reactive({
    
    withProgress(message = 'Preparing data...', value = 0, {
      
      input_mat <- score_mat
      incProgress(0.2, detail = "Selecting data based on gene...")
      
      if (is.null(input$age_select) || length(input$age_select) == 0) {
        return(input_mat)
      }
      
      if (is.null(input$tau_score_range) || length(input$tau_score_range) == 0) {
        return(input_mat)
      }
      
      # Get the selected range
      tau_range <- input$tau_score_range
      
      select_mat <- input_mat[
        input_mat[,input$age_select] >= tau_range[1] & input_mat[,input$age_select] <= tau_range[2],
        ,drop=F
      ]
      
      if (nrow(select_mat) == 0) {
        return(NULL)
      }
      
      incProgress(0.2, detail = "Performing hierarchical clustering...")
      
      hc <- flashClust(dist(select_mat), method = "ward")
      
      incProgress(0.3, detail = "Finalizing data preparation...")
      
      select_mat[hc$order,]
    })
  })
  
  # Render heatmap
  output$tau_age <- renderPlotly({
    
    input_mat <- selectData()  # Get data for the heatmap
    
    if (is.null(input_mat)) {
      return(NULL)
    }
    
    p <- plot_ly(x = colnames(input_mat), y = rownames(input_mat), z = input_mat, type = "heatmap", source = "selectHeatmap") %>%
      layout(coloraxis = list(colorbar = list(title = "Tau score")),
             yaxis = list(showticklabels = F, tickvals = NULL))  # Hide y-axis labels
    event_register(p, 'plotly_selected')
    
    return(p)
  })
  
  output$tauTable <- renderDT({
    
    select_mat <- selectData()
    
    DT::datatable(select_mat, selection = list(mode = 'single', selected = 1), extensions = c('Buttons'),
                  options = list(
                    pageLength = 5, 
                    scrollX = TRUE, 
                    dom = 'Blfrtip', 
                    buttons = list(
                      list(extend = 'csv', text = 'CSV', exportOptions = list(modifier = list(page = 'all'))),
                      list(extend = 'excel', text = 'Excel', exportOptions = list(modifier = list(page = 'all')))
                    )
                  ))
    
  }, server = FALSE)
  
  output$tau_linePlot <- renderPlotly({
    
    select_mat <- selectData()
    
    # Get the index of the selected row
    selected_row <- input$tauTable_rows_selected
    
    # Ensure a row is selected
    req(length(selected_row) == 1)
    
    # Get the selected gene name
    selected_gene <- rownames(select_mat)[selected_row]
    
    # Extract data corresponding to the selected gene
    gene_data <- select_mat[selected_gene, , drop = FALSE]
    
    # Prepare data frame for plotly
    plot_data <- data.frame(
      Age = colnames(select_mat),
      Score = as.numeric(gene_data)
    )
    
    # Render the plot
    p <- plot_ly(plot_data, x = ~Age, y = ~Score, type = 'scatter', mode = 'lines+markers',
                 line = list(shape = "linear")) %>%
      layout(title = paste("Gene:", selected_gene),
             xaxis = list(title = "Age"),
             yaxis = list(title = "Score"))
    
    return(p)
  })
  
  # -------------------------
  # By Tissue 
  # -------------------------
  
  filterData <- reactive({
    selectAge <- input$age_select
    tau_thres <- input$tau_thres_finder
    selectedTissue <- input$selectTissue  # Selected tissue name
    
    incProgress(0.2, detail = "Selecting data based on age...")
    
    selectTauMat <- tauAnno_by_age[[selectAge]]
    selectTauMat <- selectTauMat[selectTauMat$gene_biotype == "protein_coding", ]
    
    ridx <- !duplicated(selectTauMat$external_gene_name)
    # Remove duplicates
    if (any(!ridx)) {
      selectTauMat <- selectTauMat[ridx, ]
    }
    
    rownames(selectTauMat) <- selectTauMat$external_gene_name
    
    incProgress(0.3, detail = paste("Filtering data...", nrow(selectTauMat), "genes remained.."))
    
    selectTauMat <- selectTauMat[,-c(1,2,3)]
    selectTauMat <- as.matrix(selectTauMat)
    
    # Filter rows where the tau score for the selected tissue is above the threshold
    filteredData <- selectTauMat[selectTauMat[,selectedTissue] >= tau_thres, ]
    
    incProgress(0.4, detail = "Performing hierarchical clustering...")
    
    hc <- flashClust(dist(filteredData), method = "ward")
    
    incProgress(0.5, detail = "Finalizing data preparation...")
    
    filteredData <- filteredData[hc$order, ]
    
    return(filteredData)
  })
  
  preparedData <- reactive({
    withProgress(message = 'Preparing data...', value = 0, {
      
      selectAge <- input$age_select
      tau_thres <- input$tau_score
      
      selectTauMat <- filterData()
      
      incProgress(0.3, detail = paste("Filtering data...", nrow(selectTauMat), "genes remained.."))
      
      selectTauMat <- selectTauMat[,-c(1,2,3)]
      selectTauMat <- as.matrix(selectTauMat)
      
      incProgress(0.2, detail = "Performing hierarchical clustering...")
      
      hc <- flashClust(dist(selectTauMat), method = "ward")
      
      incProgress(0.3, detail = "Finalizing data preparation...")
      
      selectTauMat[hc$order,]
    })
  })
  
  output$tau_age_tissue_title <- renderText({
    req(input$selectTissue)
    paste("Tau Fraction Score:", input$age_select)
  })
  
  output$tau_age_tissue <- renderPlotly({
    selectTauMat <- filterData()  # This now calls the reactive expression
    plot_ly(x = colnames(selectTauMat), y = rownames(selectTauMat), z = selectTauMat, type = "heatmap") %>%
      layout(coloraxis = list(colorbar = list(title = "Tau score")),
             xaxis = list(tickangle = -45), yaxis = list(showticklabels = F, tickvals = NULL))
  })
  
  output$combinedTable <- renderDT({
    selectTauMat <- filterData()
    DT::datatable(selectTauMat, selection = list(mode = 'single', selected = 1),
                  options = list(pageLength = 5, scrollX = TRUE))
  }, server = FALSE)  # client-side processing
  
  selectedGene <- reactiveVal()
  # When a row in combinedTable is clicked, update the selectedGene
  observeEvent(input$combinedTable_rows_selected, {
    req(input$combinedTable_rows_selected)  # Ensure that a row is actually selected
    data <- filterData()  # Fetch the current data displayed in the combinedTable
    gene <- rownames(data)[input$combinedTable_rows_selected]  # Extract the selected gene symbol
    selectedGene(gene)  # Update the reactive value
  })
  
  # Render the tissue fraction bar plot based on the selected gene
  output$tissue_fraction_barplot <- renderPlotly({
    req(selectedGene())  # Require a selected gene to proceed
    data <- filterData()  # Get the current filtered data
    
    # Ensure the selected gene is still in the data
    if (!selectedGene() %in% rownames(data)) {
      return(NULL)
    }
    
    # Prepare the data for the plot
    plot_data <- as.data.frame(t(data[selectedGene(), , drop = FALSE]))
    colnames(plot_data) <- "score"
    
    plot_data$tissue <- rownames(plot_data)
    plot_data$score <- as.numeric(plot_data$score)
    plot_data$color <- ifelse(plot_data$tissue == input$selectTissue, 'rgb(255,127,14)', 'rgb(158,202,225)')
    
    plot_data <- plot_data %>%
      arrange(desc(score))
    
    plot_data$tissue <- factor(plot_data$tissue, levels = plot_data$tissue, ordered = TRUE)
    
    # Create the bar plot
    plot_ly(data = plot_data, x = ~tissue, y = ~score, type = 'bar', name = 'Tau Score',
            marker = list(color = plot_data$color, line = list(color = 'rgb(8,48,107)', width = 1.5))) %>%
      layout(title = paste("Gene: ", selectedGene()),
             xaxis = list(title = "Tissue", tickangle = -45),
             yaxis = list(title = "Score"),
             showlegend = FALSE)
  })

  
}

shinyApp(ui, server)


