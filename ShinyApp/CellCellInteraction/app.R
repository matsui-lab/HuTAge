library(shiny)
library(ggplot2)
library(bslib)
library(gridlayout)
library(plotly)
library(DT)
library(flashClust)
library(shinyWidgets)
library(tidyr)
library(dplyr)


cellchat_df_tis <- readRDS("cellchat_cellcell_heat_tis.rds")
cellchat_df_tis <- lapply(cellchat_df_tis,function(x)as.data.frame(x))
cellchat_df_tis <- lapply(cellchat_df_tis,function(x)rename(x, cc=source_target))
lr_prob_tis <- readRDS("cellchat_lr_table_tis.rds")


ui <- page_navbar(
  title = "Cell Cell Interaction",
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
    title = "Cell-Cell",
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
        "0.4fr",
        "1.6fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area1",
        card_header(strong("Tissue of Interest")),
        card_body(
          selectInput(
            inputId = "input_tis",
            label = "Select Tissue",
            choices = names(cellchat_df_tis)
          )
        ) # switchInput(inputId = "isFilterCCI", value = FALSE)
      ),
      grid_card(
        area = "area2",
        card_body(
          grid_container(
            layout = c(
              "area1",
              "area2"
            ),
            row_sizes = c(
              "1.1fr",
              "0.9fr"
            ),
            gap_size = "10px",
            grid_card(
              area = "area1",
              card_header(strong("Strength of Cell-Cell Interactions")),
              card_body(
                plotlyOutput(outputId = "heatmap")
              )
            ),
            grid_card(
              area = "area2",
              card_header(strong("Interaction Probability of LR Pairs Mediating Cell-Cell Interaction")),
              card_body(
                layout_column_wrap(
                  width = 1/2,
                  card(
                    card_body(
                      DTOutput(outputId = "table", width = "100%")
                    )
                  ),
                  card(
                    card_body(
                      plotlyOutput(outputId = "lineplot")
                    )
                  )
                )
              )
            )
          )
        )
      ),
      grid_card(
        area = "area3",
        card_header(strong("Cell of Interest")),
        card_body(
          radioButtons(
            inputId = "Perspective",
            label = "Data Perspective for Heatmap:",
            choices = list("All pairs"="all","Sender"="sender","Receiver"="receiver")
          ),
          selectInput(
            inputId = "Sender",
            label = "Select Sender Cell",
            choices = NULL
          ),
          selectInput(
            inputId = "Receiver",
            label = "Select Receiver Cell",
            choices = NULL
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

server <- function(input, output, session) {
  
  
  observeEvent(input$input_tis, {
    # Update sender and receiver choices when a tissue is selected
    cell_cell_name <- cellchat_df_tis[[input$input_tis]]$cc
    cell_name <- unique(sapply(strsplit(cell_cell_name,"_"),function(x)x[1]))
    updateSelectInput(session, "Sender", choices = cell_name)
    updateSelectInput(session, "Receiver", choices = cell_name, selected = cell_name[2])
  })
  
  
  # Heatmap --------------------------------------------------------------------
  
  
  output$heatmap <- renderPlotly({
    
    cellchat <- cellchat_df_tis[[input$input_tis]]
    rownames(cellchat) <- cellchat$cc
    cellchat$cc <- NULL
    cellchat[is.na(cellchat)] <- 0
    
    if (input$Perspective == "all") {
      clust <- flashClust(dist(cellchat), method = "ward")
      cellchat <- cellchat[clust$order, -c(1:2)]
    }
    if (input$Perspective == "sender") {
      cellchat <- cellchat[cellchat$sender %in% input$Sender, ]
      clust <- flashClust(dist(cellchat), method = "ward")
      cellchat <- cellchat[clust$order, -c(1:2)]
    }
    if (input$Perspective == "receiver") {
      cellchat <- cellchat[cellchat$receiver %in% input$Receiver, ]
      clust <- flashClust(dist(cellchat), method = "ward")
      cellchat <- cellchat[clust$order, -c(1:2)]
    }
    fig <- plot_ly(
      x = colnames(cellchat), # Use condition names for the x-axis
      y = rownames(cellchat), # Use source column for the y-axis
      z = as.matrix(cellchat), # Pass data as matrix to z
      type = 'heatmap',
      zmin = 0, # Minimum value for color display
      zmax = 2 # Maximum value for color display
    )
    fig <- fig %>%
      layout(
        xaxis = list(showticklabels = TRUE),
        yaxis = list(showticklabels = FALSE, tickvals = NULL))
    return(fig)
  })
  
  # Table ----------------------------------------------------------------------
  
  lr_prob_filtered <- reactive({
    req(input$input_tis,input$Sender,input$Receiver)
    input_tissue <- input$input_tis
    input_sender <- input$Sender
    input_receiver <- input$Receiver
    lr_prob <- lr_prob_tis[[input_tissue]]
    # Filtering
    lr_prob <- filter(lr_prob, sender == input_sender)
    lr_prob <- filter(lr_prob, receiver == input_receiver)
    lr_prob <- select(lr_prob,cellcell,lr,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`70-79`)
    return(lr_prob)
  })
  
  
  output$table <- renderDT({
    req(lr_prob_filtered())
    lr_prob_dt <- lr_prob_filtered()
    lr_prob_dt <- select(lr_prob_dt,-cellcell)
    lr_prob_dt <- tibble::column_to_rownames(lr_prob_dt,var = "lr")
    dt <- datatable(lr_prob_dt, 
                    caption = htmltools::tags$caption(style='caption-side: top; 
                                                      text-align: center;
                                                      color: black;
                                                      font-size: 100% ;',
                                                      paste(input$Sender," - ",input$Receiver)),
                    options = list(pageLength=5), 
                    rownames=TRUE,
                    selection=list(mode='single')) %>%
      formatRound(columns=colnames(lr_prob_dt), digits = 3)
    # formatRound(columns = "P-value", digits = 3)  # Corrected here
    dt
  })
  
  # Line Plot ------------------------------------------------------------------
  
  output$lineplot <- renderPlotly({
    req(lr_prob_filtered())
    
    input_tissue <- input$input_tis
    input_sender <- input$Sender
    input_receiver <- input$Receiver
    input_row <- input$table_rows_selected
    
    lr_prob_df <- lr_prob_filtered()
    lr_prob_df <- select(lr_prob_df,-cellcell)
    
    lr_prob_df[is.na(lr_prob_df)] <- 0
    
    # Convert data to long format
    long_df <- lr_prob_df %>%
      pivot_longer(cols = -lr, names_to = "Age", values_to = "Probability")
    
    if (!is.null(input_row)){
      input_lr <- lr_prob_df[input_row, "lr"]
      long_df <- filter(long_df, lr == input_lr)
      title <- paste(input_sender," - ",input_receiver,": ",input_lr)
    } else {
      title <- paste(input_sender," - ",input_receiver,": LR pairs")
    }
    # Create plot with plotly
    p <- plot_ly(data = long_df, x = ~Age, y = ~Probability, type = 'scatter', mode = 'lines+markers',
                 color = ~lr, text = ~lr, hoverinfo = 'text+x+y') %>%
      layout(title = title,
             xaxis = list(title = "Age"),
             yaxis = list(title = "Probability"),
             hovermode = "closest")
    p
  })
  
  
}

shinyApp(ui, server)



