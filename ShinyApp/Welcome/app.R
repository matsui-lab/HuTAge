library(shiny)
library(bslib)
library(gridlayout)
library(shinyWidgets)
library(tidyr)
library(dplyr)
library(markdown)

ui <- page_navbar(
  title = "Welcome",
  collapsible = TRUE,
  theme = bslib::bs_theme(preset = "spacelab"),
  sidebar = sidebar(
    title = "HuTAge",
    id = "tabset-default-id",
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
      inputId = "gotop_tf",
      label = "Transcriptional Factor",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/TFact/', '_self')"
    ),
    actionButton(
      inputId = "gocellcell",
      label = "Cell-Cell Interaction",
      onclick = "window.open('https://133.6.53.210:3939/HuTAge/CCI/', '_self')"
    )
  ),
  nav_panel(
    title = "",
    includeMarkdown("welcome.md")  # Ensure this file is in the correct directory
  ),
  nav_spacer(),
  nav_panel(
    title = "User Guide",
    tags$iframe(style = "height:100%; width:100%; border:none;", src = "20240729_Hutage_userguide.pdf")
  )
)

server <- function(input, output, session) {
  
}

shinyApp(ui, server)



