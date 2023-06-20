library(shiny)
library(shinyFiles)
library(DiagrammeR)
# Define the UI
ui <- fluidPage(
  titlePanel("CellTrans Shiny App"),
  fluidRow(
    column(
      width = 3,
      # Input widgets
      numericInput("cellnr", "Number of cell states", value = 2),
      uiOutput("cellTypes"),
      numericInput("timenr", "Number of time points", value = 12),
      fileInput("inputFile", "Select a text file, where the time points are stored"),
      checkboxInput("identityMatrix", label = "Identity matrix (pure initial cell compositions)", value = TRUE),
      numericInput("tauSlider", "Value of tau", value = 1),
      actionButton("loadDataBtn", "Calculate")
    ),
    column(
      width = 9, # Adjusted column width to accommodate the full width of the network plot
      plotOutput("networkPlot", width = "100%", height = "800px") # Updated width to "100%"
    )
  )
)