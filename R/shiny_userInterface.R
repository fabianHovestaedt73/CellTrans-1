library(shiny)
library(shinyFiles)
# Define the UI
ui <- fluidPage(
  titlePanel("CellTrans Shiny App"),
  # Input widgets
  numericInput("cellnr", "Number of cell states", value = 0),
  uiOutput("cellTypes"),
  selectInput("timeunits", "Select a timestep unit", choices = c("minutes", "hours", "days", "weeks", "months", "cell divisions")),
  numericInput("timenr", "Number of time points", value = 0),
  fileInput("inputFile", "Select a textfile, where the timepoints are stored"),
  checkboxInput("identityMatrix", label = "Identity matrix (pure initial cell compositions)", value = TRUE),
  shinyDirButton("dir", "Input directory", "Upload"),
  verbatimTextOutput("dir", placeholder = TRUE),
  verbatimTextOutput("dir1", placeholder = TRUE),
  #fileInput("cellDistributionMatrices", label = "Select a file or folder, where the cell state distribution matrices are stored"),
  actionButton("loadDataBtn", "Load Data"),
  
  # Output
  verbatimTextOutput("output")
)
