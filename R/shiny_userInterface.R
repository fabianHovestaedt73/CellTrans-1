# Define the UI
ui <- fluidPage(
  # Input widgets
  numericInput("cellnr", "Number of cell states", value = 0),
  uiOutput("cellTypes"),
  selectInput("timeunits", "Time step length", choices = c("minutes", "hours", "days", "weeks", "months", "cell divisions")),
  numericInput("timenr", "Number of time points", value = 0),
  fileInput("inputFile", "Select input file"),
  actionButton("loadDataBtn", "Load Data"),
  
  # Output
  verbatimTextOutput("output")
)