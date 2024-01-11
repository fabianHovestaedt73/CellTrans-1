#' CellTrans2 Shiny App UI
#'
#' UI definition for the CellTrans2 Shiny App, including input widgets for configuring the analysis.
#' @export
#' call: shinyApp(shiny_ui, shiny_server)


shiny_ui <- fluidPage(
  titlePanel("CellTrans Shiny App"),
  fluidRow(
    column(
      width = 3,
      # Input widgets
      numericInput("cellnr", "Number of cell states", value = 4),
      uiOutput("cellTypes"),
      numericInput("timenr", "Number of time points", value = 5),
      fileInput("inputFile", "Select a text file, where the time points are stored"),
      textInput("inputFolder", "Enter the path (folder), where the distribution matrices are stored"),
      checkboxInput("identityMatrix",
                    label = "Identity matrix (pure initial cell compositions)", value = TRUE),
      selectInput("dropdown", "Select a cellstate for a more detailed investigation:", choices=1),
      actionButton("loadDataBtn", "Calculate")
    ),
    column(
      width = 9,
      plotOutput("networkPlot", width = "100%", height = "1200px"),
      plotOutput("plot", width = "50%", height = "500px")
    )
  )
)
