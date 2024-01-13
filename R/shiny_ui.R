#' CellTrans2 Shiny App UI
#'
#' UI definition for the CellTrans2 Shiny App, including input widgets for configuring the analysis.
#' @export

#call: shinyApp(shiny_ui, shiny_server)
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
      numericInput("tauSlider", "Value of tau", value = 1),
      selectInput("dropdown", "Select a cellstate for a more detailed investigation:", choices=1),
      checkboxInput("showTimeToEquilibrium", label = "Show time to equilibrium", value = FALSE),
      conditionalPanel(
        condition = "input.showTimeToEquilibrium",
        textInput("initialDistributionForEquilibrium", "Enter initial distribution (comma-separated, sum = 1)"),
        textInput("tolerance", "The tolerance deviation from the equilibrium distribution")
      ),
      checkboxInput("showPlot", label = "Create celltrans plot", value = FALSE),
      checkboxInput("showPlot_PDF", label = "Create celltrans plot as PDF", value = FALSE),
      conditionalPanel(
        condition = "input.showPlot_PDF",
        textInput("pathToPlot", "Enter a path, where the files will be stored")
      ),
      actionButton("loadDataBtn", "Calculate")
    ),
    column(
      width = 9,
      plotOutput("networkPlot", width = "100%", height = "1200px"),
      plotOutput("plot", width = "50%", height = "500px")
    )
  )
)
