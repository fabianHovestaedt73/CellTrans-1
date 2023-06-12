# Server logic
server <- function(input, output, session) {
  # Dynamic UI for cell type names
  output$cellTypes <- renderUI({
    cellnr <- input$cellnr
    lapply(seq_len(cellnr), function(i) {
      textInput(paste0("cellType", i), paste("Name of cell type", i))
    })
  })
  
  # Load data when button is clicked
  observeEvent(input$loadDataBtn, {
    # Retrieve input values
    cellnr <- input$cellnr
    cellTypes <- sapply(seq_len(cellnr), function(i) {
      input[[paste0("cellType", i)]]
    })
    timeunits <- input$timeunits
    timenr <- input$timenr
    inputFile <- input$inputFile
    
    # Process the input file
    if (!is.null(inputFile)) {
      timepoints <- readLines(inputFile$datapath)
      separator <- ifelse(grepl(",", timepoints), ",", ifelse(grepl(";", timepoints), ";", " "))
      timepoints <- as.numeric(strsplit(timepoints, separator)[[1]])
      
      # Perform further data processing or store the data as needed
      
      # Output the result
      output$output <- renderPrint({
        paste("cellnr:", cellnr,
              "; cellTypes:", paste(cellTypes, collapse = ", "),
              "; timeunits:", timeunits,
              "; timenr:", timenr,
              "; timepoints:", paste(timepoints, collapse = ", "))
      })
    }
  })
}
