# Server logic
server <- function(input, output, session) {
  # Dynamic UI for cell type names
  output$cellTypes <- renderUI({
    cellnr <- input$cellnr
    lapply(seq_len(cellnr), function(i) {
      textInput(paste0("cellType", i), paste("Name of cell type", i))
    })
  })
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = '~'),
    filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
  )
  
  global <- reactiveValues(datapath = getwd())
  
  output$dir <- renderText({
    paste("test:", global$datapath)
  })
  
  output$dir1 <- renderText({
    paste(input$dir)
  })
  
  ## change smth here... if output$dir is null, display getwd() but it doesn't work
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
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
    }
    
    if(input$identityMatrix){
      expData <- matrix(0, nrow = cellnr * (timenr + 1), ncol = cellnr)
      expData[1:cellnr, ] <- diag(cellnr)
    }
    
    cellnr <- input$cellnr
    timenr <- input$timenr
    
    # Output the result
    output$output <- renderPrint({
      paste("cellnr:", cellnr,
            "; cellTypes:", paste(cellTypes, collapse = ", "),
            "; timeunits:", timeunits,
            "; timenr:", timenr,
            "; timepoints:", paste(timepoints, collapse = ", "))
    })
    
    input_path <- global$datapath
    
    naturalOrder <- function(x) {
      as.numeric(gsub("[^[:digit:]]", "", x))
    }
    
    j <- 0
    if (file.exists(input_path) && !file.info(input_path)$isdir) {
      # If input_path is a file, treat it as a single timepoint
      expData[(j*cellnr+1):((j+1)*cellnr),] <- matrix(scan(input_path, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
      while (!isTrMatrix(expData[(j*cellnr+1):((j+1)*cellnr),])) {
        dlgMessage(paste("Try again! Selected file does not contain an initial setup matrix of dimension ", cellnr, "!"))
        expData[(j*cellnr+1):((j+1)*cellnr),] <- matrix(scan(dlgOpen(title = paste0("Select cell distribution matrix at t=",timepoints[j+1],"."))$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
      }
      j <- j + 1
      } else if (dir.exists(input_path)) {
        #input_path <- dirname(input_path)
        # If input_path is a directory, read files and sort them
        files <- list.files(input_path, full.names = TRUE)
        files <- files[order(naturalOrder(files))]
        for (i in seq_along(files)) {
          t <- timepoints[i]
          j <- j + 1
          expData[(j*cellnr+1):((j+1)*cellnr),] <- matrix(scan(files[i], n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
          while (!isTrMatrix(expData[(j*cellnr+1):((j+1)*cellnr),])) {
            dlgMessage(paste("Try again! Selected file does not contain an initial setup matrix of dimension ", cellnr, "!"))
            expData[(j*cellnr+1):((j+1)*cellnr),] <- matrix(scan(dlgOpen(title = paste0("Select cell distribution matrix at t=",t,"."))$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
          }
        }
      } else {
        dlgMessage("Invalid input path! Please select a file or a directory.")
        return(NULL)
      }

    
    isTrMatrix <- function(A) {
      sumisOne=TRUE
      for (i in 1:nrow(A) ) { if (abs(sum(A[i,])-1)>1E-8) {sumisOne=FALSE}
      }
      return((all(A>=0)) & (sumisOne))
    }  
  
    })
  
    session$onSessionEnded(function() {
    stopApp()
  })
}