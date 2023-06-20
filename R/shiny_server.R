# Server logic
server <- function(input, output, session) {
  
  #set the path, where the transition matrices are stored to the working directory
  #setwd("~/OneDrive/Dokumente/Master_Angewandte_Informatik/2. Semester/FuE2/CellTrans-1/case_studies/SW620/transition_matrices")
  setwd("~/OneDrive/Dokumente/Master_Angewandte_Informatik/2. Semester/FuE2/CellTrans-1/case_studies/ADAPT_CD133CD44/transition_matrices")
  # Dynamic UI for cell type names
  output$cellTypes <- renderUI({
    cellnr <- input$cellnr
    lapply(seq_len(cellnr), function(i) {
      textInput(paste0("cellType", i), paste("Name of cell type", i))
    })
  })
  
  global <- reactiveValues(datapath = getwd())
  
  MC <- reactiveVal()
  
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
    cellnr <- input$cellnr
    timenr <- input$timenr
    cell_types <- cellTypes
    tau <- input$tauSlider
    
    isTrMatrix <- function(A) {
      sumisOne=TRUE
      for (i in 1:nrow(A) ) { if (abs(sum(A[i,])-1)>1E-8) {sumisOne=FALSE}
      }
      return((all(A>=0)) & (sumisOne))
    }
    
    QOM<-function(A) {
      n=nrow(A)
      for (i in 1:n) 
      {
        repeat 
        { 
          
          non_negative_components=sum(A[i,]>0)  
          diff=(sum(A[i,])-1)/non_negative_components
          for (j in 1:n) {
            if (A[i,j]>0) {
              A[i,j]=A[i,j]-diff
            }
          }
          if (all(A[i,]>=0)) {break}
          for (j in 1:n) {
            if (A[i,j]<0) A[i,j]=0 }
        }
      }
      return(A)
    }
    
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
    
    
    calculate_transitionMatrix1 <- function(M,t,used_timepoints, tau) {
      n=ncol(M) #Rang of Matrix
      transitionMatrix=0
      countQOM=0
      
      for (i in 1:length(t)) {
        #first submatrix in M contains initial matsrix, calculate inverse
        invInitialMatrix=solve(M[1:n,]) #
        
        #calculate steps
        # if(i == 1) {
        #   real_timeDifference <- t[i]
        # }else {
        #   real_timeDifference <- t[i] - t[i-1]
        # }
        
        real_time <- t[i]
        k <- real_time/tau
        
        #derive transition matrices beginning with second submatrix in M
        if (t[i] %in% used_timepoints ) {
          
          if (t[i]>1) {
            Ptemp <- expm( (1/(k)) * logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]) ,method="Eigen")) # berechnet für t[i]==8 die 8. Matrixwurzel, für t[i]==6 die 6. M.W. usw.
          } else {    Ptemp=(invInitialMatrix%*%M[(i*n+1):((i+1)*n), ])%^%(1/t[i])}
          
          if (isTrMatrix(Ptemp)==FALSE) {
            assign( paste("P",t[i],sep=""),QOM(Ptemp)  )
            countQOM=countQOM+1
          } else {    assign( paste("P",t[i],sep=""),Ptemp)}
          
          transitionMatrix=transitionMatrix+get( paste("P",t[i],sep=""))
        }
      }
      transitionMatrix <- (transitionMatrix/(length(used_timepoints)))/tau
      
      
      for (i in 1:n) {
        transitionMatrix[i, i] <- 1 - sum(transitionMatrix[i, -i])
      }
      return(transitionMatrix)
    }
    
    
    #function celltransitions:
      datapoints <- timepoints
      # Print the data points
      # Rest of your code using the datapoints variable
      #timepoints <- dlgList(title = "Data point(s) for estimation", multiple = TRUE, choices = input$timepoints)$res
      trMatrix <- calculate_transitionMatrix1(expData, timepoints, datapoints, tau)
      MC <- new("markovchain", states = cell_types, transitionMatrix = trMatrix, name = "Markov Chain")

      
      print("Results of CellTrans")
      print("################################")
      paste("used timepoints: ", datapoints)
      paste("tau: ", tau)
      print(MC)
      print("predicted equilibrium distribution")
      print(steadyStates(MC))
      print("##########################################")
      
      
      # Convert markovchain object to matrix format
      M <- as.matrix(MC@transitionMatrix)
      nodes <- data.frame(id = 1:cellnr)
      links <- data.frame(from = integer(0), to = integer(0))
    
      # Nested loop to generate edges 
      for (i in 1:nrow(M)) {
        for (j in 1:ncol(M)) {
          if (M[i, j] > 0) {
            edge <- data.frame(from = i, to = j, weight = M[i, j])
            links <- rbind(links, edge)
          }
        }
      }
      
      
      output$networkPlot <- renderPlot({
        net <- graph.data.frame(links, nodes, directed=T)
        #net <- simplify(net, remove.loops = F)
        V(net)$label <- cell_types
        E(net)$label <- E(net)$weight
        E(net)$label.size <- .1
        E(net)$label <- round(E(net)$label, digits = 4)
        E(net)$label.color <- "black"
        E(net)$arrow.size <- .5
        E(net)$edge.color <- "gray80"
        #E(net)$width <- E(net)$weight/6 + 0.25
        #plot(net)
        E(net)$width <- E(net)$weight * 10 / 3 + 0.5
        #line_colors <- sample(colors(), cellnr)
        V(net)$color <- colrs <- c("tomato", "gold", "green", "lightblue")
        edge.start <- get.edges(net, 1:ecount(net))[,1]
        edge.col <- V(net)$color[edge.start]
        # Perform Sugiyama layout to get fixed positions
        layout <- layout_with_sugiyama(net)
        
        # Plot the network with fixed node positions
        plot(net, layout = layout$layout[, 1:2], edge.arrow.size = 0.3, edge.curved = 0.1, edge.color = edge.col)
        })
  })
    session$onSessionEnded(function() {
    stopApp()
  })
}