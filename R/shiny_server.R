#setwd("C:/Users/Fabian/OneDrive/Dokumente/Master_Angewandte_Informatik/3. Semester/FuE/R")
#library("CellTrans")
#shinyApp(ui, server)


# Server logic
server <- function(input, output, session) {
  
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
    inputFolder <- input$inputFolder
    cellnr <- input$cellnr
    timenr <- input$timenr
    cell_types <- cellTypes
    tau <- input$tauSlider
    isTrMatrix <- function(A) {
      sumisOne=TRUE
      for (i in 1:nrow(A) ) { if (abs(sum(A[i,])-1)>1E-5) {sumisOne=FALSE}
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
    
    
    input_path <- inputFolder
    
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
        
        real_time <- t[i]
        k <- real_time/tau
        #derive transition matrices beginning with second submatrix in M
        if (t[i] %in% used_timepoints ) {
          
          if (t[i]>1) {
            #vorher: 
            Ptemp <- expm( (1/(k)) * logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]) ,method="Eigen"))
            #jetzt:
            #Ptemp <- expm( (1/(t[i])) * logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]) ,method="Eigen")) # berechnet für t[i]==8 die 8. Matrixwurzel, für t[i]==6 die 6. M.W. usw.
          } else {    Ptemp=(invInitialMatrix%*%M[(i*n+1):((i+1)*n), ])%^%(1/t[i])}
          
          if (isTrMatrix(Ptemp)==FALSE) {
            assign( paste("P",t[i],sep=""),QOM(Ptemp)  )
            countQOM=countQOM+1
          } else {    assign( paste("P",t[i],sep=""),Ptemp)}
          
          transitionMatrix=transitionMatrix+get( paste("P",t[i],sep=""))
        }
      }
      
      # P/tau = Q + I #fehlerhaft?
      #vorher: 
      #transitionMatrix <- (transitionMatrix/(length(used_timepoints)))/tau
      #jetzt:
      transitionMatrix <- transitionMatrix/(length(used_timepoints))
      
      #jetzt:
      # transitionMatrix <- (transitionMatrix/(length(used_timepoints)))
      
      
      # P/tau - I = Q ... --> Zeilensummen aber 1 und nicht 0
      #vorher:
      # for (i in 1:n) {
      #   transitionMatrix[i, i] <- 1 - sum(transitionMatrix[i, -i])
      # }
      
      # Ergebnis ist keine Ratenmatrix sondern eine Übergangsmatrix mit sehr kleinem Zeittakt?
      return(transitionMatrix)
    }
    #function celltransitions:
    datapoints <- timepoints
    # Print the data points
    # Rest of your code using the datapoints variable
    #timepoints <- dlgList(title = "Data point(s) for estimation", multiple = TRUE, choices = input$timepoints)$res
    browser()
    trMatrix <- calculate_transitionMatrix1(expData, timepoints, datapoints, tau)
    MC <- new("markovchain", states = cell_types, transitionMatrix = trMatrix, name = "Markov Chain")
    
    cat("\n")
    cat("\n")
    print("#################________Results of CellTrans_____########################")
    cat("\n")
    print(paste("_____tau_____ :", tau))
    cat("\n")
    print(MC)
    print("predicted equilibrium distribution")
    print(steadyStates(MC))
    print("##########################################################################")
    cat("\n")
    cat("\n")
    
    # I <- diag(n)
    # Ratio_M <- (1 / tau) * (transitionMatrix - I)
    
    
    # Convert markovchain object to matrix format
    M <- as.matrix(MC@transitionMatrix)
    nodes <- data.frame(id = 1:cellnr)
    links <- data.frame(from = integer(0), to = integer(0))
    
    
    #calc_P =: tau*Q
    
    
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
      n <- ncol(M)
      set.seed(123)
      net <- graph.data.frame(links, nodes, directed = T)
      V(net)$label <- cell_types
      E(net)$label <- E(net)$weight
      E(net)$arrow.size = 2.5
      E(net)$label.cex <- 2
      V(net)$label.cex <- 3
      V(net)$label.color <- "black"
      E(net)$label <- round(E(net)$label, digits = 3)
      E(net)$label.color <- "black"
      E(net)$edge.color <- "gray80"
      V(net)$color <- colrs <- c("tomato", "gold", "green", "lightblue")
      edge.start <- get.edges(net, 1:ecount(net))[, 1]
      edge.col <- V(net)$color[edge.start]
      layout <- layout.circle(net)
      E(net)$width <- E(net)$weight * 100 
      
      for (i in 1:n) {
        E(net)$width[(i - 1) * (n + 1) + 1] <- E(net)$width[(i - 1) * (n + 1) + 1] / 7.5
      }
      # self-loops
      net_no_self_loops <- net
      self_loop_list <- c()
      for (i in 1:n) {
        self_loop_ids <- get.edge.ids(net_no_self_loops, c(i, i))
        self_loop_list <- c(self_loop_list, self_loop_ids)
      }
      # Delete Zero-Edges
      m <- n * n
      delete_edges <- c()
      for (i in 1:m) {
        if (E(net)$weight[i] < 0.005) {
          delete_edges <- c(delete_edges, E(net)[i])
        }
      }
      net_no_self_loops <- net_no_self_loops - E(net_no_self_loops)[c(delete_edges, self_loop_list)]
      diagonal_indices <- seq(1, length(edge.start), by = n + 1)
      delete_edges <- c(delete_edges, diagonal_indices)
      edge.col <- V(net_no_self_loops)$color[edge.start[-delete_edges]]
      plot(net_no_self_loops, layout = layout, edge.arrow.size = 2.5, edge.curved = 0.1, edge.color = edge.col)
    })
    
    
    
    
    
      
    
      matrix_list <- list()
      tau_values <- c(2.5, 2, 1.5, 1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
      
      for (i in 1:length(tau_values)) {
        trMatrix <- calculate_transitionMatrix1(expData, timepoints, datapoints, tau_values[i])
        matrix_list[[i]] <- trMatrix
      }
      
      
      
      #selectInput for plot             
      n <- ncol(M)
      observe({
        updateSelectInput(session, "dropdown", choices = 1:n)
      })
      
      cellstate <- as.integer(input$dropdown)
      vektor <- c()
      for (j in 1:length(matrix_list)){
        vektor <- c(vektor, matrix_list[[j]][cellstate,])
      }
      
      state_values <- list()  # Create a list to store the values for each state
      for (i in 1:n) {
        state_values[[i]] <- vektor[seq(i, length(vektor), 4)]  # Assign the values to each state
      }
      
      # Erstelle eine Liste mit den gewünschten Beschriftungen
      labels <- c()
      for (i in 1:length(state_values)) {
        if (i != cellstate) {
          labels <- c(labels, paste(cellstate, "-->", i))
        }
      }
      
      
      
      
      output$plot <- renderPlot({
        tau_values <- c(2.5, 2, 1.5, 1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
        log_tau_values <- log(tau_values)  # Calculate the logarithm of tau_values
        
        plot(log_tau_values, type = "n", xlim = range(log_tau_values), ylim = range(unlist(state_values[-cellstate])),
             xlab = "log(tau)", ylab="", cex.lab = 2.0, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
        
        colors <- c("tomato", "gold", "green", "blue")  # Define the colors for the lines
        
        for (i in 1:n) {
          if (i != cellstate) {
            y <- state_values[[i]]
            lines(log_tau_values, y, type = "b", pch = 19, col = colors[i])
          }
        }
        
        colors1 <- c()
        for (i in 1:n) {
          if (i != cellstate) {
            colors1 <- c(colors1, colors[i])
          }
        }
        
        set.seed(1) # just to get the same random numbers
        par(xpd=FALSE)
        legend("left",inset=c(0, 0), legend = labels, lty = 1, col = colors1, cex = 1.5)  # Increase the size of the labels
        par(xpd=TRUE)
        par(cex.axis = 3)  # Increase the size of axis labels
        par(cex.lab = 3)  # Increase the size of xlab and ylab
      })
    
  })
  session$onSessionEnded(function() {
    stopApp()
  })
}