#' CellTrans2 Shiny App Server
#'
#' Server function for the CellTrans2 Shiny App, handling dynamic input fields,
#' data loading, transition matrix calculations, and result visualization.
#'
#' @export
#'
shiny_server <- function(input, output, session) {

  # create's input fields dynamically, depending on the input of cellnr ("Number of cell states")
  output$cellTypes <- renderUI({
    cellnr <- input$cellnr
    lapply(seq_len(cellnr), function(i) {
      textInput(paste0("cellType", i), paste("Name of cell type", i))
    })
  })

  # used for development resons, as you need the current path as working directory
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
    input_path <- input$inputFolder
    cellnr <- input$cellnr
    timenr <- input$timenr
    cell_types <- cellTypes
    tau <- input$tauSlider


    #___________________________________________________________________________
    # Check for transition matrix.
    #
    # Verification whether a matrix is a stochastic matrix of a Markov process or not, i.e. with non-negative entries and row sums equal to one.
    # param A The matrix that is checked to be a transition matrix.
    # keywords: transition matrix, Matrix process

    isTrMatrix <- function(A) {
      sumisOne=TRUE
      for (i in 1:nrow(A) ) { if (abs(sum(A[i,])-1)>1E-5) {sumisOne=FALSE}
      }
      return((all(A>=0)) & (sumisOne))
    }
    #___________________________________________________________________________

    #___________________________________________________________________________
    # Quasi-optimization of a matrix root to a stochastic matrix.
    #
    # QOM finds a transition matrix that is as close as possible, with respect to the Euclidean distance of the rows, to the fractional root of a given transition matrix. If the argument is a transition matrix, this matrix is returned unchanged.
    # parameters: A - The fractional root of a matrix that shall be regularized.
    # keywords: QOM, matrix root regularization

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
    #___________________________________________________________________________


    # read timepoints file
    if (!is.null(inputFile)) {
      timepoints <- readLines(inputFile$datapath)
      separator <- ifelse(grepl(",", timepoints), ",", ifelse(grepl(";", timepoints), ";", " "))
      timepoints <- as.numeric(strsplit(timepoints, separator)[[1]])
    }

    # In most cases, pure cell cultures were used as the starting distribution of the experiments, which is why the experiments in these cases begin with the unit matrix as the distribution matrices. The option to pass a distribution different from pure cell cultures as the default state still needs to be implemented at this point. By default, this value is always set to true.
    if(input$identityMatrix){
      expData <- matrix(0, nrow = cellnr * (timenr + 1), ncol = cellnr)
      expData[1:cellnr, ] <- diag(cellnr)
    }

    # Extract and sort numbers from a string in ascending order.
    naturalOrder <- function(x) {
      as.numeric(gsub("[^[:digit:]]", "", x))
    }


    #___________________________________________________________________________
    # read distribution matrices
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
    #___________________________________________________________________________


    #___________________________________________________________________________
    # Calculation of the transition matrix.
    #
    # This function derives the transition matrix which contains the transition probabilitites between the distinct cell states.
    # param: M - Matrix containing all experimental data: the initial experimental setup matrix and the experimental cell state distribution matrices.
    # param: t - This A vector containing all timepoints at which the cell state distributions have been experimentally measured. ---> realTime
    # param - used_timepoints A vector with the timepoints utilized to estimate the transition matrix. ---> Tau

    # keywords: transition matrix, Matrix process

    calculate_transitionMatrix <- function(M,t,used_timepoints, tau) {
      n=ncol(M) #Rang of Matrix
      transitionMatrix=0
      countQOM=0

      for (i in 1:length(t)) {
        # first submatrix in M contains initial matsrix, calculate inverse
        invInitialMatrix=solve(M[1:n,]) #

        # After the integration of tau, a different matrix root than the original value must be calculated at a certain point in time, so this value is calculated here. For example, for a value tau=0.5, the 4th matrix root is now calculated instead of the 2nd matrix root, as 2/0.5=4.
        real_time <- t[i]
        k <- real_time/tau

        # derive transition matrices beginning with second submatrix in M
        if (t[i] %in% used_timepoints ) {

          if (t[i]>1) {
            Ptemp <- expm( (1/(k)) * logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]) ,method="Eigen"))
          }
          else {    Ptemp=(invInitialMatrix%*%M[(i*n+1):((i+1)*n), ])%^%(1/t[i])}

          if (isTrMatrix(Ptemp)==FALSE) {
            assign( paste("P",t[i],sep=""),QOM(Ptemp)  )
            countQOM=countQOM+1
          }
          else {    assign( paste("P",t[i],sep=""),Ptemp)}

          transitionMatrix=transitionMatrix+get( paste("P",t[i],sep=""))
        }
      }

      transitionMatrix <- transitionMatrix/(length(used_timepoints))

      return(transitionMatrix)
    }
    #___________________________________________________________________________

    #___________________________________________________________________________
    # Printing the results of CellTrans.
    # This lines derives and prints the transition probabilities and the predicted equilibrium distribution of the cell state proportions.
    datapoints <- timepoints
    trMatrix <- calculate_transitionMatrix(expData, timepoints, datapoints, tau)
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
    cat("\n")
    #___________________________________________________________________________



    #calculate Q as ratio matrix
    n <- ncol(trMatrix)
    I <- diag(n)
    Ratio_M <- (1 / tau) * (trMatrix - I)
    for (i in 1:n) {
      Ratio_M[i, i] <- Ratio_M[i, i] * -1
    }


    # Convert markovchain object to matrix format
    M <- as.matrix(Ratio_M)
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

    # calculate and assign edges and vertices for interaction graph from calculated transition matrices
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
        trMatrix <- calculate_transitionMatrix(expData, timepoints, datapoints, tau_values[i])
        Ratio_M <- (1 / tau_values[i]) * (trMatrix - I)
        matrix_list[[i]] <- Ratio_M
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

      # Create a list to store the values for each state
      state_values <- list()
      for (i in 1:n) {
        state_values[[i]] <- vektor[seq(i, length(vektor), 4)]  # Assign the values to each state
      }

      # Create a list with the desired labelling
      labels <- c()
      for (i in 1:length(state_values)) {
        if (i != cellstate) {
          labels <- c(labels, paste(cellstate, "-->", i))
        }
      }

      #rate matrix
      print("Rate matrix:")
      print(Ratio_M)
      print("##########################################################################")
      cat("\n")
      cat("\n")


      #_________________________________________________________________________
      # # trying to create distribution matrices from calculated rate matrix
      # #
      # #calc new distribution matrices
      # final_tau <- 0.001
      # trMatrix <- calculate_transitionMatrix(expData, timepoints, datapoints, final_tau)
      # Ratio_M <- (1 / final_tau) * (trMatrix - I)
      #
      #
      # #---> todo: Anstatt aus P die Verteilungsmatrizen W zu erzeugen, werden nun aus Q die W's erzeugt. Dafür: Matrixexponential
      #
      # print("#####################___tau = 0.001____###############")
      # P <- (Ratio_M * final_tau ) + I
      # print(P)
      #
      #
      # # Create the folder if it does not exist
      # dir.create("new_distribution_matrices", showWarnings = FALSE)
      #
      # # Function for generating a unique file name
      # generateUniqueFilename <- function(matrixName) {
      #   paste("new_distribution_matrices", "/", matrixName, ".txt", sep = "")
      # }
      #
      # W5 <- as.matrix(P^5)
      # W13 <- as.matrix(P^13)
      # W23 <- as.matrix(P^23)
      # W34 <- as.matrix(P^34)
      # W68 <- as.matrix(P^68)
      #
      # # Function for saving a matrix in a text file
      # saveMatrixToFile <- function(matrix, matrixName) {
      #   filename <- generateUniqueFilename(matrixName)
      #   write.table(matrix, file = filename, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      # }
      #
      # # Save the matrices in separate text filesteien
      # saveMatrixToFile(W5, "W5")
      # saveMatrixToFile(W13, "W13")
      # saveMatrixToFile(W23, "W23")
      # saveMatrixToFile(W34, "W34")
      # saveMatrixToFile(W68, "W68")
      #_________________________________________________________________________


      #plot interaction graph
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
