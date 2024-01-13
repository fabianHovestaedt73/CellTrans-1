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
    timeunits <- "days"
    timenr <- input$timenr
    inputFile <- input$inputFile
    input_path <- input$inputFolder
    cellnr <- input$cellnr
    timenr <- input$timenr
    cell_types <- cellTypes
    tau <- input$tauSlider
    
    showEquilibrium <- input$showTimeToEquilibrium
    distributionForEquilibrium <- input$initialDistributionForEquilibrium
    tolerance <- input$tolerance
    
    showPlot <- input$showPlot
    showPlot_PDF <- input$showPlot_PDF
    pathToPDF <- input$pathToPlot


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
    # Printing the results of CellTrans.
    # This lines derives and prints the transition probabilities and the predicted equilibrium distribution of the cell state proportions.
    datapoints <- timepoints
    trMatrix <- calculate_transitionMatrix(expData, timepoints, datapoints, tau)
    print("trMatrix:")
    print(trMatrix)
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
    
    
    
    #___________________________________________________________________________
    datapoints <- timepoints
    if(showPlot){
      celltrans_plot(expData, timepoints, datapoints, tau, cellnr, cell_types, timeunits)
    }
    
    if(showEquilibrium){
      distributionForEquilibrium <- as.numeric(strsplit(distributionForEquilibrium, ", ")[[1]])
      timeToEquilibrium(expData, timepoints, datapoints, cell_types, timeunits, distributionForEquilibrium, tolerance)
    }
    
    if(showPlot_PDF){
      celltrans_plotPDF(expData, timepoints, datapoints, cell_types, timeunits, cellnr, pathToPDF, tau)
    }
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
      # trying to create distribution matrices from calculated rate matrix


      #---> todo: Anstatt aus P die Verteilungsmatrizen W zu erzeugen, werden nun aus Q die W's erzeugt. Dafür: Matrixexponential

      # #calc new distribution matrices
      # final_tau <- 0.001
      # trMatrix <- calculate_transitionMatrix(expData, timepoints, datapoints, final_tau)
      # Ratio_M <- (1 / final_tau) * (trMatrix - I)
      
      # Ratenmatrix Q
      # Q <- matrix(c(-0.034637232, 0.0153688445, 1.569748e-02, 0.003570910,
      #               0.040104199, -0.0468738520, 7.006824e-05, 0.006699584,
      #               0.083745811, 0.0006095134, -9.475599e-02, 0.010400662,
      #               0.007921515, 0.1293745739, 4.135663e-02, -0.178652717), nrow = 4, byrow = TRUE)
      # 
      # # Anfangsverteilung: Einheitsmatrix vom Rang 4
      # pi_0 <- diag(4)
      # 
      # # Iterate over time steps from 0 to 100
      # for (t in 0:80) {
      #   # Matrixexponentialfunktion berechnen
      #   expQt <- expm::expm(Q * t)
      #   
      #   # # Verteilung zum Zeitpunkt t berechnen
      #   # pi_t <- pi_0 %*% expQt
      #   
      #   # Ausgabe
      #   cat("Zeitschritt:", t, "\n")
      #   cat("Matrixexponentialfunktion:\n")
      #   print(expQt)
      #   # cat("Verteilung zum Zeitpunkt t:\n")
      #   # print(pi_t)
      #   cat("=============================================\n")
      # }
      
      #_________________________________________________________________________
      
      
      
      # #_________________________________________________________________________
      # #trying to recreate distribution matrices from P
      # 
      # # Transition probability matrix P
      # P <- matrix(c(0.96789393, 0.014941325, 0.0139548884, 0.003209861,
      #               0.03796042, 0.956398647, 0.0001849221, 0.005456010,
      #               0.07402486, 0.001521336, 0.9150223562, 0.009431451,
      #               0.01270999, 0.099836926, 0.0323102065, 0.855142879), nrow = 4, byrow = TRUE)
      # 
      # # Initial distribution
      # pi_0 <- diag(4)
      # 
      # # Time steps
      # t <- 5
      # 
      # # Calculate distribution after t time steps
      # pi_t <- pi_0 %*% (P^t)
      # 
      # # Output
      # # print("Initial Distribution:")
      # # print(pi_0)
      # print(paste("Distribution after", t, "time steps:"))
      # print(pi_t)

      browser()

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
