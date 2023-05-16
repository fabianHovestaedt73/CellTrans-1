#' Reading all necessary data.
#' 
#' This function opens a dialog box which asks for the number of cell types, the names of the cell types, the time step length and the time points of measurement. Then, the files containing the cell state proportion matrices can be selected. First, the initial experimental setup matrix can either be chosen as identity matrix (for pure initial cell populations) or a custom initial matrix can be provided. Then, the experimental cell proportion matrices are read for each time point of measurement. It is recommended to save the input into a variable for the further analysis, e.g. input <- readExperimentalData().
#' @keywords initial experimental matrix, cell distribution matrices
#' @export

library(devtools)
readExperimentalData <- function()  {
  
  dlgMessage("Welcome to CellTrans!\n Please assure that you have prepared appropriate files containing the cell state distribution matrices representing your experimental data.")
  cellnr  <- 	as.integer(dlgInput("Number of cell states")$res)
  #Read cell type names
  cell_types=rep.int(0,cellnr)
  for (i in 1:cellnr)  {
    cell_types[i] <- dlgInput(paste("Name of cell type", i))$res
  }
  
  
  
  
  
  #Ask for timeunits and timepoints
  timeunits<-dlgList(title="Time step length",c("minutes","hours","days","weeks","months","cell divisions"))$res
  timenr<-as.integer(dlgInput("Number of time points")$res)
  
  # timepoints=rep(0,timenr) #AUTOMATISIEREN
  # for (i in 1:timenr)  {
  #   timepoints[i] <-dlgInput(paste0("Timepoint ",i))$res
  # }
  # timepoints=as.numeric(timepoints)
  
  #Ask for period and starting timepoint
  period <- as.integer(dlgInput("Period of time points")$res)
  start_tp <- as.integer(dlgInput("Starting time point")$res)
  
  #Create timepoints vector based on period and starting timepoint
  timepoints <- seq(start_tp, (timenr - 1) * period + start_tp, by = period)
  
  #Print the timepoints vector  
  timepoints <- as.numeric(timepoints)
  
  print(timepoints)
  
  #Create matrix for cell distribution matrices incuding initial experimental matrix
  expData=matrix(0, nrow=cellnr*(timenr+1),ncol=cellnr)
  #Ask for initial input matrix
  res=""
  while (res=="") {
    res <- dlgList(title="Initial experimental setup matrix", choices=c("Identity matrix (pure initial cell compositions)", "Individual matrix"))$res
    if (res=="Identity matrix (pure initial cell compositions)") {
      expData[1:cellnr,]=diag(cellnr)
    } else {
      repeat {
      expData[1:cellnr,]=matrix(scan(dlgOpen(title = "Select initial experimental matrix")$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
      while (!isTrMatrix(expData[1:cellnr,])) 
        {
          dlgMessage(paste("Try again! Selected file does not contain an experimental cell state proportion matrix of dimension ",cellnr," !"))
          expData[1:cellnr,]=matrix(scan(dlgOpen(title = " Select initial experimental matrix.")$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
          while (det(expData[1:cellnr,])==0) 
          {
            dlgMessage(paste("The experimental setup matrix is not valid (not invertible)!"))
            expData[1:cellnr,]=matrix(scan(dlgOpen(title = " Select initial experimental matrix.")$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
          }  
        }
      break}
      
    }
  }
  
  
  #Ask for cell distribution matrices
  # j=0
  # for (t in timepoints) {
  #   j=j+1
  #   expData[(j*cellnr+1):((j+1)*cellnr),]=matrix(scan(dlgOpen(title = paste0("Select cell distribution matrix at t=",t,"."))$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
  #   while (!isTrMatrix(  expData[(j*cellnr+1):((j+1)*cellnr),]  ) ) {dlgMessage(paste("Try again! Selected file does not contain an initial setup matrix of dimension ",cellnr,"!"))
  #     expData[(j*cellnr+1):((j+1)*cellnr),]=matrix(scan(dlgOpen(title = paste0("Select cell distribution matrix at t=",t,"."))$res, n = cellnr*cellnr), cellnr, cellnr, byrow = TRUE)
  #   }
  # }
  # return(list("cellnr"=cellnr, "cell_types"=cell_types, "timeunits"=timeunits, "timenr"=timenr, "timepoints"=timepoints, "experimentalData"=expData))
  # 
  
  # Öffne einen Dialog zum Auswählen einer Datei oder eines Ordners
  # dlg <- dlgFile(title = "Select file or folder", filter = c("All files (*.*)|*.*"), multiple = FALSE)
  # 
  # # Extrahiere den Pfad aus dem Rückgabewert des Dialogs
  # if (dlg$was_cancelled) {
  #   stop("No file or folder selected.")
  # } else if (dlg$is_folder) {
  #   input_path <- dlg$path
  # } else {
  #   input_path <- dirname(dlg$path)
  # }
  
  naturalOrder <- function(x) {
    as.numeric(gsub("[^[:digit:]]", "", x))
  }
  
  library(tcltk)
  
  select_file_or_dir <- function() {
    # Abfrage, ob Datei oder Ordner ausgewählt werden soll
    choice <- utils::menu(c("Datei", "Ordner"), title = "Auswahl treffen")
    
    if (choice == "Datei") {
      # Datei auswählen
      input_path <- choose.file()
    } else {
      # Ordner auswählen
      input_path <- choose.dir()
      print("directory")
      print(input_path)
    }
    
    # Rückgabe des ausgewählten Pfads
    return(input_path)
  }
  
  #select_file_or_dir()
  
  # Prompt user to select a file or folder
  input_path <- "C:/Users/Fabian/OneDrive/Dokumente/Master_Angewandte_Informatik/2. Semester/CellTrans/Supplementary Material_712241/SW620"
  # Prompt user to select a file or folder
  
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
  
  return(list("cellnr"=cellnr, "cell_types"=cell_types, "timeunits"=timeunits, "timenr"=timenr, "timepoints"=timepoints, "experimentalData"=expData))
}
