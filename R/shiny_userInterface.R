library(shiny)
library(shinyFiles)
library(DiagrammeR)
library(igraph)
library(devtools)
library(readr)
library(tcltk)
library(tcltk2)
library(visNetwork)
library(plotly)
library(expm)

ui <- fluidPage(
  titlePanel("CellTrans Shiny App"),
  fluidRow(
    column(
      width = 3,
      # Input widgets
      numericInput("cellnr", "Number of cell states", value = 4),
      uiOutput("cellTypes"),
      numericInput("timenr", "Number of time points", value = 5),
      fileInput("inputFile", "Select a text file, where the time points are stored"),
      textInput("inputFolder", "Enter the path (folder), where the transition matrices are stored"),
      checkboxInput("identityMatrix",
                    label = "Identity matrix (pure initial cell compositions)", value = TRUE),
      numericInput("tauSlider", "Value of tau", value = 1),
      # Hier wird das Dropdown-Menü hinzugefügt
      selectInput("dropdown", "Select a cellstate for a more detailed investigation:", choices=1),
      # Hier wird die Schaltfläche "Calculate" hinzugefügt
      actionButton("loadDataBtn", "Calculate")
    ),
    column(
      width = 9,
      plotOutput("networkPlot", width = "100%", height = "1200px"),
      plotOutput("plot", width = "50%", height = "500px")
    )
  )
)
