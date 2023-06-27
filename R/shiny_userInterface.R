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

# Define the UI
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
      checkboxInput("identityMatrix", label = "Identity matrix (pure initial cell compositions)", value = TRUE),
      numericInput("tauSlider", "Value of tau", value = 1),
      selectInput("dropdown", "Select a cellstate for a more detailed investigation:", choices = 1),
      actionButton("loadDataBtn", "Calculate")
    ),
    column(
      width = 9,
      plotOutput("networkPlot", width = "100%", height = "900px"), # Increased height to 1000px
      plotOutput("plot", width = "100%", height = "500px"),
      plotOutput("plot1", width = "100%", height = "500px")
      )
  )
)
