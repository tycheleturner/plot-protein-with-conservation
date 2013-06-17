#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Protein Plotting Script with Conservation Shiny UI Script
#Programming Language: R
#Updated 06/14/2013

#Description: This script is the ui.R script required for the Shiny application of Plot Protein With Conservation: Visualization of Mutations

library(shiny)

# Define UI
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Plot Protein With Conservation: Visualization of Mutations"),

	#Sidebar panel
  sidebarPanel(
  	 
	fileInput(inputId="mutationFile", label="Mutation File (required):"),
				 	
	fileInput(inputId="proteinArchitectureFile", label="Protein Architecture File (required):"),

	fileInput(inputId="postTranslationalModificationFile", label="Post Translational Modification File (required):"),
				 
	fileInput(inputId="alignmentFile", label="Alignment File (required):"),
	
	textInput("referenceSequencePositionInFile", "Reference Sequence Position In File (required):", "Data"),
	
	textInput("nameOfQuery", "Name of Mutation File Data:", "Query"),
	
	textInput("tickSize", "Tick Size:", 10),	
	
	checkboxInput(inputId = "labels",
      label = strong("Show Labels?"),
      value = FALSE),
      
    checkboxInput(inputId = "showConservationScore",
      label = strong("Show Conservation Score?"),
      value = FALSE),
  
    checkboxInput(inputId = "showReferenceSequence",
      label = strong("Show Reference Sequence?"),
      value = FALSE),
  
    checkboxInput(inputId = "showGridsAtTicks",
      label = strong("Show Gridlines at Ticks?"),
      value = FALSE),
  
    checkboxInput(inputId = "zoom",
      label = strong("Zoom In?"),
      value = FALSE),

	textInput("zoomStart", "Zoom Start:", "Starting Position"),	

	textInput("zoomEnd", "Zoom End:", "Ending Position"),
	
	checkboxInput(inputId = "second",
      label = strong("Second Mutation File For Plot?"),
      value = FALSE),

	fileInput(inputId="secondmutationFile", label="Second Mutation File (optional):"),			 
	
	textInput("nameOfQuery2", "Name of Second Mutation File Data:", "Second Query")


  ),

  # Show a tabset with plot
  mainPanel(
	h5("Tychele N. Turner, tycheleturner@gmail.com"),
	h5("Laboratory of Aravinda Chakravarti, Ph.D."),
	
    tabsetPanel(
      tabPanel("Plot", plotOutput("plot")), 
      tabPanel("Table", tableOutput("table"))
    )
  )
))
