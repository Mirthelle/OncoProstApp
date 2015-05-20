library(shiny)
source("server.R")

#########################################
## Generar la interfaz de usuario (UI) ##
#########################################

shinyUI(fluidPage(
  
  ## Título de la App
  titlePanel("OncoNirvana 2.0"),
  
  navbarPage("MENU"),
  
  sidebarLayout(
    sidebarPanel(
      #el objeto uiOutuput lo obtengo del servidor
      uiOutput("choose_inputs")     
      ),
    
    mainPanel(   
      tabsetPanel(             
        tabPanel("BoxPlot",                 
                 plotOutput("TaylorBoxPlot", height = "400px")
                 ),
        tabPanel("Descriptive statistics",
                 verbatimTextOutput("descriptive_statistic"),
                 tags$style(type='text/css', "statistics {background-color: rgba(0,0,255,0.10); color: black;}")
                 ),
        tabPanel("T test",
                 verbatimTextOutput("ttest")
                 ),
        tabPanel("Summary", tableOutput("summary"))
      )
          #h3(textOutput("caption")),
          #plotOutput("TaylorBoxPlot",height = "400px")
      
    )
)
))