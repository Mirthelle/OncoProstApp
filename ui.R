library(shiny)
source("server.R")

#########################################
## Generar la interfaz de usuario (UI) ##
#########################################

shinyUI(fluidPage(
  
  ## App title
  titlePanel("OncoProstApp"),
  
  navbarPage("MENU",
             
             ##############################################
             tabPanel("Gene - miRNA correlation analyisis",
             ##############################################
             sidebarLayout(
               sidebarPanel(
                 selectInput("miRNA",
                             "Select your miRNA of interest:",
                             choices = list.gnames("taylor_miRNA")
                             ),
                 selectInput("gene",
                             "Select your gene of interest:",
                             choices = sign.list.gnames("taylor_genes")
                             )
               ),
               mainPanel(
                     plotOutput("scatterplot", width=500, height=500),
                     tags$h4("T-test for the correlation between miRNA and gene:"),
                     verbatimTextOutput("ttest")

               )
             )
             ),
             
             ###########################################
             tabPanel("Diferential expression analysis",
             ###########################################
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("database2", 
                                       "Prostate Cancer Database",
                                       c("Taylor (miRNAs)" = "taylor_miRNA",
                                         "Taylor (genes)" = "taylor_genes"
                                       )
                          )
                        ),
                        mainPanel(
                          tabsetPanel(
                            #                             tabPanel("table",
                            #                                      tags$h4("Genes differentially expressed: "),
                            #                                      tableOutput("difftable")
                            #                                      ),
                            tabPanel("QQplots",
                                     plotOutput("qqplots", width=800, height=900)
                            ),
                            tabPanel("Histograms",
                                     plotOutput("histpvalues", width=800, height=900)
                            )
                          )
                          
                        )
                      )
             ),
             
             #############################
             tabPanel("Survival analysis",
             #############################
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("database4", 
                                       "Prostate Cancer Database",
                                       c("Taylor (miRNAs)"="taylor_miRNA",
                                         "Taylor (genes)"="taylor_genes"
                                       )
                          )
                        ),
                        mainPanel(
                          tags$h4("Kapplan-Meyer Survival Plots:"),
                          plotOutput("survival", width="800px", height = "2000px")
                        )
                      )
                      )
             
             )
  )
)