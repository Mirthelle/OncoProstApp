library(shiny)
source("server.R")

#########################################
## Generar la interfaz de usuario (UI) ##
#########################################

shinyUI(fluidPage(
  
  ## Apps title
  titlePanel("OncoProstApp"),
  
  navbarPage("MENU",
             #######################################
             tabPanel("Single Gene/miRNA analysis",
             #######################################
                       sidebarLayout(
                         sidebarPanel(
                           
                           # Database selection
                           radioButtons("database", 
                                        "Select a Prostate Cancer Database:",
                                        c("Taylor (miRNAs)"="taylor_miRNA",
                                          "Taylor (genes)"="taylor_genes"
                                          ),
                                        selected="taylor_miRNA"
                                        ),
                           # Gene/miRNA selection
                           wellPanel(uiOutput("gnames_list")),
                           
                           # Group by selection
                           wellPanel(radioButtons("group_by",
                                                  "Show results by:",
                                                  choices = c("Disease Status" = "disease_status",
                                                              "Gleason Grade" = "gleason_grade_T",
                                                              "TNM Stage" = "TNM"
                                                              )
                                                  )
                                     )
                           ),
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Boxplot", 
                                      tags$h4("Boxplot for expression values:"),
                                      plotOutput("boxplot", height = "800px")
                                      ),
                             tabPanel("Data Summary",
                                      tags$h4("Summary of the data:"),
                                      verbatimTextOutput("summary")
                                      ),
                             tabPanel("Linear Model Summary",
                                      tags$h4("Summary of the linear model:"),
                                      verbatimTextOutput("lm")
                                      ),
                             tabPanel("ANOVA test",
                                      tags$h4("ANOVA test: "),
                                      verbatimTextOutput("anova.test")
                                      ),
                             tabPanel("tests",
                                      tags$p("Dynamic input value:"),
                                      verbatimTextOutput("text")
                                      )
                           )
                         )
                       )
                      ),
             
             ###########################################
             tabPanel("Diferential expression analysis",
             ###########################################
                      sidebarLayout(
                        sidebarPanel(
                           radioButtons("database", 
                                        "Prostate Cancer Database",
                                        c("Taylor (miRNAs)"="taylor_miRNA",
                                          "Taylor (genes)"="taylor_genes"
                                          )
                          )
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("table",
                                     tags$h4("Genes differentially expressed: "),
                                     tableOutput("difftable")
                                     ),
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
             
             ##############################################
             tabPanel("Gene - miRNA correlation analyisis",
             ##############################################
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("database", 
                                       "Prostate Cancer Database",
                                       c("Grasso"="grasso",
                                         "Taylor"="taylor",
                                         "Tomlins"= "tomlins"
                                       )
                          )
                        ),
                        mainPanel(
                        )
                      )
                      ),
             
             #############################
             tabPanel("Survival analysis",
             #############################
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("database", 
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