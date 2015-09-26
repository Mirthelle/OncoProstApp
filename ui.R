library(shiny)
source("server.R")

#########################################
## Generar la interfaz de usuario (UI) ##
#########################################

shinyUI(fluidPage(
  
  ## App title
  titlePanel("OncoProstApp"),
  
  navbarPage("MENU",
             
             tabPanel("Home",
                      navlistPanel(
                        tabPanel("What is OncoProstApp",
                                 h3("Welcome to OncoProstApp!"),
                                 p("OncoProstApp is a R-based web application which performs several predifined statistical 
                                  analyses over a selected gene expression prostate cancer database."),
                                 p("The database is formed by four different gene expression databases in prostate cancer. 
                                   The gene expression values were obtained using micro-arrays technology and the data is 
                                   publicly available at GEO portal (for the GEO accession codes, visit the 'About OncoProstApp' 
                                   section."),
                                 p("This web was developed as a master thesis for the MSc in Bioinformatics, at Universidad de Murcia.
                                   ")
                                 ),
                        
                        tabPanel("How to use OncoProstApp",
                                 h3("How to use OncoProstApp"),
                                 p("In the navigation bar you can find three types of analysis that can be performed:"),
                                 tags$ul(
                                   tags$li('Correlation between a gene and a miRNA. This test performs a correlation analysis between a 
                                           user-selected gene and miRNA. Returns a scatterplot and the result of a t-test.'),
                                   tags$li('Differential expression analysis. Returns various plots and data tables which show which genes 
                                           are differentially expressed in the database of interest, adjusting by gleason grade, lymph node
                                           invasion and metastasis.'),
                                   tags$li('Survival analysis. This analysis uses phenotipe data and survival time to draw different 
                                           Kapplay-Meyer plots for type of tumour, gleason grade and TNM stage.')
                                 )
                                 ),
                        
                        tabPanel("About OncoProstApp",
                                 h3("About OncoProstApp"),
                                 h4('Bibliography:'),
                                 tags$ul(
                                   tags$li('Barry S Taylor, Nikolaus Schultz, Haley Hieronymus, Anuradha Gopalan, Yonghong Xiao, Brett S Carver, 
                                            Vivek K Arora, Poorvi Kaushik, Ethan Cerami, Yevgeniy Antipin, Nicholas Mitsiades, Thomas Landers, Igor 
                                            Dolgalev, E John, Manda Wilson, Nicholas D Socci, Alex E Lash, Adriana Heguy, A James, Howard I Scher, 
                                            Victor E Reuter, Peter T Scardino, Chris Sander, L Sawyers, William L Gerald, Mskcc Prostate, and Cancer 
                                            Oncogenome. Integrative genomic proﬁling of human prostate cancer. Cancer cell, 18(1):11–22, 2011. 
                                            doi: 10.1016/j.ccr.2010.05.026.Integrative.'),
                                   tags$li('Catherine S Grasso, Yi-mi Wu, Dan R Robinson, Xuhong Cao, M Saravana, Amjad P Khan, Michael J Quist, 
                                            Xiaojun Jing, J Robert, J Chad Brenner, Irfan a Asangani, Bushra Ateeq, Sang Y Chun, Lee Sam, Matt 
                                            nstett, Rohit Mehra, John R Prensner, Gregory a Ryslik, Fabio Vandin, Benjamin J Raphael, P Lakshmi, 
                                            Daniel R Rhodes, Kenneth J Pienta, and Arul M Chinnaiyan. The Mutational Landscape of Lethal Castrate 
                                            Resistant Prostate Cancer. Nature, 487(7406):239–243, 2013. doi: 10.1038/ nature11125.The.'),
                                   tags$li('Scott a Tomlins, Rohit Mehra, Daniel R Rhodes, Xuhong Cao, Lei Wang, Saravana M Dhanasekaran, Shanker 
                                            Kalyana-Sundaram, John T Wei, Mark a Rubin, Kenneth J Pienta, Rajal B Shah, and Arul M Chinnaiyan. 
                                            Integrative molecular concept modeling of prostate cancer progression. Nature genetics, 39(1):41–51, 
                                            2007. ISSN 1061-4036. doi: 10.1038/ng1935.')
                                 ),
                                 h4('Source Code'),
                                 p('You can download and consult the source code for this Shiny web app in the following link: '),
                                 tags$a("https://github.com/Mirthelle/OncoProstApp")
                                )
                        )
                      ),
             
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
                                         "Taylor (genes)" = "taylor_genes", 
                                         "Grasso" = "grasso", 
                                         "Tomlins" = "tomlins"
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
                            ),
                            tabPanel("Differentially Expressed Genes",
                                      downloadButton('downDiffTable', 'Download'),
                                      tableOutput("restable")
                                     ),
                            tabPanel("Heat maps"
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
                                         "Taylor (genes)"="taylor_genes",
                                         "Grasso" = "grasso",
                                         "Tomlins" = "tomlins"
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