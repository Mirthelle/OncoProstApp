# PAQUETES Y SCRIPTS NECESARIOS
library(shiny)
library(RMySQL)
source("functions_library.R")
source("connection.R")
source("list.gnames.R")

#######################################
##     ONCOPROSTAPP SHINY CODE       ##
#######################################
## Creating getting database connection (function in conection.R)
con  <- getConnection()

shinyServer(function(input, output) {

  ## Gene list selection for SINGLE GENE ANALYSIS
  output$gnames_list <- renderUI ({
    if (is.null(input$database))
      return()
    
    switch(input$database,
           "none" = selectInput("gnames",
                                "Select gene name:",
                                choices = c()),
           "grasso_GPL6480_feature" = selectInput("gnames", 
                                                  "Select gene name:",
                                                  choices = list.gnames("grasso_GPL6480_feature")),
           "taylor_GPL10264_feature" = selectInput("gnames",
                                                   "Select gene name:",
                                                   choices = list.gnames("taylor_GPL10264_feature")),
           "taylor_GPL8227_feature" = selectInput("gnames",
                                                  "Select miRNA name:",
                                                  choices = list.gnames("taylor_GPL8227_feature")),
           "tomlins_GPL2013_feature" = selectInput("gnames",
                                                   "Select gene name:",
                                                   choices = list.gnames("tomlins_GPL2013_feature"))
           )
          
  })
  
  ## DRAWING BOXPLOT
  df_bp <- reactive ({

    # Getting database name
    table <- unlist(strsplit(input$database, "[_]"))
    table <- paste(table[1], "_", table[2], sep='')
    
    # Query for obtaining groups and expression values
    queryBP <- paste0("SELECT ", input$group_by, ", expr_value FROM ", table, "_pheno, ", 
                      table, "_expr WHERE geo_accession=gsm_id AND spot_id IN (SELECT probe_id FROM ", 
                      table, "_feature WHERE gene_symbol LIKE '", input$gnames, "');")
    expr_boxplot <- dbGetQuery(con, queryBP)
  
    df_boxplot <- df_boxplot.build_df_boxplot(expr_boxplot, group_by=input$group_by)
    return(df_boxplot)
  })
  
  output$boxplot <- renderPlot ({
    
#     table <- unlist(strsplit(input$database, "[_]"))
#     table <- paste(table[1], "_", table[2], sep='')
#     
#     # Query for obtaining groups and expression values
#     queryBP <- paste0("SELECT ", input$group_by, ", expr_value FROM ", table, "_pheno, ", 
#                       table, "_expr WHERE geo_accession=gsm_id AND spot_id IN (SELECT probe_id FROM ", 
#                       table, "_feature WHERE gene_symbol LIKE '", input$gnames, "') ORDER BY ", 
#                       input$group_by, " ASC")
#     expr_boxplot <- dbGetQuery(con, queryBP)
#     
#     df_boxplot <- df_boxplot.build_df_boxplot(expr_boxplot, group_by=input$group_by)
    
    boxplot(df_bp(), col="violetred4", main=paste(input$gnames, "by", input$group_by), 
            ylab="Expression values", las=2)
  })
  
#   output$text <- renderPrint({
#     input$gnames
#   })
#   
  output$summary <- renderPrint({
    summary(df_bp())
  })
 
})

