# PAQUETES Y SCRIPTS NECESARIOS
library(shiny)
library(RMySQL)
library(genefilter)
source("functions_library.R")
source("connection.R")
source("list.gnames.R")

#######################################
##     ONCOPROSTAPP SHINY CODE       ##
#######################################
## Creating getting database connection (function in conection.R)
con  <- getConnection()

shinyServer(function(input, output) {

  ## Gene list selection for SINGLE GENE ANALYSIS ## 
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
                                                   choices = sign.gene.list("taylor_GPL10264")),
           "taylor_GPL8227_feature" = selectInput("gnames",
                                                  "Select miRNA name:",
                                                  choices = list.gnames("taylor_GPL8227_feature")),
           "tomlins_GPL2013_feature" = selectInput("gnames",
                                                   "Select gene name:",
                                                   choices = list.gnames("tomlins_GPL2013_feature"))
           )
          
  })
  
  ## REACTIVE ENVIRONMENTS ##
  expr_df <- reactive ({
    # Getting database name
    table <- unlist(strsplit(input$database, "[_]"))
    table <- paste(table[1], "_", table[2], sep='')
    
    # Query for obtaining groups and expression values
    queryBP <- paste0("SELECT ", input$group_by, ", expr_value FROM ", table, "_pheno, ", 
                      table, "_expr WHERE geo_accession=gsm_id AND spot_id IN (SELECT probe_id FROM ", 
                      table, "_feature WHERE gene_symbol LIKE '", input$gnames, "');")
    expr_boxplot <- as.data.frame(dbGetQuery(con, queryBP))
    return(expr_boxplot)
  })
  
  all_expr_df <- reactive ({
    # Getting database name
    table <- unlist(strsplit(input$database, "[_]"))
    table <- paste(table[1], "_", table[2], sep='')
    
    # Query for obtaining groups and expression values
    queryBP <- paste0("SELECT gleason_grade_2, pathological_stage, expr_value FROM ", table, "_pheno, ", 
                      table, "_expr WHERE geo_accession=gsm_id AND spot_id IN (SELECT probe_id FROM ", 
                      table, "_feature WHERE gene_symbol LIKE '", input$gnames, "');")
    expr_boxplot <- as.data.frame(dbGetQuery(con, queryBP))
    return(expr_boxplot)
  })
  
  df_bp <- reactive ({
      df_boxplot <- df_boxplot.build_df_boxplot(expr_df(), group_by=input$group_by)
      df_boxplot <- df_boxplot[, order(colnames(df_boxplot))]
      return(df_boxplot)
  })
  
  ## OUTPUTS ##
  ## Drawing boxplot
  output$boxplot <- renderPlot ({
  
    boxplot(df_bp(), col="violetred4", main=paste(input$gnames, "by", input$group_by), 
            ylab="Expression values", las=2)
  })
  
  ## Descriptive summary
  output$summary <- renderPrint({
    summary(df_bp())
  })
  
  ## T test summary
  output$lm <- renderPrint ({
    fit <- lm(expr_value ~ get(input$group_by), data=expr_df())
    summary(fit)
  })
  
  ## Anova summary  
  output$anova.test <- renderPrint({
    fit <- lm(expr_value ~ gleason_grade_2 + pathological_stage, data=all_expr_df())
    anova(fit)

  })
  
#   output$text <- renderPrint({
#     df_bp()
#   })
  

 
})

