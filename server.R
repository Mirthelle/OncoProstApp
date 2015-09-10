#######################################
##        PACKAGES AND SCRIPTS       ##
#######################################
library("shiny")
library("RMySQL")
#library("genefilter")
library("oligo")
source("functions_library.R")
source("connection.R")
#source("list.gnames.R")

#######################################
##     ONCOPROSTAPP SHINY CODE       ##
#######################################
## Creating getting database connection (function in conection.R)
con  <- getConnection()



shinyServer(function(input, output) {

###########################################################################################################
#### CREATING SELECTIONS
###########################################################################################################
  ## Gene list selection for SINGLE GENE ANALYSIS ## 
  output$gnames_list <- renderUI ({
    if (is.null(input$database))
      return()
    
    switch(input$database,
           "none" = selectInput("gnames",
                                "Select gene name:",
                                choices = c()),
           "taylor_genes" = selectInput("gnames",
                                        "Select gene name:",
                                        choices = sign.gene.list("taylor_genes")),
           "taylor_miRNA" = selectInput("gnames",
                                        "Select miRNA name:",
                                        choices = list.gnames("taylor_miRNA_feature")),
           )
          
  })
  
###########################################################################################################
#### REACTIVE ENVIRONMENTS
###########################################################################################################

  ###################################################
  expr_df <- reactive ({
  ###################################################
  # Creates a data frame with all data needed for   #
  # generating a boxplot of a single gene/miRNA     #
  # grouped by GG, pathological stage or disease    #
  # type.                                           #
  ###################################################
    # Query for obtaining groups and expression values
    queryBP <- paste0("SELECT ", input$group_by, ", expr_value FROM ", input$database, "_pheno, ", 
                      input$database, "_expr WHERE geo_accession=gsm_id AND probe_id IN (SELECT probe_id FROM ", 
                      input$database, "_feature WHERE gene_symbol LIKE '", input$gnames, "');")
    expr_boxplot <- as.data.frame(dbGetQuery(con, queryBP))
    return(expr_boxplot)
  })

  ###################################################
  all_expr_df <- reactive ({
  ###################################################
  # Obtains a data frame for a single gene/miRNA    #
  # including all groups                            #
  ###################################################
    # Query for obtaining groups and expression values
    queryBP <- paste0("SELECT gleason_grade_T, pathological_stage, expr_value FROM ", input$database, "_pheno, ", 
                      input$database, "_expr WHERE geo_accession=gsm_id AND probe_id IN (SELECT probe_id FROM ", 
                      input$database, "_feature WHERE gene_symbol LIKE '", input$gnames, "');")
    expr_boxplot <- as.data.frame(dbGetQuery(con, queryBP))
    return(expr_boxplot)
  })

  ###################################################
  df_bp <- reactive ({
  ###################################################
  # Using expr_df() reactive environment builds the #
  # final data frame to be used in creating the     #
  # boxplot.                                        #
  ###################################################
      df_boxplot <- df_boxplot.build_df_boxplot(expr_df(), input$group_by)
      df_boxplot <- df_boxplot[, order(colnames(df_boxplot))]
      return(df_boxplot)
  })

###########################################################################################################
#### OUTPUTS
###########################################################################################################

###########################################################################################################
## 1. SINGLE GENE/miRNA ANALYSIS

  ## Drawing boxplot
  output$boxplot <- renderPlot ({
  
    boxplot(df_bp(), col=darkColors(ncol(df_bp())), main=paste(input$gnames, "by", input$group_by), 
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

