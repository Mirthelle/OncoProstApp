#######################################
##        PACKAGES AND SCRIPTS       ##
#######################################
library("shiny")
library("RMySQL")
#library("genefilter")
library("oligo")
library("affy")
library("affyPLM")
library("limma")
#library("MDA")
library("survival")
library("gplots")
#library("hwriter")
library("mclust")
library("IDPmisc")
library("lme4")
library("coin")
library("amap")
source("functions_library.R")
source("connection.R")

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
    queryBP <- paste0("SELECT disease_status, gleason_grade_T, TNM, expr_value FROM ", input$database, "_pheno, ", 
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
    par(mar = c(20, 4, 4, 2) + 0.1)
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
    fit <- lm(expr_value ~ disease_status + gleason_grade_T + TNM, data=all_expr_df())
    anova(fit)

  })
  
  output$text <- renderPrint({
    class(df_bp())
    #attach(df_bp())
    table(colnames(df_bp()))
  })
  

  ###########################################################################################################
  ## 4. SURVIVAL ANALYSIS
  output$survival <- renderPlot ({
    p_query <- sprintf(paste0("SELECT * FROM ", input$database, "_pheno"))
    p <- dbGetQuery(con, p_query)
    
    par(mfrow=c(3,1))
    
    N1 <- length(unique(p$disease_status))
    plot(survfit(Surv(survival_time, event) ~ disease_status, data = p), 
         main = "Survival according to disease status", lty = 1:N1, 
         col = 1:N1, ylab = "Probability", xlab = "Survival Time in Months")
    legend("bottomleft", legend=unique(p$disease_status), lty=1:N1, col=1:N1, horiz=FALSE, bty='n')

    N2 <- length(unique(p$gleason_grade_T))
    plot(survfit(Surv(survival_time, event) ~ gleason_grade_T, data = p), 
         main = "Survival according to Gleason Grade", lty = 1:N2, 
         col = 1:N2, ylab = "Probability", xlab = "Survival Time in Months")
    legend("bottomleft", legend=unique(p$gleason_grade_T), lty=1:N2, col=1:N2, horiz=FALSE, bty='n')
    
    N3 <- length(unique(p$TNM))
    plot(survfit(Surv(survival_time, event) ~ TNM, data = p), 
         main = "Survival according to TNM stage", lty = 1:N3, 
         col = 1:N3, ylab = "Probability", xlab = "Survival Time in Months")
    legend("bottomleft", legend=unique(p$TNM), lty=1:N3, col=1:N3, horiz=FALSE, bty='n')
    
    })

 
})

