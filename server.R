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
library("hwriter")
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

eb <- reactive ({
  ## Getting expression and pheno data
  e <- expr.matrix(input$database)
  p_query <- sprintf(paste0("SELECT * FROM ", input$database, "_pheno"))
  p <- dbGetQuery(con, p_query)
  
  # Comparisons using model specification
  
  p_2 <- na.omit(p) # Delete NA rows in pheno
  
  gg2 <- subset(p_2, p_2$gleason_grade_S <= 6, select=colnames(p_2))
  gg3 <- subset(p_2, p_2$gleason_grade_S > 6, select=colnames(p_2))
  
  
  gg_2 <- rep('<= 6', nrow(gg2))
  gg_3 <- rep('> 6', nrow(gg3))
  gg2 <- cbind(gg2, gg_2)
  gg3 <- cbind(gg3, gg_3)
  colnames(gg2)[19] <- "gg_group"
  colnames(gg3)[19] <- "gg_group"
  
  p_3 <- rbind(gg2, gg3)
  
  design <- model.matrix(~ gg_group + LNI + metastasis, data=p_3);
  
  # Delete NA rows for expressions
  e_3 <- e[,p_3$geo_accession]
  
  # Contrasts
  lmf <- lmFit(e_3, design);
  betas <- coef(lmf);
  
  cont.mat <- cbind(
    c(0, 1, 0, 0, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 1)
  )
  
  colnames(cont.mat) <- c("gleason_grade", "LNI", "Metastasis")
  
  ### Moderated t-test
  lmc <- contrasts.fit(lmf, cont.mat);
  eb.0 <- eBayes(lmc);
  return(eb.0)
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
## 2. DIFFERENTIAL EXPRESSION ANALYSIS

output$difftable <- renderTable ({
  results <- decideTests(eb(), method="global")
  summary(results)
})

output$qqplots <- renderPlot ({
  qqpval<-function(p, pch=16, col=4, ...)
  {
    p<-p[!is.na(p)]
    n<-length(p)
    pexp<-(1:n)/(n+1)
    plot(-log(pexp,10), -log(sort(p),10), xlab="-log(expected P value)",
         ylab="-log(observed P value)", pch=pch, col=col, ...)
    abline(0,1,col=2)
  }
  
  eb <- eb()
  
  par(mfrow=c(2, 2));
  qqpval(eb$p.value[, 1], main='Gleason Grade');
  qqpval(eb$p.value[, 2], main='LNI');
  qqpval(eb$p.value[, 3], main='Metastasis');
})

output$histpvalues  <- renderPlot ({
  eb <- eb()
  
  par(mfrow=c(2, 2));
  hist(eb$p.value[, 1], main='Gleason Grade');
  hist(eb$p.value[, 2], main='LNI');
  hist(eb$p.value[, 3], main='Metastasis');
})

output$restable <- renderDatatable ({
  e <- expr.matrix(input$database)
  
  # Results table
  
  ### Anotation
  d_query <- sprintf(paste0("SELECT * FROM ", table, "_feature"))
  d <- dbGetQuery(con, d_query)
  rownames(d)  <- d$probe_id
  
  ### Some global descriptives
  mns <- apply(e, 1, mean);
  sds <- apply(e, 1, sd);
  
  sm <- mclapply(as.list(as.data.frame(t(e_2))), FUN=function(o)
  {
    r <- do.call(rbind, tapply(o, p_3$gg_group, function(o) c(mean(o), sd(o))));
    rn <- rownames(r);
    r <- as.vector(t(r));
    names(r) <- paste(rep(rn, each=2), c(".mean", ".sd"), sep='');
    r;
  }, mc.cores=7);
  sm <- do.call(rbind, sm);
  
  d <- cbind(d, global.mean=mns[rownames(d)], global.sd=sds[rownames(d)],
             sm[rownames(d), ], sm[rownames(d), ]);
  
  ### Tests results
  
  for (j in 1:ncol(cont.mat))
  {
    dr <- topTable(eb.0, adjust="BH", number=nrow(e_2), coef=colnames(cont.mat)[j], confint=T);
    rownames(dr) <- dr$probeset;
    dr <- dr[, c("logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")];
    colnames(dr) <- tolower(colnames(dr));
    colnames(dr)[colnames(dr) == 'logfc'] <- "beta";
    colnames(dr)[colnames(dr) == 'ci.025'] <- "beta.l95";
    colnames(dr)[colnames(dr) == 'ci.975'] <- "beta.u95";
    colnames(dr) <- paste(colnames(dr), colnames(cont.mat)[j], sep='.');
    d <- cbind(d, dr[rownames(d), ]);
  }
  
  colnames(d) <- gsub("p.value", "pval", colnames(d));
  colnames(d) <- gsub("p.val", "pval", colnames(d));
  
  ### Save Rdata and csv
  d <- d[order(d$adj.pval.gleason_grade), ];
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

