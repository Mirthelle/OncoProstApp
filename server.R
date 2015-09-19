#######################################
##        PACKAGES AND SCRIPTS       ##
#######################################
library("shiny")
library("RMySQL")
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
library("shinysky")
source("functions_library.R")
source("connection.R")
load(file="data/taylor_genes_gnames.RData")

#######################################
##     ONCOPROSTAPP SHINY CODE       ##
#######################################
## Creating getting database connection (function in conection.R)
con  <- getConnection()



shinyServer(function(input, output) {

###########################################################################################################
#### CREATING SELECTIONS
###########################################################################################################

  
###########################################################################################################
#### REACTIVE ENVIRONMENTS
###########################################################################################################

eb <- reactive ({
  ## Getting expression and pheno data
  e <- expr.matrix(input$database2)
  p_query <- sprintf(paste0("SELECT * FROM ", input$database2, "_pheno"))
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
  eb <- eBayes(lmc);
  return(eb)
})

data_corr <- reactive ({
  e.gene <- expr.matrix("taylor_genes")
  e.miRNA <- expr.matrix("taylor_miRNA")
  
  query <- "SELECT gene_symbol FROM taylor_genes_feature"
  feat.genes <- dbGetQuery(con, query)
  rownames(e.gene) <- feat.genes[,1]
  
  query2 <- "SELECT sample_id, geo_accession FROM taylor_genes_pheno"
  query3 <- "SELECT sample_id, geo_accession FROM taylor_miRNA_pheno"
  list_genes <- dbGetQuery(con, query2)
  list_mirna <- dbGetQuery(con, query3)
  
  list_genes <- list_genes[order(list_genes$geo_accession),]
  list_mirna <- list_mirna[order(list_mirna$geo_accession),]
  
  e.gene <- e.gene[, order(colnames(e.gene))]
  e.miRNA <- e.miRNA[, order(colnames(e.miRNA))]
  
  colnames(e.gene) <- list_genes[,1]
  colnames(e.miRNA) <-  list_mirna[,1]

  data.gene <- t(t(e.gene[input$gene,]))
  colnames(data.gene) <- 'expr_gene'
  sample_g <- rownames(data.gene)
  data.gene <- data.frame(cbind(sample_g, data.gene))
  
  data.mirna <- t(t(e.miRNA[input$miRNA,]))
  colnames(data.mirna) <- 'expr_mirna'
  sample_m <- rownames(data.mirna)
  data.mirna <- data.frame(cbind(sample_m, data.mirna))
  
  data <- merge(data.mirna, data.gene, by.x = "sample_m", by.y = "sample_g")
  data$expr_mirna <- as.numeric(data$expr_mirna)
  data$expr_gene <- as.numeric(data$expr_gene)
  return(data)
})


###########################################################################################################
#### OUTPUTS
###########################################################################################################ç

###########################################################################################################
## 2. CORRELATION BETWEEN GENES AND MIRNAS
output$scatterplot <- renderPlot ({
  data <- data_corr()
  reg <- lm(data$expr_mirna~data$expr_gene)
  plot(data$expr_mirna, data$expr_gene, type="p", 
       main = paste("Correlation between ", input$gene, " and ", input$miRNA),
       xlab = "miRNA expression values", ylab = "Gene expression value", col="darkblue")
  abline(reg, col='red')
})

output$ttest <- renderText ({
  data <- data_corr()
  t.test(as.vector(data[,2]), as.vector(data[,3]))
})

  
###########################################################################################################
## 2. DIFFERENTIAL EXPRESSION ANALYSIS

# output$difftable <- renderTable ({
#   results <- decideTests(eb(), method="global")
#   summary(results)
# })

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

output$restable <- renderDataTable ({
  e <- expr.matrix(input$database2)
  
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
    p_query <- sprintf(paste0("SELECT * FROM ", input$database4, "_pheno"))
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

