## P R U E B A S ##
library("RMySQL")
library("oligo")
source("OncoProstApp/functions_library.R")
source("OncoProstApp/connection.R")
#source("OncoProstApp/list.gnames.R")
con  <- getConnection()

## VARIABLES

table <- 'taylor_miRNA'
group <- 'gleason_grade_T'
gene <- 'hsa-let-7a'



### Getting tables
queryBP <- paste0("SELECT ", group, ", expr_value FROM ", table, "_pheno, ", 
                  table, "_expr WHERE geo_accession=gsm_id AND probe_id IN (SELECT probe_id FROM ", 
                  table, "_feature WHERE gene_symbol LIKE '", gene, "');")
expr <- as.data.frame(dbGetQuery(con, queryBP))


### DF boxplot
df<-expr_boxplot
# Important to define df$group_by as df[,group_by]
colss <- unique(df[,group])
#rm(df2)
df2 <- data.frame(NA)     # empty data.frame to fill with data
for (i in 1:length(colss)) {
  v <- df[which(df[,group]==colss[i]), 2]  # Pair gene with data expression in the same row
  
  df2 <- cbind.fill(df2, as.data.frame(v))
}

df2 <- df2[, colSums(is.na(df2)) != nrow(df2)]
colss2 <- colss[!is.na(colss)]
colnames(df2) <- colss2


## Correlation between genes and miRNA
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

mirna <- 'hsa-let-7a'
gene <- 'A2M'

data.gene <- t(t(e.gene[gene,]))
colnames(data.gene) <- 'expr_gene'
sample_g <- rownames(data.gene)
data.gene <- as.data.frame(cbind(sample_g, data.gene))

data.mirna <- t(t(e.miRNA[mirna,]))
colnames(data.mirna) <- 'expr_mirna'
sample_m <- rownames(data.mirna)
data.mirna <- as.data.frame(cbind(sample_m, data.mirna))

data <- merge(data.mirna, data.gene, by.x = "sample_m", by.y = "sample_g")
data$expr_mirna <- as.numeric(data$expr_mirna)
data$expr_gene <- as.numeric(data$expr_gene)

reg <- lm(data$expr_mirna~data$expr_gene)
plot(data$expr_mirna, data$expr_gene)
abline(reg, col='red')

cor(data$expr_mirna, data$expr_gene)
summary(t.test(data$expr_mirna, data$expr_gene))

