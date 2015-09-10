## P R U E B A S ##
library("RMySQL")
library("genefilter")
library("oligo")
source("OncoProstApp/functions_library.R")
source("OncoProstApp/connection.R")
#source("OncoProstApp/list.gnames.R")

## VARIABLES
con  <- getConnection()
table <- 'taylor_miRNA'
group <- 'gleason_grade_T'
gene <- 'hsa-let-7a'



### Getting tables
queryBP <- paste0("SELECT ", group, ", expr_value FROM ", table, "_pheno, ", 
                  table, "_expr WHERE geo_accession=gsm_id AND probe_id IN (SELECT probe_id FROM ", 
                  table, "_feature WHERE gene_symbol LIKE '", gene, "');")
expr_boxplot <- as.data.frame(dbGetQuery(con, queryBP))


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


