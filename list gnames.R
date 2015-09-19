
library("RMySQL")

source("OncoProstApp/functions_library.R")
source("OncoProstApp/connection.R")

con  <- getConnection()
table <- "taylor_genes"
query_gnames <- sprintf(paste0("SELECT gene_symbol, gb_acc FROM ", table, "_feature"))
list_gnames <- dbGetQuery(con, query_gnames)

taylor_genes_gnames <- list_gnames
save(taylor_genes_gnames, file="OncoProstApp/data/taylor_genes_gnames.RData")

query_gnames <- sprintf(paste0("SELECT distinct(gene_symbol) FROM ", table, "_feature"))
list_gnames <- as.vector(dbGetQuery(con, query_gnames))
list_gnames <- as.vector(list_gnames[,1])
list <- sample(list_gnames, 300)
