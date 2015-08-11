source("connection.R")

list.gnames <- function(table) {
  query_gnames <- sprintf(paste("SELECT distinct gene_symbol FROM ", table, " ORDER BY gene_symbol ASC"))
  list_gnames <- dbGetQuery(con, query_gnames)
  return(list_gnames[,1])
}
