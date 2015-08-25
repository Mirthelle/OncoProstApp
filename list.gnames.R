source("connection.R")

list.gnames <- function(table) {
  query_gnames <- sprintf(paste("SELECT distinct gene_symbol FROM ", table, " ORDER BY gene_symbol ASC"))
  list_gnames <- dbGetQuery(con, query_gnames)
  return(list_gnames[,1])
}

# source("OncoProstApp/connection.R")
# table <- "taylor_GPL10264"
# table <- "taylor_GPL8227"
# group <- "gleason_grade_2"

expr.matrix <- function(table) {
  query_expr <- sprintf(paste0("SELECT * FROM ", table, "_expr"))
  
  e <- dbGetQuery(con, query_expr)
  expr <- e$expr_value
  
  new_expr <- matrix(expr, length(unique(e$spot_id)), length(unique(e$gsm_id)), byrow=F)
  rownames(new_expr) <- unique(e$spot_id)
  colnames(new_expr) <- unique(e$gsm_id)
  
  return(new_expr)
}

sign.gene.list <- function(table, group) {
  e <- expr.matrix(table)
  query_pheno <- sprintf(paste0("SELECT geo_accession, biopsy_gleason_grade FROM ", table, "_pheno"))
  p <- dbGetQuery(con, query_pheno)
  
  query_feature <- sprintf(paste0("SELECT probe_id, gene_symbol FROM ", table, "_feature"))
  f <- dbGetQuery(con, query_feature)
  
  pvals <- rowFtests(e, factor(p$biopsy_gleason_grade))
  head(pvals)
  
  fdr <- p.adjust(pvals$p.value, method="fdr")
  names(fdr)  <- rownames(e)
  sig <- fdr[fdr<0.0001]
  genes <- names(sig)
  index <- which(f$probe_id %in% genes)
  gnames <- f[index, 2]
  gnames <- gnames[!is.na(gnames)]
  return(gnames)
}

