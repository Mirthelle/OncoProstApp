##                                                                 ##
##   F U N C T I O N S  L I B R A R Y   O N C O P R O S T A P P    ##
##                                                                 ##
#####################################################################
#source("connection.R")

###############################################
##              CBIND FUNCTION               ##
###############################################
cbind.fill <- function(...){
###############################################
# Used in df_nirvana_boxplot.R                #
# in order to fill the data.frame.            #
# https://gist.github.com/abelsonlive/4112423 #
###############################################
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#####################################
##      BOXPLOT DATA.FRAME         ##
#####################################
df_boxplot.build_df_boxplot <- function(result, group_by)  {
#####################################
# This frunction creates the needed #
# data.frame for drawing a boxplot  #
#####################################

  df <- result
  # Important to define df$group_by as df[,group_by]
  colss <- unique(df[,group_by])
  #rm(df2)
  df2 <- data.frame(NA)     # empty data.frame to fill with data
  
for (i in 1:length(colss)) {
    v <- df[which(df[,group_by]==colss[i]), 2]  # Pair gene with data expression in the same row
    
    df2 <- cbind.fill(df2, as.data.frame(v)) # Fill new data.frame
  }
  
  df2 <- df2[, colSums(is.na(df2)) != nrow(df2)]
  colss2 <- colss[!is.na(colss)]
  colnames(df2) <- colss2

  return(df2)
}

#####################################
expr.matrix <- function(table) {
#####################################
# Obtains database expression table #
# and converts it into an expression# 
# matrix                            #
#####################################
  query_expr <- sprintf(paste0("SELECT * FROM ", table, "_expr"))
  
  e <- dbGetQuery(con, query_expr)
  expr <- e$expr_value
  
  new_expr <- matrix(expr, length(unique(e$probe_id)), length(unique(e$gsm_id)), byrow=F)
  rownames(new_expr) <- unique(e$probe_id)
  colnames(new_expr) <- unique(e$gsm_id)
  
  return(new_expr)
}

#####################################
list.gnames <- function(table) {
#####################################
# Selects gene names from a given   #
# database and returns them as a    #
# character vector.                 #
#####################################
  query_gnames <- sprintf(paste0("SELECT distinct gene_symbol FROM ", table, "_feature ORDER BY gene_symbol ASC"))
  list_gnames <- dbGetQuery(con, query_gnames)
  return(list_gnames[,1])
}

sign.list.gnames <- function(table) {
  query_gnames <- paste0("SELECT distinct gene_symbol FROM ", table, "_feature")
  list_gnames <- dbGetQuery(con, query_gnames)
  list_small <- sample(list_gnames[,1], 1000)
  list_small <- list_small[order(list_small)]
  return(list_small)
  #return(list_gnames[,1])
}