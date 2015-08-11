######################################################################
##                                                                  ##
##   F U N C T I O N S  L I B R A R Y   4   O N C O N I R V A N A   ##
##                                                                  ##
######################################################################

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
# data.frame for drawing a boxplot #
#####################################

  df<-result
  # Important to define df$group_by as df[,group_by]
  colss <- unique(df[,group_by])
  rm(df2)
  df2 <- data.frame(NA)     # empty data.frame to fill with data
  for (i in 1:length(colss)) {
    v <- df[which(df[,group_by]==colss[i]),c(2)]  # Pair gene with data expression in the same row
    
    df2 <- cbind.fill(df2, as.data.frame(v,colnames=colss[i])) # Fill new data.frame
  }
  
  df2 <- df2[,-1]    # Deleting first column with NA
  
  colnames(df2) <- colss  # Gene names as colnames
  return(df2)
}
