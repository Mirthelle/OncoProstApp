######################################################################
##                                                                  ##
##   F U N C T I O N S  L I B R A R Y   4   O N C O N I R V A N A   ##
##                                                                  ##
######################################################################

###############################################
##              FUNCIÓN CBIND                ##
###############################################
cbind.fill <- function(...){
###############################################
# Usada por df_nirvana_boxplot.R              #
# para completar el data.frame.               #
# https://gist.github.com/abelsonlive/4112423 #
###############################################
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

A <- c("a", "b", "c", "d")
B <- c("e", "f", "g", "h")

prueba <- cbind.fill(A, B)

#####################################
##    DATA.FRAME PARA BOXPLOT      ##
#####################################
df_nirvana_boxplot.build_df_boxplot <- function(result,group_by)  {
#####################################
# Esta función crea el data.frame   #
# necesario para dibujar un boxplot #
# con Shiny en server.R.            #
#####################################

  df<-result
  #es muy importante definir df$group_by con df[,group_by]
  #porque de la primera forma no evalúa el parámetro group_by
  colss<-unique(df[,group_by])  #como si hiciera un distinct
  colss
  rm(df2)
  df2<-data.frame(NA)     #creo un data frame vacío (NA) para volcar los datos.
  for(i in 1:length(colss)){
    v <- df[which(df[,group_by]==colss[i]),c(2)]  #emparejo para cada gen todos sus datos de expresión
    #en una misma fila  
    
    df2 <- cbind.fill(df2, as.data.frame(v,colnames=colss[i])) #relleno el data.frame nuevo
  }
  
  df2 <- df2[,-c(1)]    #quito la primera columna de NA
  
  colnames(df2) <- colss  #pongo los nombres de los genes como nombres de las columnas
  return(df2)
}
