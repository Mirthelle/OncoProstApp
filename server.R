# PAQUETES Y SCRIPTS NECESARIOS
library(shiny)
library(RMySQL)
source("functions_library.R")
source("connectionDB.R")

#######################################
## CÓDIGO DE LA APLICACIÓN WEB SHINY ##
#######################################

## Crear o usar conexión abierta de MySQL (función en global.R)
connection <- getConnection()

## Código de las aplicaciones y funciones que se ejecutarán.
shinyServer(function(input, output) {
    #Obtener los tipos de tumores y convertirlos en una lista.
    statement_tumor<-sprintf("SELECT distinct tumor_type FROM taylor_21036")    
    tumor_type <- dbGetQuery(connection, statement_tumor)
    
    # Obteber todos los distintos nombres de miRNAs y convertirlos en una lista.
    statement_id_gpl <- sprintf("SELECT distinct id_gpl FROM taylor_exp_21036 order by id_gpl asc")  
    did_gpl <- dbGetQuery(connection, statement_id_gpl)
    #como devuelve un data.frame, cojo sólo la 1ª columna para construir la select input
    
    # Crear una lista con todos los select input que se quieren mostrar.
    #se la paso a renderUI, ya que sólo puedo llamarlo una vez
    input_list <- list(selectizeInput("gen", "Choose microRNA", choices = did_gpl[,1], multiple = FALSE, options=NULL),
                       radioButtons("result_by", label = h4("Show results by"),
                                    choices = list("Disease Status" = "disease_status", 
                                                   "Gleason Grade" = "biopsy_gleason_grade",
                                                   "Pathological Stage" = "pathological_stage"
                                                   ), 
                                    selected = "disease_status")
                       )
#                        checkboxGroupInput("tumor_type", label = h4("Tipos de tumor a analizar"),
#                                           choices = list("Primary tumor" = "Primary tumor",
#                                                          "Metastasis" = "Metastasis",
#                                                          "Vcap" = "Vcap"
#                                                          )

                                                              
    
    output$choose_inputs <- renderUI({input_list})   # Guardamos todos los inputs en una lisa y lo pasamos a renderUI
                                                     # porque solo se puede llamar 1 vez.
   
    #dibujamos el boxplot con las opciones elegidas
    formulaText <- reactive({
      paste(input$gen,"~",input$result_by)       
   })
    
       # Devuelve formula text para mostrarlo como caption en el gráfico
    output$caption <- renderText({
     formulaText()
    })

   # Generate a plot of the requested variable      
    output$TaylorBoxPlot <- renderPlot({
      if (is.null(input$gen))
        return()
      #browser()
      #View(exp_boxplot())
      #browser()      
      #View(input$gen)
  
     # Generate a summary of the data
     make_summaryQuery <- reactive({
     summary_query<-paste0("select geo_accession,last_update_date,contact_name,platform_id,description,data_processing from taylor_21036,taylor_exp_21036 where geo_accession=id_geoacc and id_gpl like '",input$gen,"' limit 5")
     summary_query_get<-dbGetQuery(connection,summary_query)    
     })
     
     output$summary <- renderTable({      
     make_summaryQuery()
      })
     
     query<-paste0("select ",input$result_by,",expvalue from taylor_21036,taylor_exp_21036 where geo_accession=id_geoacc and id_gpl like '",input$gen,"' order by ",input$result_by,"")                    
     #View(query)   
     exp_boxplot<-dbGetQuery(connection,query)
     class(exp_boxplot)
    df_final_boxplot<-df_nirvana_boxplot.build_df_boxplot(exp_boxplot,group_by=input$result_by)
     #browser()
    
  # Estadísticos descriptivos
  output$descriptive_statistic <- renderPrint({
        summary(df_final_boxplot)
     })
  
  # Análisis t student (t test)
    grupos  <- c("disease_status", "clint_stage", "biopsy_gleason_grade" )
    output$ttest <- renderPrint({
      for (i in grupos){
        df_ttest <- df_nirvana_boxplot.build_df_boxplot(exp_boxplot, i)
        i
        df_ttest  <- df_ttest[,1]
        df_test <- na.exclude(df_ttest)
        # Añadir if para que solo seleccione columnas con más de 1 elemento, pq si no da error.
        t.test(df_test)
      }
    })
  
     par(mar = c(14,3,3,2))
      boxplot_title<-paste(input$gen,"~",input$result_by) 
      boxplot(df_final_boxplot,col="lightblue",main=boxplot_title,las=2)
     
      
    })    
   
    
    #dbDisconnect(connection)
})

