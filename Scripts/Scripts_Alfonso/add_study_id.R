
add_study_id <- function(dataframe_input){
  
  col_names_dataframe <- colnames(dataframe_input)
  
  dataframe_input$study_id <- dataframe_input$Network_id
  
  dataframe_input$study_id <- as.character(dataframe_input$study_id)
  
  dataframe_input$study_id[grepl("peralta_2006", dataframe_input$study_id)] <- "peralta_2006"
  dataframe_input$study_id[grepl("small_1976", dataframe_input$study_id)] <- "small_1976"
  dataframe_input$study_id[grepl("arroyo_correa_new_zealand", dataframe_input$study_id)] <- "arroyo_correa_2019"
  dataframe_input$study_id[grepl("fang_2008", dataframe_input$study_id)] <- "fang_2008"
  dataframe_input$study_id[grepl("kaiser-bunbury_2017", dataframe_input$study_id)] <- "kaiser-bunbury_2017"
  dataframe_input$study_id[grepl("inouye_1988", dataframe_input$study_id)] <- "inouye_1988"
  dataframe_input$study_id[grepl("kaiser-bunbury_2010", dataframe_input$study_id)] <- "kaiser-bunbury_2010"
  dataframe_input$study_id[grepl("kaiser-bunbury_2011", dataframe_input$study_id)] <- "kaiser-bunbury_2011"
  dataframe_input$study_id[grepl("burkle_usa_2013", dataframe_input$study_id)] <- "burkle_2013"
  dataframe_input$study_id[grepl("dicks_2002", dataframe_input$study_id)] <- "dicks_2002"
  dataframe_input$study_id[grepl("dupont_2009", dataframe_input$study_id)] <- "dupont_2009"
  dataframe_input$study_id[grepl("bartomeus_spain_2008_medca", dataframe_input$study_id)] <- "bartomeus_2008"
  dataframe_input$study_id[grepl("bartomeus_spain_2008_batca", dataframe_input$study_id)] <- "bartomeus_2008"
  dataframe_input$study_id[grepl("bartomeus_spain_2008_selop", dataframe_input$study_id)] <- "bartomeus_2008"
  dataframe_input$study_id[grepl("bartomeus_spain_2008_miqop", dataframe_input$study_id)] <- "bartomeus_2008"
  dataframe_input$study_id[grepl("bartomeus_spain_2008_fraop", dataframe_input$study_id)] <- "bartomeus_2008"
  dataframe_input$study_id[grepl("lundgren_2005", dataframe_input$study_id)] <- "lundgren_2005"
  dataframe_input$study_id[grepl("olesen_2002_mauritius", dataframe_input$study_id)] <- "olesen_2002_mauritius"
  dataframe_input$study_id[grepl("olesen_2002_azores", dataframe_input$study_id)] <- "olesen_2002_azores"
  dataframe_input$study_id[grepl("bartomeus_spain_2008", dataframe_input$study_id)] <- "bartomeus_spain_2008"
  dataframe_input$study_id[grepl("bundgaard_2003_denmark", dataframe_input$study_id)] <- "bundgaard_2003"
  dataframe_input$study_id[grepl("elberling_sweeden_1999", dataframe_input$study_id)] <- "elberling_1999"
  
  dataframe_input <- dataframe_input[,c("study_id",col_names_dataframe)]
  
  return(dataframe_input)
  
}

