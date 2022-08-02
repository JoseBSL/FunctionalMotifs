
#Upload files to google drive
#install packages
remotes::install_github("claudiozandonella/trackdown",
                        build_vignettes = TRUE, force=T)

install.packages("googledrive")

#set google doc
googledrive::drive_auth()

#load libraries
library(googledrive)
library(trackdown)

#upload file
upload_file(
  file = "Manuscript/Main_text/Main_text_Journal_of_Animal_Ecol.Rmd", 
  gfile = "FunctionalMotifs"
)


#download file 
download_file(
  file = "Manuscript/Main_text/Main_text_Journal_of_Animal_Ecol.Rmd", 
  gfile = "FunctionalMotifs"
)
