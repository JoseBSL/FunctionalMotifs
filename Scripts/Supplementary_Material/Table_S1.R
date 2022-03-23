##############################
#TABLE S1. List of studies
##############################

#Load libraries
library(kableExtra)
library(dplyr)

#Create data.frame with references

#vector with number of netw per study
`Number of networks` <- c("6", "2", "3", "1", "1", "1", "1", "2", "1", "1", "8", "16", "6", "2", "4", "1", "3", "1")

#`Network structure` <- c(rep("Web", 22), rep("Metaweb",6))

`First author` <- c("Bartomeus", "Dicks", "Dupont", "Elberling", "Fang", "Inouye", "Lundgren", "Olesen", "Small", "Souza", "Kaiser-Bunbury", "Bartomeus", 
                    "Kaiser-Bunbury", "Kaiser-Bunbury", "Peralta", "Burkle", "Arroyo-Correa", "Bundgaard")

Year <- c(2008,2002,2003,1999,2008,1988, 2005,2002,1976,2017,2017,2008,2011,2010,2006,2013,2019,2003)

#Country <- c("Spain", "England", "Denmark", "Sweden", "China", "Australia", "Greenland", "Mauritius and Azores", #"Canada", "Brasil", "Seychelles", "Spain", "Seychelles",
#"Mauritius", "Argentina", "USA", "New Zealand", "Denmark",  "Denmark", "Japan", "Canada", "Venezuela", "Japan", #"Ecuador", "New Zealand", "Venezuela",
#"USA", "Ecuador")

DOI <- c("https://doi.org/10.1007/s00442-007-0946-1", "https://doi.org/10.1046/j.0021-8790.2001.00572.x", "https://doi.org/10.1111/j.1365-2656.2008.01501.x",
         "https://doi.org/10.1111/j.1600-0587.1999.tb00507.x", "https://doi.org/10.1111/1749-4877.12190", "https://doi.org/10.1111/j.1442-9993.1988.tb00968.x",
         "https://doi.org/10.1657/1523-0430(2005)037[0514:TDAHCW]2.0.CO;2", "https://doi.org/10.1046/j.1472-4642.2002.00148.x",
         "/13960/t4km08d21", "https://doi.org/10.1111/1365-2745.12978", "https://doi.org/10.1038/nature21071", "https://github.com/ibartomeus/BeeFunData",
         "https://doi.org/10.1111/j.1365-2745.2010.01732.x", "https://doi.org/10.1016/j.ppees.2009.04.001", "https://doi.org/10.1111/ele.13510",
         "https://doi.org/10.1126/science.1232728", "https://doi.org/10.1111/1365-2745.13332", "Unpublished, Master thesis")


longitude <- c("3.296797",
               "1.575532",
               "9.275556", 
               "18.5", 
               "99.63806", 
               "148.266667",
               "-52", 
               "-31", "57.43",
               "-75.5",
               "-57.885",
               "55.447778",
               "-6.16895",
               "55.433333",
               "57.733333",
               "57.4024",
               "-68.015892",
               "-89.8968771",
               "171.756111",
               "10.233333")

latitude <- c("42.315336",
              "52.762395",
              "56.0725",
              "68.35", 
              "27.90139", 
              "-36.45",
              "71", 
              "39.4", "-20.25",
              "45.4",
              "-21.701111",
              "-4.67",
              "37.234966",
              "-4.666667",
              "-22.7",
              "-20.4538",
              "-32.008985",
              "39.278958",
              "-43.035889",
              "56.066667")

id <- c("Bartomeus1", 
        "Dicks",
        "Dupont", 
        "Elberling", 
        "Fang", 
        "Inouye", 
        "Lundgren", 
        "Olesen", "Olesen",
        "Small",
        "Souza",
        "Kaiser-Bunbury1",
        "Bartomeus2",
        "Kaiser-Bunbury2",
        "Kaiser-Bunbury3",
        "Peralta",
        "Burkle",
        "Arroyo-Correa",
        "Bundgaard")

references <- data.frame(`First author`, Year, `Number of networks`,DOI)

colnames(references) <- c("First author", "Year", "Number of networks",  "DOI")

#Check number of studies
#nrow(references) # 18 studies
#check number of networks
#sum(as.numeric(references$`Number of networks`)) #60

references %>%
  arrange(`First author`)%>%
  kable( longtable = T, booktabs = T,linesep = "\\addlinespace",align = c("cccc")) %>%
  kable_styling(latex_options = c("repeat_header","striped"), font_size = 12, full_width=F,position = c("left"))
