if(!"rmarkdown" %in% installed.packages())
  install.packages("rmarkdown")
if(!"markdown" %in% installed.packages())
  install.packages("markdown")
if(!"formatR" %in% installed.packages())
  install.packages("formatR")
if(!"knitr" %in% installed.packages())
  install.packages("knitr")
if(!"pander" %in% installed.packages())
  install.packages("pander")

library(rmarkdown)

args <- commandArgs(trailingOnly=TRUE)
patientID <- args[1]
sdir <- args[2]
  
render("Report.Rmd", "pdf_document", params = list(ID=patientID, script_dir=sdir))
