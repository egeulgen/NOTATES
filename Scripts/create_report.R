if(!require(rmarkdown))
  install.packages("rmarkdown")
if(!require(markdown))
  install.packages("markdown")
if(!require(formatR))
  install.packages("formatR")
if(!require(knitr))
  install.packages("knitr")

library(rmarkdown)

args <- commandArgs(trailingOnly=TRUE)
patientID <- args[1]
sdir <- args[2]
  
render("Report.Rmd", "pdf_document", params = list(ID=patientID, script_dir=sdir))
