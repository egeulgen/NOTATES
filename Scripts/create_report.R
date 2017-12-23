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
render("Report.Rmd", "pdf_document", params = list(ID=args[1], script_dir=args[2], exome_length=args[3]))