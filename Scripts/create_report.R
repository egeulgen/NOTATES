##################################################
## Project: NOTATES
## Script purpose: Wrapper to create PDF report
## Date: Oct 26, 2019
## Author: Ege Ulgen
##################################################

# Packages and arg.s ------------------------------------------------------
if(!suppressMessages(require("rmarkdown")))
  install.packages("rmarkdown")
if(!suppressPackageStartupMessages(require("markdown")))
  install.packages("markdown")
if(!suppressPackageStartupMessages(require(knitr)))
  install.packages("knitr")
if(!suppressPackageStartupMessages(require("ggplot2")))
  install.packages("ggplot2")

suppressPackageStartupMessages(library(rmarkdown))
args <- commandArgs(trailingOnly=TRUE)

# Knit PDF report ---------------------------------------------------------
## Copy the template RMD file
path2reportRMD <- file.path(paste0("Report_", args[1], ".Rmd"))
file.copy(file.path(args[2], "Report.Rmd"),
          path2reportRMD)

render(path2reportRMD, "pdf_document", params = list(ID=args[1], script_dir=args[2], exome_length=args[3], type=args[4]))

# cleanup
unlink(path2reportRMD)