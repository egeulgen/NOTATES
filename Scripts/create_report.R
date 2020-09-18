#!/usr/bin/env Rscript --vanilla

##################################################
## Project: NOTATES
## Script purpose: Wrapper to create PDF report
## Date: Nov 2, 2019
## Author: Ege Ulgen
##################################################

# Packages and arg.s ------------------------------------------------------
if(!suppressPackageStartupMessages(require(rmarkdown)))
  install.packages("rmarkdown")
if(!suppressPackageStartupMessages(require(knitr)))
  install.packages("knitr")
if(!suppressPackageStartupMessages(require(kableExtra)))
  install.packages("kableExtra")
if(!suppressPackageStartupMessages(require(ggplot2)))
  install.packages("ggplot2")

suppressPackageStartupMessages(library(rmarkdown))
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# Knit PDF summary --------------------------------------------------------
## Copy the template RMD file
summaryRMD_fname <- paste0("Summary_", args[1], ".Rmd")
file.copy(file.path(args[2], "Summary_report.Rmd"),
          summaryRMD_fname, 
          overwrite = TRUE)

## Render report
render(input = summaryRMD_fname, 
       params = list(ID = args[1], 
                     script_dir = path.expand(args[2]), 
                     exome_length = as.numeric(args[3]),
                     type = args[5],
                     primary_tm = as.logical(args[6]),
                     tumor_sample = args[7]),
       output_format = "pdf_document")

# Clean up
unlink(summaryRMD_fname)

# Knit PDF report ---------------------------------------------------------
## Copy the template RMD file
reportRMD_fname <- paste0("Report_", args[1], ".Rmd")
file.copy(file.path(args[2], "Report.Rmd"),
          reportRMD_fname, 
          overwrite = TRUE)

## Render report
render(input = reportRMD_fname, 
       params = list(ID = args[1], 
                     script_dir = path.expand(args[2]), 
                     exome_length = as.numeric(args[3]),
                     exome_bed = args[4],
                     type = args[5],
                     primary_tm = as.logical(args[6]),
                     tumor_sample = args[7]),
       output_format = "pdf_document")

# Clean up
unlink(reportRMD_fname)
