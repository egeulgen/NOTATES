if(!require(rmarkdown))
  install.packages("rmarkdown")
if(!require(markdown))
  install.packages("markdown")
if(!require(formatR))
  install.packages("formatR")
if(!require(knitr))
  install.packages("knitr")
if(!"ggplot2" %in% installed.packages())
  install.packages("ggplot2")
library(rmarkdown)

args <- commandArgs(trailingOnly=TRUE)

## Copy the template RMD file
path2reportRMD <- file.path(paste0("Report_", args[1], ".Rmd"))
file.copy(file.path(args[2], "Report.Rmd"),
          path2reportRMD)

render(path2reportRMD, "pdf_document", params = list(ID=args[1], script_dir=args[2], exome_length=args[3], type=args[4]))

# cleanup
unlink(path2reportRMD)