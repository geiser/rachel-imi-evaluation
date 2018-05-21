wants <- c('readr', 'dplyr', 'devtools','readxl')
has <- wants %in% rownames(installed.packages())
if (any(!has)) install.packages(wants[!has])
if (!any(rownames(installed.packages()) %in% c('careless'))) {
  devtools::install_github('ryentes/careless')
}

library(daff)
library(readr)
library(dplyr)
library(careless)
library(readxl)

SourceIMI <- read_excel("data/IMI-rachel.xlsx", sheet = "IMI_Todos_Alunos")

##########################################################################
## Removing Careless in Motivation Survey                               ##
##########################################################################

resp <- select(SourceIMI, starts_with("UserID"), starts_with("Item"))

(careless_info <- careless::longstring(select(resp, -starts_with("UserID")), na=T))
respIMI <- resp[careless_info <= 12 & complete.cases(resp),]

filename <- 'data/SourceIMI.csv'
if (!file.exists(filename)) {
  write_csv(respIMI, filename)
}

## write in latex
render_diff(ddIMI <- diff_data(resp, respIMI))
filename <- 'report/latex/careless-IMI.tex'
write_careless_in_latex(
  ddIMI, filename, in_title = "in the IMI data collected over the pilot empirical study")

