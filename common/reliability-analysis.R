
wants <- c('lavaan', 'xlsx', 'reshape', 'rJava', 'psych', 'dplyr', 'readr', 'readxl', 'Hmisc', 'devtools', 'MVN')
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
#library(devtools)
#install_github("kassambara/r2excel")
# FIX for (Mac)OS X Sierra:
#  sudo ln -f -s $(/usr/libexec/java_home)/lib/server/libjvm.dylib /usr/local/lib

library(dplyr)
library(readr)
library(psych)
library(lavaan)
library(r2excel)

## get fit measures as data.frame
get_fitMeasures <- function(dat, mdl, estimator = "ML", plotFile = NULL) {
  fit <-  cfa(mdl, data = dat, std.lv = T, estimator = estimator)
  
  if (!is.null(plotFile)) {
    png(filename = plotFile, width = 420, height = 960)
    semPaths(fit, layout = "tree2", rotation = 2, curvePivot = T, intercepts = F, residuals = F, reorder = T)
    dev.off()
  }
  
  return(data.frame(
    "df" = tryCatch(fitMeasures(fit, "df.scaled"), error = function(e) NA)
    , "chisq" = tryCatch(fitMeasures(fit, "chisq.scaled"), error = function(e) NA)
    , "AGFI" = tryCatch(fitMeasures(fit, "agfi"), error = function(e) NA)
    , "TLI" = tryCatch(fitMeasures(fit, "tli.scaled"), error = function(e) NA)
    , "CFI" = tryCatch(fitMeasures(fit, "cfi.scaled"), error = function(e) NA)
    , "RMSEA" = tryCatch(fitMeasures(fit, "rmsea.scaled"), error = function(e) NA)
    , "CI.lwr" = tryCatch(fitMeasures(fit, "rmsea.ci.lower.scaled"), error = function(e) NA)
    , "CI.upr" = tryCatch(fitMeasures(fit, "rmsea.ci.upper.scaled"), error = function(e) NA)
  ))
}

## function to write kmo module as sheet
write_kmo_in_workbook <- function(kmo_mod, wb) {
  library(r2excel)
  
  sheet <- xlsx::createSheet(wb, sheetName = "MSA")
  xlsx.addHeader(wb, sheet, "Measure Sampling Adequacy Using KMO Factor Adequacy", startCol = 1)
  
  xlsx.addLineBreak(sheet, 2)
  adequacy <- 'unacceptable'
  if (kmo_mod$MSA >= 0.5 && kmo_mod$MSA < 0.6) adequacy <- 'miserable'
  if (kmo_mod$MSA >= 0.6 && kmo_mod$MSA < 0.7) adequacy <- 'mediocre'
  if (kmo_mod$MSA >= 0.7 && kmo_mod$MSA < 0.8) adequacy <- 'middling'
  if (kmo_mod$MSA >= 0.8 && kmo_mod$MSA < 0.9) adequacy <- 'meritorious'
  if (kmo_mod$MSA >= 0.9) adequacy <- 'marvelous'
  xlsx.addTable(wb, sheet, data.frame("OverallMSA"=kmo_mod$MSA, adequacy = adequacy), startCol = 1, row.names = F)
  
  
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "MSA for each item", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet,  t(data.frame(kmo_mod$MSAi)), startCol = 1, row.names = F)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Image matrix", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, kmo_mod$Image, startCol = 1)
}

## function to write cfa module as sheet
write_cfa_in_workbook <- function(cfa_mod, wb) {
  library(r2excel)
  
  sheet <- xlsx::createSheet(wb, sheetName = "CFA")
  xlsx.addHeader(wb, sheet, "Confirmatory Factor Analysis Using Minimum Residual", startCol = 1)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Fit Measures for Model", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, t(data.frame(fitMeasures(cfa_mod))), startCol = 1, row.names = F)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Parameter Estimates", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(parameterEstimates(cfa_mod, standardized = T, fmi = T)), startCol = 1, row.names = F)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized Model Parameters (lambda)", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(lavInspect(cfa_mod, "std")$lambda), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized Model Parameters (theta)", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(lavInspect(cfa_mod, "std")$theta), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized Model Parameters (psi)", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(lavInspect(cfa_mod, "std")$psi), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized Model Parameters (nu)", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(lavInspect(cfa_mod, "std")$nu), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized Model Parameters (alpha)", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(lavInspect(cfa_mod, "std")$alpha), startCol = 1, row.names = T)
}

## function to write efa module as sheet
write_fa_in_workbook <- function(fa_mod, wb) {
  library(r2excel)
  
  sheet <- xlsx::createSheet(wb, sheetName = "EFA")
  xlsx.addHeader(wb, sheet, "Exploratory Factor Analysis Using Minimum Residual", startCol = 1)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Standardized loadings", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(unclass(fa_mod$loadings)), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Measures of factor score adequacy", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(print(fa_mod)), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Item complexity", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(Value=fa_mod$complexity), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Extra information", level = 2, startCol = 1)
  ext_info <- paste0('Mean item complexity = ', mean(fa_mod$complexity), '\n')
  ext_info <- paste0(ext_info, 'The degrees of freedom for the model are = ', fa_mod$dof, '\n')
  ext_info <- paste0(ext_info, 'The objective function was = ', fa_mod$objective, '\n')
  ext_info <- paste0(ext_info, 'The Chi Square of the model is = ', fa_mod$chi, '\n')
  ext_info <- paste0(ext_info, 'The root mean square of the residuals (RMSR) is = ', fa_mod$rms, '\n')
  ext_info <- paste0(ext_info, 'Tucker Lewis Index of factoring reliability is =  ', fa_mod$TLI, '\n')
  ext_info <- paste0(ext_info, 'RMSEA index =  ', fa_mod$RMSEA[[1]], '\n')
  ext_info <- paste0(ext_info, 'The 90% confidence intervals of RMSEA are lower = ', fa_mod$RMSEA[[2]], ' and upper = ', fa_mod$RMSEA[[3]],'\n')
  ext_info <- paste0(ext_info, 'BIC = ', fa_mod$BIC, '\n')
  xlsx.addParagraph(wb, sheet, value = ext_info, startCol = 1)
}

## function to write alpha module as sheet
write_alpha_in_workbook <- function(alpha_mod, wb, name, short_name = NULL) {
  library(r2excel)
  
  if (is.null(short_name)) short_name <- paste0(" (", substr(name, 1, 12), ")")
  
  sheetName <- paste0("RelAnlys", ifelse(!is.null(name), short_name, ""))
  sheet <- xlsx::createSheet(wb, sheetName = sheetName)
  xlsx.addHeader(wb, sheet, paste0("Reliability Analysis", ifelse(!is.null(name), paste0(" of ", name))), startCol = 1)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Summary", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, as.data.frame(summary(alpha_mod)), startCol = 1, row.names = F)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Reliability if an item is dropped", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(alpha_mod$alpha.drop), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Item statistics", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(alpha_mod$item.stats), startCol = 1, row.names = T)
  
  xlsx.addLineBreak(sheet, 2)
  xlsx.addHeader(wb, sheet, "Frequency of each item response", level = 2, startCol = 1)
  xlsx.addTable(wb, sheet, data.frame(alpha_mod$response.freq), startCol = 1, row.names = T)
  
}
