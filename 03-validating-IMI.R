library(readr)
library(dplyr)
library(psych)
library(lavaan)
library(ggraph)
library(semPlot)

library(MVN)
library(daff)
library(robustHD)

wdat <- read_csv('data/WinsorizedIMI.csv')

############################################################################
## Check Assumptions to Reliability Analysis                              ##
############################################################################

png(filename = "report/validation-IMI/univariate-histogram.png", width = 840, height = 840)
(mvn_mod <- mvn(select(wdat, starts_with("Item")), multivariateOutlierMethod = "quan", univariatePlot = "histogram", showOutliers = T))
dev.off()

estimator_cfa <- "WLSMVS" # for non-normal and ordinal data

## kmo factor adequacy
(kmo_mod <- KMO(cor(select(wdat, starts_with("Item")))))

## factorial analysis
factanal(~Item01+Item02+Item06+Item09+Item11+Item12+Item14+Item16+Item17
         +Item18+Item19+Item21+Item22+Item23+Item24+Item25+Item38+Item26
         +Item27+Item39+Item28+Item29
         , factors=4, data=wdat, rotation = "varimax")

factanal(~Item21+Item09+Item19+Item38+Item12+Item22
         +Item06+Item17+Item23
         +Item28+Item27+Item25+Item29+Item26
         +Item14+Item01+Item16+Item18+Item11
         +Item02+Item24+Item39
         , factors=4, data=wdat, rotation = "varimax")

# removing items with loading less than < 0.4
factanal(~Item21+Item09+Item19+Item38+Item12+Item22
         +Item06+Item17+Item23
         +Item28+Item27+Item25+Item29+Item26
         +Item14+Item01+Item16+Item18+Item11
         +Item24+Item39
         , factors=4, data=wdat, rotation = "varimax")

# removing crosloading items with loading less than < 0.2
factanal(~Item21+Item09+Item19+Item38+Item12+Item22
         +Item28+Item27+Item25+Item29+Item26
         +Item14+Item01+Item16+Item18+Item11
         +Item06+Item17+Item23
         , factors=4, data=wdat, rotation = "varimax")

# first eingenvalue between 3 to 5 then
# we removed items that don't fit by meaning based on consult with psychometrics
factanal(~Item21+Item09+Item19+Item38+Item12+Item22 # IE: Interest/Enjoyment
         +Item28+Item27+Item25+Item29+Item26 # PC2: Perceived Competence
         +Item14+Item01+Item16+Item18+Item11 # PT: Pressure/Tension
         +Item06+Item17+Item23 # PC: Perceived Choice
         , factors=4, data=wdat, rotation = "varimax")

(fa_mod <- fa(select( 
  wdat
  , starts_with("Item21"), starts_with("Item09"), starts_with("Item19"), starts_with("Item38"), starts_with("Item12"), starts_with("Item22")
  , starts_with("Item28"), starts_with("Item27"), starts_with("Item25"), starts_with("Item29"), starts_with("Item26")
  , starts_with("Item14"), starts_with("Item01"), starts_with("Item16"), starts_with("Item18"), starts_with("Item11")
  , starts_with("Item06"), starts_with("Item17"), starts_with("Item23")
), nfactors = 4, rotate = "varimax"))

## validating models where orthogonal means no correlation between factors

multi_mdl <- '
IE =~ Item21 + Item09 + Item19 + Item38 + Item12 + Item22
PC2 =~ Item28 + Item27 + Item25 + Item29 + Item26
PT =~ Item14 + Item01 + Item16 + Item18 + Item11
PC =~ Item06 + Item17 + Item23
'
second_mdl <- '
IE =~ Item21 + Item09 + Item19 + Item38 + Item12 + Item22
PC2 =~ Item28 + Item27 + Item25 + Item29 + Item26
PT =~ Item14 + Item01 + Item16 + Item18 + Item11
PC =~ Item06 + Item17 + Item23
IM =~ NA*IE + PC2 + PT + PC
IM ~~ 1*IM
'
bifactor_mdl <- '
IE =~ a*Item21 + a*Item09 + a*Item19 + a*Item38 + a*Item12 + a*Item22
PC2 =~ b*Item28 + b*Item27 + b*Item25 + b*Item29 + b*Item26
PT =~ c*Item14 + c*Item01 + c*Item16 + c*Item18 + c*Item11
PC =~ d*Item06 + d*Item17 + d*Item23
IM =~ Item21 + Item09 + Item19 + Item38 + Item12 + Item22 +
Item28 + Item27 + Item25 + Item29 + Item26 +
Item14 + Item01 + Item16 + Item18 + Item11 +
Item06 + Item17 + Item23
'

(fitMeasures_df <- do.call(rbind, lapply(
  list(
    "Global sample" = list(
      dat = wdat,
      mdls = list(
        "Multidimensional model" = list(mdl = multi_mdl, plotFile = "multidimensional-model.png")
        , "Second-order model" = list(dat = wdat, mdl = second_mdl, plotFile = "second-order-factor-model.png")
        , "Bi-factor model" = list(dat = wdat, mdl = bifactor_mdl, plotFile = "bi-factor-model.png")))
  )
  , FUN = function(s) {
    fit_df <- do.call(rbind, lapply(
      s$mdls
      , FUN = function(x) {
        return(get_fitMeasures(dat = s$dat, mdl = x$mdl, estimator = estimator_cfa))
      }
    ))
    return(rbind(c(NA), fit_df))
  }))
)

# select second-order model to measure intrinsic motivation
(cfa_mod <- cfa(second_mdl, data = wdat, std.lv = T, estimator = estimator_cfa))


############################################################################
## Reliability Analysis Using Cronbach's alpha                            ##
############################################################################

rdat <- select(wdat, starts_with("UserID"))

rdat["Item21IE"] <- wdat["Item21"]
rdat["Item09IE"] <- wdat["Item09"]
rdat["Item19IE"] <- wdat["Item19"]
rdat["Item38IE"] <- wdat["Item38"]
rdat["Item12IE"] <- wdat["Item12"]
rdat["Item22IE"] <- wdat["Item22"]

rdat["Item28PC2"] <- wdat["Item28"]
rdat["Item27PC2"] <- wdat["Item27"]
rdat["Item25PC2"] <- wdat["Item25"]
rdat["Item29PC2"] <- wdat["Item29"]
rdat["Item26PC2"] <- wdat["Item26"]

rdat["Item14PT"] <- wdat["Item14"]
rdat["Item01PT"] <- wdat["Item01"]
rdat["Item16PT"] <- wdat["Item16"]
rdat["Item18PT"] <- wdat["Item18"]
rdat["Item11PT"] <- wdat["Item11"]

rdat["Item06PC"] <- wdat["Item06"]
rdat["Item17PC"] <- wdat["Item17"]
rdat["Item23PC"] <- wdat["Item23"]

rdat <- dplyr::mutate(
  rdat
  , `Interest/Enjoyment` = (rdat$Item21IE+rdat$Item09IE+rdat$Item38IE+rdat$Item12IE+rdat$Item22IE+8-rdat$Item19IE)/6
  , `Perceived Competence` = (rdat$Item28PC2+rdat$Item27PC2+rdat$Item25PC2+rdat$Item29PC2+rdat$Item26PC2)/5
  , `Pressure/Tension` = (rdat$Item14PT+rdat$Item16PT+rdat$Item18PT+16-(rdat$Item01PT+rdat$Item11PT))/5
  , `Perceived Choice` = (rdat$Item23PC+16-(rdat$Item06PC+rdat$Item17PC))/3
)
rdat <- dplyr::mutate(
  rdat
  , `Intrinsic Motivation` = (rdat$`Interest/Enjoyment`
                              +rdat$`Perceived Competence`
                              +rdat$`Perceived Choice`
                              +8-(rdat$`Pressure/Tension`))/4
)
if (!file.exists('data/IMI.csv')) {
  write_csv(rdat, path = 'data/IMI.csv')
}

alpha_mods <- lapply(
  list(
    "IM" = list(
      id = "IM",
      dat = select(rdat, starts_with("Item")),
      lbl = "Intrinsic Motivation",
      inv_keys = c("Item19IE"
                   , "Item14PT", "Item16PT", "Item18PT"
                   , "Item06PC", "Item17PC")),
    "IE" = list(
      id = "IE",
      dat = select(rdat, ends_with("IE")),
      lbl = "Interest/Enjoyment",
      inv_keys = c("Item19IE")),
    "PC2" = list(
      id = "PC2",
      dat = select(rdat, ends_with("PC2")),
      lbl = "Perceived Competence",
      inv_keys = c()),
    "PT" = list(
      id = "PT",
      dat = select(rdat, ends_with("PT")),
      lbl = "Pressure/Tension",
      inv_keys = c("Item01PT", "Item11PT")),
    "PC" = list(
      id = "PC",
      dat = select(rdat, ends_with("PC")),
      lbl = "Perceived Choice",
      inv_keys = c("Item06PC", "Item17PC"))
    )
  , FUN = function(x) {
    alpha_mod <- psych::alpha(select(x$dat, starts_with("Item")), keys = x$inv_keys)
    
    cat("\n... ", x$lbl, " ...\n"); summary(alpha_mod)
    
    return(list(id = x$id, lbl = x$lbl, all = alpha_mod))
  })

## Write results in an Excel Workbook
if (!file.exists("report/validation-IMI/RelAnalysis.xlsx")) {
  filename <- "report/validation-IMI/RelAnalysis.xlsx"
  wb <- createWorkbook(type="xlsx")
  
  write_kmo_in_workbook(kmo_mod, wb)
  write_fa_in_workbook(fa_mod, wb)
  lapply(alpha_mods, FUN = function(mod){
    write_alpha_in_workbook(mod$all, wb, mod$lbl, mod$id)
  })
  
  xlsx::saveWorkbook(wb, filename)
}

# Export summaries in latex format
write_cfa_model_fits_in_latex(
  fitMeasures_df
  , in_title = "in the validation of the adapted Portuguese IMI"
  , filename = "report/latex/cfa-model-fit.tex")

filename <- "report/latex/reliability-analysis.tex"
if (!file.exists(filename)) {
  write_rel_analysis_in_latex(
    fa_mod, cfa_mod, alpha_mods
    , in_title = "adapted Portuguese IMI"
    , filename = filename
    , key_labels = list('Global'='all')
    , robust = T
  )
}

