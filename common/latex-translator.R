

## Function to obtain a dataframe with the summary of descriptive and post-hoc 
## information from a parametric statistic result
get_descritive_and_post_hoc_dataframe <- function(parametric_result) {
  library(Hmisc)
  
  t_test_pairs <- parametric_result$t.test.pairs
  constrasts_set <- parametric_result$post.hoc$contrasts
  lsmeans_set <- parametric_result$post.hoc$lsmeans
  tukey_set <- parametric_result$post.hoc$tukey
  
  return (do.call(
    rbind,
    lapply(names(constrasts_set), FUN = function(src) {
      constrasts <- constrasts_set[[src]]
      lsmeans <- lsmeans_set[[src]]
      tukeys <- tukey_set[[src]]
      t_test_pair <- t_test_pairs[[src]]
      
      #for (i in 1:nrow(constrasts)) {
      to_return <- do.call(rbind, lapply(c(1:nrow(constrasts)), FUN = function (i) {
        constr_info <- as.list(constrasts[i,])
        
        constr_keys <- strsplit(as.character(constr_info$contrast), ' - ')[[1]]
        pairs <- as.vector(sapply(constr_keys, FUN = function(x) {
          columns_to_pair <- colnames(constrasts)[
            !colnames(constrasts) %in% c('contrast','estimate','SE','df','t.ratio','p.value')]
          if (!is.null(columns_to_pair) && !is.na(columns_to_pair) && length(columns_to_pair) > 0) {
            return(paste(na.omit(
              c(x, as.character(constrasts[i, columns_to_pair]))), collapse = '.'))
          } else return(x)
        }))
        
        ## get tukey info
        tukey_inv <- F
        tukey_key <- gsub("\\.",":", paste0(pairs, collapse = "-"))
        if (!any(rownames(tukeys) %in% tukey_key)) {
          tukey_inv <- T
          tukey_key <- gsub("\\.",":", paste0(rev(pairs), collapse = "-"))  
        }
        tukey_info <- as.list(tukeys[tukey_key,])
        if (tukey_inv) {
          tukey_info$diff <- -1*tukey_info$diff
          tukey_info_lwr <- tukey_info$lwr
          tukey_info$lwr <- -1*tukey_info$upr
          tukey_info$upr <- -1*tukey_info_lwr
        }
        
        ## get t.test info
        t_test_key <- paste0(pairs, collapse = ':')
        if (!any(names(t_test_pair) %in% t_test_key)) {
          t_test_key <- paste0(rev(pairs), collapse = ':')
        }
        t_test_info <- t_test_pair[[t_test_key]]$two.sided$result
        
        
        ## get descriptive part from lsmeans
        to_add <- lsmeans[lsmeans$Pairs %in% as.vector(pairs),]
        to_add <- dplyr::mutate(to_add, t.ratio = NA, p.value = NA, p.ajd = NA, g = NA)
        to_add <- merge(to_add, t_test_info, by.x = "Pairs", by.y = "Group")
        to_add <- to_add[,c('Pairs','N','Mean','lsmean','SE','df.x','lower.CL','upper.CL')]
        colnames(to_add) <- c('Group','N','mean','lsmean','SE', 'df','lwr.CI','upr.CI')
        
        ## obtaining the data result
        result_df <- plyr::rbind.fill(to_add, data.frame(
          "Group" = paste0(pairs, collapse = " - ")
          , "N" = sum(to_add$N)
          , "mean" = tukey_info$diff
          , "lsmean" = constr_info$estimate
          , "SE" = constr_info$SE
          , "df" = NA
          , "lwr.CI" = tukey_info$lwr
          , "upr.CI" = tukey_info$upr
          , "t.ratio" = constr_info$t.ratio
          , "p.value" = constr_info$p.value
          , "p.ajd" = tukey_info$`p adj`
          , "g" = t_test_info$g[[1]]
        ))
        row_names <- result_df$Group
        
        result_df <- dplyr::mutate(
          round(result_df[,-1], 3)
          , sig=ifelse(result_df$p.ajd > 0.05, NA
                       , ifelse(result_df$p.ajd > 0.01, '*', '**'))
          , mag=ifelse((!is.null(result_df$p.ajd)
                        & !is.na(result_df$p.ajd)
                        & result_df$p.ajd <= 0.05)
                       , as.character(t_test_info$magnitude[[1]]), NA))
        rownames(result_df) <- row_names
        return(result_df)
      }
      ))
      
      return(to_return)
    })
  )
  )
}


write_kruskal_statistics_analysis_in_latex <- function(
  nonparametric_results, info, filename, in_title = NULL, append = F) {
  
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
  }
  
  write(paste0("\\section{Summaries of Nonparametric Statistics Analysis"
               , in_title, "}"), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  list_info <- as.list(names(info))
  names(list_info) <- names(info)
  
  result_df <- do.call(rbind, lapply(list_info, FUN = function(name) {
    n_result <- nonparametric_results[[name]]
    mod_df <- round(n_result$mod.df, 3)
    
    Sig <- sapply(mod_df$p.value, FUN = function(x) {
      return(ifelse(x > 0.05, NA, ifelse(x >0.01, '*', '**')))
    })
    
    mod_df <- cbind(mod_df, Sig)
    #mod_df <- rbind(c(NA), mod_df)
    return(mod_df)
  }))
  
  ##
  latex(result_df
        , caption = paste("Summary of Kruskal-Wallis rank test results", in_title)
        #, insert.bottom = 
        , size = "small", longtable = T, ctable=F, landscape = F
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(
    c("\\raggedleft{", "\\scriptsize{"
      , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
      , "}}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
  ## wilcoxon dataframe test
  result_wilcoxon_df <- do.call(rbind, lapply(list_info, FUN = function(name) {
    n_result <- nonparametric_results[[name]]
    
    set_wilcox_mods_df <- do.call(rbind, lapply(
      n_result$wilcox.pairs, FUN = function(wilcox_mods){
        wilcox_mods_df <- do.call(rbind, lapply(
          wilcox_mods, FUN = function(mods) {
            mods$dat <- NULL
            mods_df <- do.call(rbind, lapply(
              mods, FUN = function(mod) {
                if (any(mod$result$p.value <= 0.05)) {
                  return(mod$result)
                } else {
                  return(NULL)
                }
              }))
            return(mods_df)
          }))
        return(wilcox_mods_df)
      }))
    return(set_wilcox_mods_df)
  }))
  
  result_wilcoxon_df$Median <- round(result_wilcoxon_df$Median, 2)
  result_wilcoxon_df$Mean.Ranks <- round(result_wilcoxon_df$Mean.Ranks, 2)
  result_wilcoxon_df$Sum.Ranks <- round(result_wilcoxon_df$Sum.Ranks, 2)
  result_wilcoxon_df$Z <- round(result_wilcoxon_df$Z, 2)
  result_wilcoxon_df$p.value <- round(result_wilcoxon_df$p.value, 3)
  result_wilcoxon_df$r <- round(result_wilcoxon_df$r,3)
  
  write(paste0("\\section{Wilcoxon Pairs Statistics Analysis"
               , in_title, "}"), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  latex(result_wilcoxon_df
        , caption = paste("Descriptive statistic of the pair wilcoxon analysis", in_title)
        #, insert.bottom = 
        , size = "small", longtable = T, ctable=F, landscape = T
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(
    c("\\raggedleft{", "\\scriptsize{"
      , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
      , "}}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}

##################################################################


## Function to write rel_analysis as latex
write_rel_analysis_in_latex <- function(
  fa_mod, cfa_mod, alpha_mods, in_title = "none", filename = "RelAnalysis.tex", key_labels = list(), robust = F, append = F) {
  
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
    write("\\section{Summaries of Reliability Analysis}", file = filename, append = T)
  }
  
  # summary of exploratory/confirmatory factor anlaysis
  write("\\subsection{Exploratory and Confirmatory Factor Analysis}", file = filename, append = T)
  write("", file = filename, append = T)
  
  fa_summary_df <- round(data.frame(unclass(fa_mod$loadings)), 3)
  fa_aux_df <- round(data.frame(print(fa_mod))[c(1,3,4),], 3)
  colnames(fa_aux_df) <- colnames(fa_summary_df)
  fa_summary_df <- rbind(fa_summary_df, fa_aux_df)
  
  fit_measures <- as.list(round(fitMeasures(cfa_mod, fit.measures = "all"), 3))
  cfa_summary <- paste0(
    c("\\begin{flushright}{", "\\scriptsize{"
      , paste0(c(paste0("CFI: ", if(robust) fit_measures$cfi.scaled else fit_measures$cfi)
                 , paste0("TLI: ", if(robust) fit_measures$tli.scaled else fit_measures$tli)
                 , paste0("df: ", if(robust) fit_measures$df.scaled else fit_measures$df)
                 , paste0("$chi^2$: ", if(robust) fit_measures$chisq.scaled else fit_measures$chisq)
                 , paste0("p-value: ", if(robust) fit_measures$pvalue.scaled else fit_measures$pvalue)
                 , paste0("RMSEA: ", paste0(if(robust) fit_measures$rmsea.scaled else fit_measures$rmsea, " "
                                                     , "[", if(robust) fit_measures$rmsea.ci.lower.scaled else fit_measures$rmsea.ci.lower
                                                     , ", ", if(robust) fit_measures$rmsea.ci.upper.scaled else fit_measures$rmsea.ci.upper, "]"))
                 , paste0("p-value(RMSEA): ", fit_measures$rmsea.pvalue))
               , collapse = "; ")
      , "}}\\end{flushright}", "\n"), collapse = " ")
  latex(fa_summary_df
        , caption = paste("Summary of exploratory and confirmatory factor analysis for the", in_title)
        , size = "small", longtable = T, ctable=F, landscape = F
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(cfa_summary, file = filename, append = T)
  write("", file = filename, append = T)
  
  # summary of reliability analysis
  write("\\subsection{Reliability Analysis}", file = filename, append = T)
  write("", file = filename, append = T)
  
  rel_summary_df <- data.frame(lapply(key_labels, FUN = function(k) {
    return(round(unname(sapply(alpha_mods, FUN = function(x) {
      return(x[[k]]$total$std.alpha)
    })), 3))
  }))
  row.names(rel_summary_df) <- unname(sapply(alpha_mods, FUN = function(x) {
    return(x$lbl)
  }))
  
  latex(rel_summary_df
        , caption = paste("Result of reliability analysis for the", in_title)
        , size = "small", longtable = T, ctable=F, landscape = F
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write("", file = filename, append = T)
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}

## Function to write in latex a param statistic analysis
write_param_statistics_analysis_in_latex <- function(
  parametric_results, dvs, filename, in_title = NULL, append = F) {
  
  library(Hmisc)
  
  write(paste0("\\section{Summaries of Parametric Statistics Analysis"
               , in_title, "}"), file = filename, append = append)
  write("", file = filename, append = T)
  
  ##
  result_df <- do.call(rbind, lapply(dvs, FUN = function(dv) {
    p_result <- parametric_results[[dv]]
    aov_df <- round(p_result$ezAov$Anova, 3)
    
    Sig <- sapply(aov_df$`Pr(>F)`, FUN = function(x) {
      return(ifelse(x > 0.05, NA, ifelse(x >0.01, '*', '**')))
    })
    Sig[[1]] <- NA
    
    aov_df <- cbind(aov_df, Sig)
    aov_df <- rbind(c(NA), aov_df)
    return(aov_df)
  }))
  rownames(result_df) <- NULL
  
  rownames(result_df) <- c(unlist(lapply(dvs, FUN = function(dv) {
    p_result <- parametric_results[[dv]]
    #sapply(rownames(p_result$ezAov$Anova), paste, dv, sep=".")
    return(c(dv, paste(dv, rownames(p_result$ezAov$Anova), sep = ":")))
  })))
  
  ##
  latex(result_df
        , caption = paste("Summary of two-way ANOVA results", in_title)
        , size = "small", longtable = T, ctable=F, landscape = F
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(c("\\begin{flushright}{", "\\scriptsize{"
                 , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
                 , "}}\\end{flushright}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  post_hoc_df <- do.call(rbind, lapply(dvs, FUN = function(dv) {
    parametric_result <- parametric_results[[dv]]
    result_df <- get_descritive_and_post_hoc_dataframe(parametric_result)
    row_names <- rownames(result_df)
    result_df <- rbind(rep(NA, ncol(result_df)), result_df)
    rownames(result_df) <- c(dv, paste(dv, row_names, sep=":"))
    return(result_df)
  }))
  
  latex(post_hoc_df
        , caption = paste("Descriptive statistics and Tukey post-hoc test results", in_title)
        , size = "small", longtable = T, ctable=F, landscape = T
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(c("\\begin{flushright}{", "\\scriptsize{"
                 , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
                 , "}}\\end{flushright}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
}

write_nonparam_statistics_analysis_in_latex <- function(
  nonparametric_results, dvs, filename, in_title = NULL, append = F) {
  
  library(Hmisc)
  
  write(paste0("\\section{Summaries of Nonparametric Statistics Analysis"
               , in_title, "}"), file = filename, append = append)
  write("", file = filename, append = append)
  
  ##
  result_df <- do.call(rbind, lapply(dvs, FUN = function(dv) {
    n_result <- nonparametric_results[[dv]]
    sch_df <- round(n_result$sch, 3)
    
    Sig <- sapply(sch_df$p.value, FUN = function(x) {
      return(ifelse(x > 0.05, NA, ifelse(x >0.01, '*', '**')))
    })
    
    sch_df <- cbind(sch_df, Sig)
    sch_df <- rbind(c(NA), sch_df)
    return(sch_df)
  }))
  rownames(result_df) <- NULL
  
  rownames(result_df) <- c(unlist(lapply(dvs, FUN = function(dv) {
    n_result <- nonparametric_results[[dv]]
    return(c(dv, paste(dv, rownames(n_result$sch), sep = ":")))
  })))
  
  ##
  latex(result_df
        , caption = paste("Summary of Scheirer-Ray-Hare results", in_title)
        , size = "small", longtable = T, ctable=F, landscape = F
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(c("\\begin{flushright}{", "\\scriptsize{"
                 , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
                 , "}}\\end{flushright}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
  ## wilcoxon dataframe test
  list_dvs <- as.list(dvs)
  names(list_dvs) <- dvs
  result_wilcoxon_df <- do.call(rbind, lapply(list_dvs, FUN = function(dv) {
    n_result <- nonparametric_results[[dv]]
    
    set_wilcox_mods_df <- do.call(rbind, lapply(
      n_result$wilcox.pairs, FUN = function(wilcox_mods){
      wilcox_mods_df <- do.call(rbind, lapply(
        wilcox_mods, FUN = function(mods) {
          mods$dat <- NULL
          mods_df <- do.call(rbind, lapply(
            mods, FUN = function(mod) {
              if (any(mod$result$p.value <= 0.05)) {
                return(mod$result)
              } else {
                return(NULL)
              }
            }))
          return(mods_df)
        }))
      return(wilcox_mods_df)
    }))
    return(set_wilcox_mods_df)
  }))
  
  result_wilcoxon_df$Median <- round(result_wilcoxon_df$Median, 2)
  result_wilcoxon_df$Mean.Ranks <- round(result_wilcoxon_df$Mean.Ranks, 2)
  result_wilcoxon_df$Sum.Ranks <- round(result_wilcoxon_df$Sum.Ranks, 2)
  result_wilcoxon_df$Z <- round(result_wilcoxon_df$Z, 2)
  result_wilcoxon_df$p.value <- round(result_wilcoxon_df$p.value, 3)
  result_wilcoxon_df$r <- round(result_wilcoxon_df$r,3)
  
  write(paste0("\\section{Wilcoxon Pairs Statistics Analysis"
               , in_title, "}"), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  latex(result_wilcoxon_df
        , caption = paste("Descripteive statistic of the pair wilcoxon analysis", in_title)
        , size = "small", longtable = T, ctable=F, landscape = T
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(c("\\begin{flushright}{", "\\scriptsize{"
                 , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
                 , "}}\\end{flushright}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
}


write_param_and_nonparam_statistics_analysis_in_latex  <- function(
  parametric_results, nonparametric_results, dvs, filename, in_title = NULL, append = F) {
  
  library(Hmisc)
  
  write(paste0("\\section{Summaries of Parametric and Nonparametric Statistics Analysis"
               , in_title, "}"), file = filename, append = append)
  write("", file = filename, append = T)
  
  ##
  result_df <- do.call(rbind, lapply(dvs, FUN = function(dv) {
    p_result <- parametric_results[[dv]]
    n_result <- nonparametric_results[[dv]]
    
    aov_df <- round(p_result$ezAov$Anova, 3)
    sch_df <- rbind(c(NA), round(n_result$sch, 3))
    
    Sig <- sapply(aov_df$`Pr(>F)`, FUN = function(x) {
      return(ifelse(x > 0.05, NA, ifelse(x >0.01, '*', '**')))
    })
    Sig[[1]] <- NA
    aov_sch_df <- cbind(aov_df, Sig)
    
    Sig <- sapply(sch_df$p.value, FUN = function(x) {
      return(ifelse(x > 0.05, NA, ifelse(x >0.01, '*', '**')))
    })
    Sig[[1]] <- NA
    aov_sch_df <- cbind(aov_sch_df, sch_df)
    aov_sch_df <- cbind(aov_sch_df, Sig)
    
    aov_sch_df <- rbind(c(NA), aov_sch_df)
    return(aov_sch_df)
  }))
  rownames(result_df) <- NULL
  
  rownames(result_df) <- c(unlist(lapply(dvs, FUN = function(dv) {
    p_result <- parametric_results[[dv]]
    return(c(dv, paste(dv, rownames(p_result$ezAov$Anova), sep = ":")))
  })))
  
  ##
  latex(result_df
        , caption = paste("Summary of two-way ANOVA and Scheirer-Ray-Hare results", in_title)
        , size = "small", longtable = T, ctable=F, landscape = T
        , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste0(c("\\begin{flushright}{", "\\scriptsize{"
                 , "Signif. codes: ", "0 ``**'' 0.01 ``*'' 0.05"
                 , "}}\\end{flushright}", "\n"), collapse = " "), file = filename, append = T)
  write("", file = filename, append = T)
  
  write_param_statistics_analysis_in_latex(
    parametric_results, dvs, filename, in_title = in_title, append = T)
  write_nonparam_statistics_analysis_in_latex(
    nonparametric_results, dvs, filename, in_title = in_title, append = T)
}

## Summary of careless responses 
write_careless_in_latex <- function(dd, filename, in_title = NULL, append = F) {
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
    write("\\section{Summary of Carless Responses}", file = filename, append = T)
  }
  
  if (!is.null(dd$get_data()) && length(dd$get_data()) > 0) {
  latex(
    as.data.frame(dd$get_data())
    , rowname = NULL
    , caption = paste("Summary of careless responses", in_title)
    , size = "small", longtable = T, ctable=F, landscape = F
    , rowlabel = "", where='!htbp', file = filename, append = T)
  }
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}


## Summary of careless responses 
write_winsorized_in_latex <- function(dd, filename, in_title = NULL, append = F) {
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
    write("\\section{Summary of Winzorized Responses}", file = filename, append = T)
  }
  
  if (!is.null(dd$get_data()) && length(dd$get_data()) > 0) {
    latex(
      as.data.frame(dd$get_data())
      , rowname = NULL
      , caption = paste("Summary of Winsorized responses", in_title)
      , size = "scriptsize", longtable = T, ctable=F, landscape = T
      , rowlabel = "", where='!htbp', file = filename, append = T)
  }
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}

## Write goodness of fit statistic fo the models
write_cfa_model_fits_in_latex <- function(fitMeasures_df, filename, in_title = NULL, append = F) {
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
    write("\\section{Goodness of Fit Statistics}", file = filename, append = T)
  }
  
  latex(
    round(fitMeasures_df, 2)
    , caption = paste("Goodness of fit statistics", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = F
    , rowlabel = "", where='!htbp', file = filename, append = T)
  write(paste(
    "\\begin{flushright}{\\tiny"
    , "df: degree of freedom; AGFI: Adjusted Goodness of Fit Index;"
    , "CFI: Comparative Fit Index; TLI: Tucker-Lewis Index;"
    , "RMSEA: Root Mean Square Error of Approximation"
    , "}\\end{flushright}\n"), file = filename, append = T)
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}

## Write rating scale model as latex
write_rsm_in_latex <- function(rsm_summaries, filename, in_title = NULL, append = F) {
  library(Hmisc)
  write("", file = filename, append = append)
  if (!append) {
    write(paste("\\documentclass[6pt]{article}"
                ,"\\usepackage{longtable}"
                ,"\\usepackage{rotating}"
                ,"\\usepackage{lscape}"
                ,"\\usepackage{ctable}"
                ,"\\begin{document}", sep = "\n")
          , file = filename, append = T)
    write("\\title{Summary of Rating Scale Model}", file = filename, append = T)
    write("\\maketitle", file = filename, append = T)
  }
  
  ##
  data_df <- rsm_summaries$model_fit
  colnames(data_df) <- sub("_", ".", colnames(data_df))
  latex(
    round(data_df, 2)
    , caption = paste("Fit statistics of RSMs", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = F
    , where='!htbp', file = filename, append = T)
  write(paste(
    "\\begin{flushright}{\\tiny"
    #, "aQ3: adjusted Q3 statistic;"
    #, "MADaQ3: Median Absolute Deviation in aQ3;"
    , "SRMSR: Standardized Root Mean Squared Residual;"
    , "WLE.rel: Weighted Likelihood Estimate reliability;"
    , "EAP.rel: Expected A Posteriori reliability;"
    , "logLik: Log Likelihood;"
    , "AIC: Akaike's Information Criterion; BIC: Bayesian Information Criterion;"
    , "}\\end{flushright}\n"), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  data_df <- rsm_summaries$estimate_params
  colnames(data_df) <- sub("_", ".", colnames(data_df))
  latex(
    round(data_df, 3)
    , caption = paste("Items parameters of the RSMs", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = T
    , where='!htbp', file = filename, append = T)
  write(paste(
    "\\begin{flushright}{\\tiny"
    , "xsi.item: Delta parameter in the item difficulty;"
    , "AXsi: Negative estimated item itercepts;"
    , "}\\end{flushright}\n"), file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  data_df <- rsm_summaries$thurstonian_thresholds
  colnames(data_df) <- sub("_", ".", colnames(data_df))
  latex(
    round(data_df, 3)
    , caption = paste("Thurstonian thresholds for the RSMs", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = T
    , where='!htbp', file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  data_df <- rsm_summaries$infit_outfit
  colnames(data_df) <- sub("_", ".", colnames(data_df))
  latex(
    round(data_df, 3)
    , caption = paste("Item infit/outfit statistic for the RSMs", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = F
    , where='!htbp', file = filename, append = T)
  write("", file = filename, append = T)
  
  ##
  data_df <- rsm_summaries$person_ability
  colnames(data_df) <- sub("_", ".", colnames(data_df))
  latex(
    round(data_df, 3)
    , caption = paste("Estimated latent traits and person model fits for the RSMs", in_title)
    , size = "scriptsize", longtable = T, ctable=F, landscape = T
    , rowlabel = "", where='!htbp', file = filename, append = T)
  write("", file = filename, append = T)
  
  if (!append) {
    write("\\end{document}", file = filename, append = T)
  }
}

