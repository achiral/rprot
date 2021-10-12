
#limma
#########################################################
#LIMMA-pipeline-proteomics/Mitochondrial_Loop/limma_helper_functions.R
#########################################################
## This is a pipeline to analyze proteiomic data for the paper Molecular wiring of a mitochondrial translational feedback loop (accepted in Molecular Cell)
## Author:Wasim Aftab
##
## Clear variables
rm(list = ls())
## Clear screen
cat('\014')
## Chdir to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
#setwd(dirname(path))
setwd("/Users/user/Dropbox/0_Work/R/limma") #作業ディレクトリ設定
cur_dir <- getwd()

## Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("limma", "qvalue")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

## Installing CRAN packages
list.of.packages <- c("dplyr", "stringr", "MASS", "plotly", "htmlwidgets")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Loading packages
library(dplyr)
library(stringr)
library(MASS)
library(plotly)
library(limma)
library(qvalue)
library(htmlwidgets)

## Load the helper functions
source("limma_helper_functions.R")

## Load the proteingroups file
myFilePath <- file.choose()
proteingroups <-
  as.data.frame(read.table(myFilePath,
                           header = TRUE,
                           sep = "\t"))

## Kill the code if proteingroups does not contain crap columns
temp <-
  select(
    proteingroups,
    matches("(Reverse|Potential.contaminant|Only.identified.by.site)")
  )
if (!nrow(temp) * ncol(temp)) {
  stop("File error, It does not contain crap...enter another file with crap")
}

## Display data to faciliate choice of treatment and control
temp <- select(proteingroups, matches("(ibaq|lfq)"))
print(names(temp))

## Remove "+" identifeid rows from proteingroups
idx <- NULL
temp <-
  select(
    proteingroups,
    matches("(Only.identified.by.site|Reverse|Potential.contaminant)")
  )
for (i in 1:ncol(temp)) {
  index <- which(unlist(!is.na(match(temp[, i], "+"))))
  idx <- union(idx, index)
}
proteingroups <- proteingroups[-idx, ] # removing indexed rows

## Remove proteins with Sequence coverage [%] < 5 [%] and Unique peptides == 1
temp <-
  select(proteingroups,
         matches("(^Sequence.coverage(.){4}$|^Unique.peptides$)"))
proteingroups <-
  proteingroups[-union(which(temp[, 1] == 1), which(temp[, 2] < 5)) , ] # removing indexed rows

## Extract Uniprot and gene symbols
Uniprot <- character(length = nrow(proteingroups))
Symbol <- character(length = nrow(proteingroups))
for (i in 1:nrow(proteingroups)) {
  temp <- as.character(proteingroups$Fasta.headers[i])
  splits <- unlist(strsplit(temp, '\\|'))
  Uniprot[i] <- splits[2]
  splits <- unlist(str_match(splits[3], "GN=(.*?) PE="))
  Symbol[i] <- splits[2]
}

## Extract data for Limma
treatment <- 
  readline('Enter treatment name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
control <- 
  readline('Enter control name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
ibaq <- 
  readinteger_binary('Enter 1 for iBAQ or 0 for LFQ= ')

if (ibaq) {
  temp <-
    select(proteingroups, matches(paste('^.*', "ibaq", '.*$', sep = '')))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol)
} else {
  temp <-
    select(proteingroups, matches(paste('^.*', "lfq", '.*$', sep = '')))
  treatment_reps <- select(temp, matches(treatment))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol)
}

## Find out Blank rows, i.e. proteins with all zeros in treatment and in control, see followig example
## iBAQ.Mrpl40_1 iBAQ.Mrpl40_2 iBAQ.Mrpl40_3 iBAQ.Kgd4_1 iBAQ.Kgd4_2 iBAQ.Kgd4_3  id Uniprot Symbol
## -------------------------------------------------------------------------------------------------
##       0             0             0           0           0           0        84  Q02888  INA17
print(names(data))
rep_treats <- 
  readinteger("Enter the number of treatment replicates=")
rep_conts <- 
  readinteger("Enter the number of control replicates=")
FC_Cutoff <- readfloat("Enter the fold change cut off=")
temp <-
  as.matrix(rowSums(apply(data[, 1:(rep_treats + rep_conts)], 2, as.numeric)))
idx <- which(temp == 0)
if (length(idx)) {
  data <- data[-idx,] # removing blank rows
}

## Find out outliers, i.e. proteins with all zeros in treatment and all/some values in control, see followig example
## Uniprot  Symbol  treat_1   treat_2   treat_3   contrl_1    contrl_2    contrl_3
## -----------------------------------------------------------------------------------
## P25554   SGF29	    0	        0	        0        2810900	     0	        0
## -----------------------------------------------------------------------------------
temp <-
  as.matrix(rowSums(apply(data[, 1:rep_treats], 2, as.numeric)))
idx <- which(temp == 0)
if (length(idx)) {
  outliers <- data[idx,]
  filename_outliers <-
    paste("Outliers_treatment_", treatment, "_", control, sep = "")
  data <- data[-idx, ] # removing indexed rows
}
## Find out outliers, i.e. proteins with all zeros in control and all/some values in treatment, see followig example
## iBAQ.Mrpl40_1 iBAQ.Mrpl40_2 iBAQ.Mrpl40_3 iBAQ.Kgd4_1 iBAQ.Kgd4_2 iBAQ.Kgd4_3  id   Uniprot   Symbol
## -----------------------------------------------------------------------------------------------------
##     662810        505600        559130        0           0           0        79   P38845    CRP1
## -----------------------------------------------------------------------------------------------------
temp <-
  as.matrix(rowSums(apply(data[, (rep_treats + 1):(rep_conts + rep_treats)], 2, as.numeric)))

idx <- which(temp == 0)
if (length(idx)) {
  outliers_control <- data[idx,]
  filename_outliers_control <-
    paste("Outliers_control_", treatment, "_", control, sep = "")
  data <- data[-idx, ] # removing indexed rows
}

## Impute missing values
data_limma <-
  log2(apply(data[c(1:(rep_treats + rep_conts))], 2, as.numeric))
data_limma[is.infinite(data_limma)] <- NA
nan_idx <- which(is.na(data_limma))

fit <- fitdistr(c(na.exclude(data_limma)), "normal")
mu <- as.double(fit$estimate[1])
sigma <- as.double(fit$estimate[2])
sigma_cutoff <- 6
new_width_cutoff <- 0.3
downshift <- 1.8
width <- sigma_cutoff * sigma
new_width <- width * new_width_cutoff
new_sigma <- new_width / sigma_cutoff
new_mean <- mu - downshift * sigma
imputed_vals_my = rnorm(length(nan_idx), new_mean, new_sigma)
data_limma[nan_idx] <- imputed_vals_my

## Read a List of Experimentally Determined Mitochondrial Proteins 
Mito_proteins <-
  as.data.frame(read.table(
    "Mitochondrial _Proteins/Experimentally_Determined_Mitochondrion_Proteins.txt",
    header = TRUE
  ))
Mito_proteins <- as.matrix(Mito_proteins$Gene_Standard_Name)

## Check and continue further only if Mitochondrial proteins are present in your list
idx <- which(data$Symbol %in% Mito_proteins)
if (length(idx)){
  data_limma_bak <- data_limma
  data_bak <- data
  data <- data[idx, ]
  data_limma <- data_limma[idx, ]
  
  ## Limma main code
  design <-
    model.matrix( ~ factor(c(rep(2, rep_treats), rep(1, rep_conts))))
  colnames(design) <- c("Intercept", "Diff")
  res.eb <- eb.fit(data_limma, design, data$Symbol)
  Sig_FC_idx <-
    union(which(res.eb$logFC < (-FC_Cutoff)), which(res.eb$logFC > FC_Cutoff))
  Sig_Pval_mod_idx <- which(res.eb$p.mod < 0.05)
  Sig_Pval_ord_idx <- which(res.eb$p.ord < 0.05)
  Sig_mod_idx <-  intersect(Sig_FC_idx, Sig_Pval_mod_idx)
  Sig_ord_idx <-  intersect(Sig_FC_idx, Sig_Pval_ord_idx)
  categ_Ord <- rep(c("Not Significant"), times = length(data$Symbol))
  categ_Mod <- categ_Ord
  categ_Mod[Sig_mod_idx] <- "Significant"
  categ_Ord[Sig_ord_idx] <- "Significant"
  dat <-
    cbind(
      res.eb,
      categ_Ord,
      categ_Mod,
      NegLogPvalMod = (-log10(res.eb$p.mod)),
      NegLogPvalOrd = (-log10(res.eb$p.ord))
    )
  
  ## Save the data files
  final_data <- cbind(data, dat)
  filename_final_data <- "final_data"
  
  ## Create plotly object and save plot as html
  filename_mod <- "limma_plot"
  filename_ord <- "ord_plot"
  result_dir <- paste("Result_", treatment, "_", control)
  display_plotly_figs(dat, FC_Cutoff, filename_mod, filename_ord, result_dir)
  
  write.table(
    final_data,
    paste(filename_final_data, '.tsv', sep = ''),
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
  )
  
  ## Write outliers in treatment
  write.table(
    outliers,
    paste(filename_outliers, '.tsv', sep = ''),
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
  )
  
  ## Write outliers in control
  write.table(
    outliers_control,
    paste(filename_outliers_control, '.tsv', sep = ''),
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
  )
  
  ## Wipe out Control Enrichments
  idx_treatment_enrich <- which(dat$logFC >= 0)
  dat <- dat[idx_treatment_enrich,]
  display_plotly_figs_half(dat, FC_Cutoff, filename_mod)
  
  ## Reset the current dir
  setwd(cur_dir)
} else {
  print('No mitochondrial proteins in your list, CODE IS TERMINATED')
}
##########Function Blocks##################################################################################
eb.fit <- function(dat, design, gene) {
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma) ^ 2
  s2.post <- fit.eb$s2.post
  t.ord <-
    fit.eb$coefficients[, 2] / fit.eb$sigma / fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <-
    data.frame(logFC,
               t.ord,
               t.mod,
               p.ord,
               p.mod,
               q.ord,
               q.mod,
               df.r,
               df.0,
               s2.0,
               s2,
               s2.post,
               gene)
  return(results.eb)
}

readinteger <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  if (is.na(n) || (n <= 0)) {
    print("Enter positive integer only")
    n <- readinteger(str)
  }
  return(n)
}

readfloat <- function(str)
{
  n <- readline(prompt = str)
  n <- as.double(n)
  if (is.na(n) || (n <= 0)) {
    print("Enter a positive number only")
    n <- readfloat(str)
  }
  return(n)
}

readinteger_binary <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  if (is.na(n) || (n < 0) || (n > 1)) {
    print("Enter positive integer only")
    n <- readinteger_binary(str)
  }
  return(n)
}

data_sanity_check <- function(temp, exprement, exp_str) {
  exprement_reps <-
    select(temp, matches(paste('^.*', exp_str, '.*$', sep = '')))
  row <- nrow(exprement_reps)
  col <- ncol(exprement_reps)
  if (!row * col || suppressWarnings(!is.na(as.integer(exp_str)))) {
    cat('Check', exprement, 'name and enter correct one\n')
    exp_str <-
      readline(
        cat(
          'Enter',
          exprement,
          'name(case insensitive) as it appeared in the iBAQ/LFQ column= '
        )
      )
    exprement_reps <- data_sanity_check(temp, 'treatment', exp_str)
  }
  return(exprement_reps)
}

display_plotly_figs <-
  function(dat,
           FC_Cutoff,
           filename_mod,
           filename_ord,
           results_dir) {
    f <- list(family = "Arial, sans-serif",
              size = 18,
              color = "#7f7f7f")
    x_axis <- list(
      constraintoward = "bottom",
      title = "log2 fold change",
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 1,
      tickwidth = 1,
      position = -1,
      rangemode = "tozero"
    )
    y_axis <- list(
      title = "-log10 p-value",
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 0.5,
      tickwidth = 1,
      position = -1,
      rangemode = "nonnegative"
    )
    
    p1 <-
      plot_ly(
        dat,
        x = ~ logFC,
        y = ~ NegLogPvalMod,
        type = "scatter",
        mode = "markers",
        color = ~ categ_Mod,
        colors = c('#0C4B8E', '#BF382A'),
        hoverinfo = 'text',
        text = ~ paste(
          "Gene:",
          dat$gene,
          "</br>",
          "</br>Fold Change:",
          logFC,
          "</br>-log 10[p-value]:",
          NegLogPvalMod
        )
      ) %>%
      layout(
        shapes = list(
          list(
            type = 'line',
            x0 = FC_Cutoff,
            x1 = FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            line = list(dash = 'dot', width = 1)
          ),
          list(
            type = 'line',
            x0 = -FC_Cutoff,
            x1 = -FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            line = list(dash = 'dot', width = 1)
          ),
          list(
            type = 'line',
            x0 = -ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            x1 = ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            y0 = -log10(0.05),
            y1 = -log10(0.05),
            line = list(dash = 'dot', width = 1)
          )
        ),
        title = paste("volcano plot of moderated p-values"),
        xaxis = x_axis,
        yaxis = y_axis,
        showlegend = TRUE
      )
    
    p2 <-
      plot_ly(
        dat,
        x = ~ logFC,
        y = ~ NegLogPvalOrd,
        type = "scatter",
        mode = "markers",
        color = ~ categ_Ord,
        colors = c('#0C4B8E', '#BF382A'),
        hoverinfo = 'text',
        text = ~ paste(
          "Gene:",
          dat$gene,
          "</br>",
          "</br>Fold Change:",
          logFC,
          "</br>-log 10[p-value]:",
          NegLogPvalOrd
        )
      ) %>%
      layout(
        shapes = list(
          list(
            type = 'line',
            x0 = FC_Cutoff,
            x1 = FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalOrd)),
            line = list(dash = 'dot', width = 1)
          ),
          list(
            type = 'line',
            x0 = -FC_Cutoff,
            x1 = -FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalOrd)),
            line = list(dash = 'dot', width = 1)
          ),
          list(
            type = 'line',
            x0 = -ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            x1 = ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            y0 = -log10(0.05),
            y1 = -log10(0.05),
            line = list(dash = 'dot', width = 1)
          )
        ),
        title = paste("volcano plot of ordinary p-values"),
        xaxis = x_axis,
        yaxis = y_axis,
        showlegend = TRUE
      )
    
    # subDir <- "Results"
    subDir <- results_dir
    mainDir <- getwd()
    
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))
    
    # ## Save as .html
    htmlwidgets::saveWidget(as_widget(p1), paste(filename_mod, '.html', sep = ''))
    htmlwidgets::saveWidget(as_widget(p2), paste(filename_ord, '.html', sep = ''))
    
    # ## Save as .pdf
    # export(p1, file = paste(filename_mod, '.pdf', sep = ''))
    # export(p2, file = paste(filename_ord, '.pdf', sep = ''))
    
    orca(p1, file = paste(filename_mod, '.pdf', sep = ''))
    orca(p2, file = paste(filename_ord, '.pdf', sep = ''))
    
  }

##display_plotly_figs_half
display_plotly_figs_half <-
  function(dat,
           FC_Cutoff,
           filename_mod) {
    f <- list(family = "Arial, sans-serif",
              size = 18,
              color = "#7f7f7f")
    x_axis <- list(
      constraintoward = "bottom",
      title = "log2 fold change",
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 1,
      tickwidth = 1,
      position = -1,
      rangemode = "tozero"
    )
    y_axis <- list(
      title = "-log10 p-value",
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 0.5,
      tickwidth = 1,
      position = -1,
      rangemode = "nonnegative"
    )
    
    p1 <-
      plot_ly(
        dat,
        x = ~ logFC,
        y = ~ NegLogPvalMod,
        type = "scatter",
        mode = "markers",
        color = ~ categ_Mod,
        colors = c('#0C4B8E', '#BF382A'),
        hoverinfo = 'text',
        text = ~ paste(
          "Gene:",
          dat$gene,
          "</br>",
          "</br>Fold Change:",
          logFC,
          "</br>-log 10[p-value]:",
          NegLogPvalMod
        )
      ) %>%
      layout(
        shapes = list(
          list(
            type = 'line',
            x0 = FC_Cutoff,
            x1 = FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            line = list(dash = 'dot', width = 1)
          ),
          list(
            type = 'line',
            x0 = 0,
            x1 = ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            y0 = -log10(0.05),
            y1 = -log10(0.05),
            line = list(dash = 'dot', width = 1)
          )
        ),
        title = paste("volcano plot of moderated p-values"),
        xaxis = x_axis,
        yaxis = y_axis,
        showlegend = TRUE
      )
    
    # ## Save as .html
    htmlwidgets::saveWidget(as_widget(p1), paste("Half_", filename_mod, '.html', sep = ''))
    
    # ## Save as .pdf
    # export(p1, file = paste("Half_", filename_mod, '.pdf', sep = ''))
    orca(p1, file = paste("Half_", filename_mod, '.pdf', sep = ''))
  }
###########################################################################################################
###########################################################################################################
###########################################################################################################