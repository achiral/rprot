

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstatsQCgui")
BiocManager::install("devtools", dependencies = TRUE)
deps <- c("devtools")
devtools::install_github("eralpdogu/MSstatsQCgui")
