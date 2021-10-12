# new_script
################################################################################
setwd("/home/rstudio/project")
getwd()
################################################################################
# install_packages("pkg_name")
# require("pkg_name")
# pacman::p_load(pkg_name, pkg_name,...,pkg_name)
# library(pkg_name)
# ex)  library(palmerpenguins)
#      head(penguins)
library(DEP)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
################################################################################