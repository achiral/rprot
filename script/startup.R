# Startup
################################################################################
# set working directory
setwd("/home/rstudio/rproject")
################################################################################
# do analysis session
## Session > Restart R
## File > Open Project > /home/rstudio/project/new_project.Rproj
## File > New Project
## File > New File > R Script
################################################################################
# load package
# pacman::p_load(MSstats, SummarizedExperiment, tidyverse)
# require(biomaRt)            # biomaRt data base
require(DESeq2)
require(dplyr)                # select function
require(GO.db)                # gene ontology
require(MSstats)
require(multcomp)             # multiple comparison
require(openxlsx)             # in/output xlsx (write.xlsx)
require(org.Hs.eg.db)         # human
require(org.Mm.eg.db)         # mouse
require(PANTHER.db)           # PANTHER data base
require(readxl)               # input xlsx (read_excel)
require(SummarizedExperiment) # DEP
require(tidyverse)            # ggplot2,dplyr
require(writexl)              # output
require(xlsx)                 # input xlsx
################################################################################




