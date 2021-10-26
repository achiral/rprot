# Startup
################################################################################
# set working directory
# setwd("/home/rstudio/rproject")
################################################################################
# do analysis session
## Session > Restart R
## File > Open Project > /home/rstudio/project/new_project.Rproj
## File > New Project
## File > New File > R Script
################################################################################
# load package
# pacman::p_load(MSstats, SummarizedExperiment, tidyverse)
require(AnnotationHub)
# require(biomaRt)            # biomaRt data base
require(clusterProfiler)      # enrichment analysis
require(DESeq2)
require(DOSE)
require(dplyr)                # select function
require(plyr)
require(europepmc)

# require(fontregisterer)       # (https://ill-identified.hatenablog.com/entry/2020/10/03/200618#ux3069ux306eosux3067ux3082ux3060ux3044ux305fux3044ux540cux3058ux3088ux3046ux306bux3067ux304dux308bux65b9ux6cd5)
# require(ggplot2)
require(ggnewscale)           # enrichment map plot
require(ggupset)
require(GO.db)                # gene ontology
require(gridExtra)            # svg
require(Hmisc)
require(MeSHDbi)
require(meshes)
require(MeSH.Hsa.eg.db) # BioC 2.14-3.13
require(MeSH.Mmu.eg.db) # BioC 2.14-3.13
require(msigdbr)
require(MSstats)
require(multcomp)             # multiple comparison
require(openxlsx)             # in/output xlsx (write.xlsx)
require(org.Hs.eg.db)         # human
require(org.Mm.eg.db)         # mouse
require(PANTHER.db)           # PANTHER data base
require(pathview)
require(ReactomePA)
require(readxl)               # input xlsx (read_excel)
# require(rvg)                # error
# require(stringr)
require(SummarizedExperiment) # DEP
# require(systemfonts)          # error
# require(svglite)            # error
require(tidyverse)            # ggplot2,dplyr
require(writexl)              # output
require(xlsx)                 # input xlsx
################################################################################
if(Sys.info()["sysname"] == "Windows"){
  if(as.integer(str_extract(Sys.info()["release"], "^[0-9]+")) >=8){
    family_sans <- "Yu Gothic"
    family_serif <- "Yu Mincho"
  } else {
    family_sans <- "MS Gothic"
    family_serif <- "MS Mincho"
  }
} else if(Sys.info()["sysname"] == "Linux") {
  family_sans <- "Noto Sans CJK JP"
  family_serif <- "Noto Serif CJK JP"
} else if(Sys.info()["sysname"] == "Darwin"){
  family_serif <- "Hiragino Mincho ProN"
  family_sans <- "Hiragino Sans"
} else {
  # インストールすればとりあえず動く
  family_sans <- "Noto Sans CJK JP"
  family_serif <- "Noto Serif CJK JP"
}



