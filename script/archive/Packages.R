# R packages
#############################################################
# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
BiocManager::install()                              # update
BiocManager::install("airway")
BiocManager::install("aLFQ")
BiocManager::install("BaylorEdPsych")               # NA
BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")             # enrichmant analysis
BiocManager::install("ComplexHeatmap")              # heatmap
BiocManager::install("DEP")                       # NA, proteomics
BiocManager::install("DESeq")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
biocManager::install("DOSE")
BiocManager::install("genefilter")                  # heatmap
BiocManager::install("gmm")                         # NA, DEP
BiocManager::install("GO.db")
BiocManager::install("gplots")                      # heatmap
BiocManager::install("imsbInfer")                   # NA
BiocManager::install("limma")                       # DEP
BiocManager::install("loadTransitonsMSExperiment")  # NA
BiocManager::install("MeSH.Hsa.eg.db")              # BioC 2.14-3.13
BiocManager::install("MeSH.Mmu.eg.db")              # BioC 2.14-3.13
BiocManager::install("MeSHDbi")                     # BioC 3.14 (Nov. 2021, with R-4.2.0)
BiocManager::install("meshes")
BiocManager::install("mouse4302.db")
BiocManager::install("msigdbr")
BiocManager::install("MSstats")
BiocManager::install("mzR")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("PANTHER.db")
BiocManager::install("pathview")
BiocManager::install("PECA")
BiocManager::install("RCyjs")   # cyREST, Cytoscape
BiocManager::install("ReactomePA")
BiocManager::install("RforProteomics")
BiocManager::install("readinteger_binary")
BiocManager::install("SummarizedExperiment")      # DEP
BiocManager::install("SWATH2stats")

devtools::install_github("GuangchuangYu/DOSE")    # enrichMap

remotes::install_github("Gedevan-Aleksizde/fontregisterer", repos = NULL, type = "source")


install.packages("agricolae")                       
install.packages("BH")
install.packages("BiocManager")                   # done
install.packages("cowplot")
install.packages("data.table")                      # MSstat
install.packages("devtools")
install.packages("doSNOW")                          # MSstat

install.packages("ddply")
install.packages("dplyr")                         # tidyverse, MSstat

install.packages("europepmc")

install.packages("foreach")                         # MSstat
install.packages("ggcorrplot") 
install.packages("ggnewscale")                      # enrichment analysis
install.packages("ggrepel")                         # MSstat
install.packages("ggplot2")                         # MSstat
install.packages("ggThemeAssist")

install.packages("ggupset")

install.packages("githubinstall")
install.packages("gplots")                          # MSstat
install.packages("gridExtra")                       # svg
install.packages("Hmisc")                           # llist()
install.packages("imputeMissings")
install.packages("lattice")                         # done, limma
install.packages("lme4")                            # MSstat
install.packages("mgcv")                            # done, limma
install.packages("mice")
install.packages("minpack.lm")                      # MSstat
install.packages("missForest")
install.packages("mlbench")
install.packages("multcomp")                      # done, multiple comparison
install.packages("openxlsx")                      # done, exlx input
install.packages("pacman")
install.packages("palmerpenguins")                # demo
install.packages("randomForest")                    # MSstat
install.packages("ranger")
install.packages("remotes")
install.packages("readxl")                        # done, exlx input, read_excel
install.packages("reshape")                         # MSstat
install.packages("reshape2")                        # MSstat
install.packages("rJava")                         # done
install.packages("rsvg")                            # svg
# install.packages('RSVGTipsDevice')                # removed from the CRAN repository
install.packages("rvg")                             # svg（https://www.karada-good.net/analyticsr/r-382）
install.packages("S4Vectors")                     # done, DEP/SummarizedExperiment
install.packages("sandwich")                      # done, gmm
install.packages("scales")                          # muted()
install.packages("sets")                            # set operation
install.packages("sgof")                            # bh(), Multiple Hypothesis Testing
install.packages("snow")                            # MSstat
install.packages("stringr")                         # MSstat
install.packages("survival")                        # MSstat
install.packages("svglite")                       # svg(https://ill-identified.hatenablog.com/entry/2020/10/03/200618#svg-ux3067ux306eux4fddux5b58)
install.packages("systemfonts")                   # rvg()
install.packages("tablaxlsx")                     # NA, xlsx table output
install.packages("tidyr")                         # MSstat
install.packages("tidyverse")                     # done, 
install.packages("VIM")
install.packages("writexl")                       # done, writexl output
install.packages("xlsx")                          # done, xlsx output
# install.packages("xlsx2")                         # NA, xlsx output
install.packages("XLConnect")                       # xlsx in/output
#############################################################
# Load BioConductor Packages
library(airway)
library(aLFQ)
# library(BaylorEdPsych)                            # NA
# library(BiocManager)                              # done
# library(biomaRt)
library(ComplexHeatmap)                             # heatmap
# library(DEP)                                      # NA, proteomics
library(DESeq)
library(DESeq2)
library(EnhancedVolcano)
library(genefilter)                                 # heatmap
# library(gmm)                                      # NA, DEP
# library(GO.db)
library(gplots)                                     # heatmap
library(imsbInfer)                                  # NA
library(loadTransitonsMSExperiment)                 # NA
# library(mouse4302.db)
# library(MSstats)                                  # NA
library(mzR)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(PANTHER.db)
library(PECA)
library(RCyjs)                                      # cyREST, Cytoscape
library(readinteger_binary)
# library(SummarizedExperiment)                     # DEP
library(SummarizedExperiment)
library(SWATH2stats)

#############################################################
# Load BaseR Packages
library(agricolae)                                  
library(BH)                                         # FDR
library(cowplot)
# library(data.table)                               # DEP, MSstats
library(devtools)
library(doSNOW)                                     # MSstats
# library(dplyr)                                    # tidyverse, MSstats
require(fontregisterer)       # (https://ill-identified.hatenablog.com/entry/2020/10/03/200618#ux3069ux306eosux3067ux3082ux3060ux3044ux305fux3044ux540cux3058ux3088ux3046ux306bux3067ux304dux308bux65b9ux6cd5)
# library(foreach)                                  # DEP, MSstats
library(ggcorrplot)                                 
library(ggrepel)                                    # MSstats
library(ggplot2)                                    # MSstats
library(ggThemeAssist)
# library(gplots)                                   # MSstats, heatmap
# library(gridExtra)                                # svg, MSstats
library(imputeMissings)
# library(lattice)                                  # done, DEP, limma
# library(lme4)                                     # MSstats
# library(magrittr)                                 # tidyverse?
# library(mgcv)                                     # done, DEP?, limma
library(mice)
library(minpack.lm)                                 # MSstats
library(missForest)
library(mlbench)
# library(multcomp)                                 # done, multiple comparison
# library(openxlsx)                                 # done, exlx output(write.xlsx)
# library(palmerpenguins)                           # demo
# library(RColorBrewer)                             # color, tidyverse?
# library(readxl)                                   # base?, exlx input(read_excel)
library(reshape)                                    # MSstats
library(reshape2)                                   # MSstats
library(randomForest)                               # MSstats
# library(rJava)                                    # done
library(rsvg)                                       # svg
# library(RSVGTipsDevice)                             # error
library(rvg)                                        # svg（https://www.karada-good.net/analyticsr/r-382）
# library(S4Vectors)                                # DEP/SummarizedExperiment
# library(sandwich)                                 # DEP?, gmm?
# library(scales)                                   # muted(), tidyverse?
library(sets)                                       # set operation
library(sgof)                                       # bh(), Multiple Hypothesis Testing
# library(snow)                                     # MSstats
# library(stringr)                                  # MSstats
# library(survival)                                 # MSstats
# library(tablaxlsx)                                # NA, xlsx table output
# library(tidyr)                                    # MSstats
# library(tidyverse)                                # ggplot2, dplyr
library(VIM)
library(writexl)                                    # xlsx output
library(xlsx)                                       # xlsx output
# library(xlsx2)                                    # NA, xlsx output
library(XLConnect)                                  # NA(JAVA8-11), xlsx in/output
#############################################################
installed.packages()
sessionInfo() 
#############################################################
