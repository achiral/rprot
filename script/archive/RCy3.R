
#RCy3:Cytoscape
#https://f1000research.com/articles/8-1774/v3
##############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RCy3")
library(RCy3)
##############################################################
cytoscapePing()
#help(package=RCy3)
#browseVignettes("RCy3")

