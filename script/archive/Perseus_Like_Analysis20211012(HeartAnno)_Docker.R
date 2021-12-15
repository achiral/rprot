#Perseus_Like_Analysis
#Make_Annotation_List
#log2-Impute(MNAR)-Subtract(Median):like a Perseus
################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db", "mouse4302.db","GO.db",
#                       "PANTHER.db", "biomaRt"))
################################################################################
rm(list = ls(all = TRUE))
source("/home/rstudio/rproject/script/archive/functions.R")
detach_all()
source("/home/rstudio/rproject/script/startup.R")
################################################################################
# dirname = dir.choose()
# filename = file.choose()
################################################################################

## $$$$$ set data table directory $$$$$ #####
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/AMY")
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/Heart")
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/HIP")
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/NAc")
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/PFC")
setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/STR")


getwd()
setwd("/home/rstudio/rproject")
# annotation table
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/AMY")
dat_a <- read_excel("SWATH.xlsx", 2)
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/HIP")
dat_h <- read_excel("SWATH.xlsx", 2)
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/NAc")
dat_n <- read_excel("SWATH.xlsx", 2)
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")
dat_p <- read_excel("SWATH.xlsx", 2)
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/STR")
dat_s <- read_excel("SWATH.xlsx", 2)
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/Heart")
dat_ht <- read_excel("SWATH.xlsx", 2)
################################################################################
#DIR <- "/home/rstudio/rproject/data/Perseus_Like_Analysis/Other3"  # set path
#dir.create(DIR, showWarnings = T, recursive = T)                   # make directory
# unlink(DIR, recursive = T)                                       # remove directory
# file.create("temp1.txt", showWarnings = T)                       # make file
# file.copy("temp1.txt","temp2.txt", overwrite = T)                # copy file o/w
# file.rename("temp2.txt", "temp3.txt")                            # rename file
# file.remove("temp1.txt","temp2.txt","temp3.txt")                 # remove file
################################################################################
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/Other")
getwd()
dir() 

t(colnames(dat_a))
num <- grep("(Peak Name|Group)",colnames(dat_a))
x <- rbind(dat_a[,num],dat_h[,num],dat_ht[,num],dat_n[,num],dat_p[,num],dat_s[,num])
#split,extract
split_pn <- data.frame(str_split(x$`Peak Name`, pattern = "\\|", simplify = TRUE))
colnames(split_pn) <- c("sp", "Protein.IDs", "GeneName") #列名変更
Protein.IDs <- data.frame(str_sub(split_pn$`Protein.IDs`, start = 1, end = 6)) #`Protein.IDs`列の1-6文字目(Protein.IDs)抽
Gene.names <- data.frame(str_sub(split_pn$`GeneName`, start = 1, end = -7)) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
Species <- data.frame(str_sub(split_pn$`GeneName`, start = -5, end = -1)) #`GeneName`列の-5〜-1文字目(Species)抽出
split_pn2 <- cbind(Protein.IDs, Gene.names, Species)
colnames(split_pn2) <- c("Protein.IDs", "GeneName", "Species") #列名変更
split_gr <- data.frame(str_split(x$`Group`, pattern = ".OS=|.GN=|.PE=|.SV=", simplify = TRUE))
colnames(split_gr) <- c("Description", "OS", "GN", "PE", "SV") #列名変更
xx <- cbind(x, split_pn2, split_gr)
# Remove duplication
xxx <- xx %>% distinct(Protein.IDs,.keep_all=TRUE)
# Search Duplication
xxx$Protein.IDs %>% duplicated() %>% any()
# Duplication table
xxx %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
################################################################################
################################################################################
################################################################################
# output annotation table
# write_xlsx(xxx, "anno.xlsx", format_headers = FALSE)
################################################################################
################################################################################
################################################################################
# SWATHのAnnotation情報にEntrezIDなど追加
anno <- xxx
#生物種レベルのアノテーション（OrgDb）
id <- anno$`Protein.IDs`
GN <- anno$GN
#GeneName <- anno$GeneName
res_id <- select(org.Mm.eg.db, keys = as.character(id), keytype = "UNIPROT", 
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN <- select(org.Mm.eg.db, keys = as.character(GN), keytype = "SYMBOL",
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN <- res_GN[,c(6,2,3,4,5,1)] #arrange columns
#rbind
res_id_GN <- rbind(res_id, res_GN)
#remove duplicates
ex_id <- res_id_GN %>% distinct(UNIPROT, .keep_all = T)
ex_GN <- res_id_GN %>% distinct(SYMBOL, .keep_all = T)
ex_res_id_GN <- rbind(ex_id, ex_GN) %>% filter(!is.na(ENTREZID)) %>% filter(!is.na(UNIPROT)) %>% distinct(UNIPROT, .keep_all = T)
ex_res_id_GN_Other <- rbind(ex_id, ex_GN) %>% filter(!is.na(ENTREZID)) %>% filter(is.na(UNIPROT)) %>% distinct(SYMBOL, .keep_all = T)
#left_join
anno_id <- left_join(anno, ex_res_id_GN, by = c("Protein.IDs" = "UNIPROT"))
anno_GN <- left_join(anno, ex_res_id_GN, by = c("GN" = "SYMBOL"))
anno_id_GN <- rbind(anno_id[1:14], anno_GN[-11]) %>% filter(!is.na(ENTREZID)) %>% distinct(Protein.IDs, .keep_all = T)

anno_id_Other <- left_join(anno, ex_res_id_GN_Other, by = c("Protein.IDs" = "UNIPROT"))
anno_GN_Other <- left_join(anno, ex_res_id_GN_Other, by = c("GN" = "SYMBOL"))
anno_id_GN_Other <- rbind(anno_id_Other[1:14], anno_GN_Other[-11]) %>% filter(!is.na(ENTREZID)) %>% distinct(Protein.IDs, .keep_all = T)

anno2 <- left_join(anno, rbind(anno_id_GN[,c(3,11:14)], anno_id_GN_Other[,c(3,11:14)]), by = "Protein.IDs")
#not NA value
anno2_notNA <- anno2 %>% filter(!is.na(ENTREZID))
#NA value
anno2_NA <- anno2 %>% filter(is.na(ENTREZID))
anno2_NA_Mm <- anno2_NA %>% filter(Species == "MOUSE")
anno2_NA_Other <- anno2_NA %>% filter(Species != "MOUSE")

#remove libraries
detach_all()
library(org.Mm.eg.db)
#entrezID searched from internet
ent <- c("18563", "234695", "100503185", "14467", "14070")
res_ent <- select(org.Mm.eg.db, keys = ent, keytype = "ENTREZID",
                  columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
#remove duplicates
library(tidyverse)
res_ent <- res_ent %>% filter(!is.na(ENTREZID)) %>% distinct(ENTREZID, .keep_all = T)
res_ent[1,]
res_ent <- res_ent[,c(2,1,3,4,5)]
#cbind
anno2_NA_Mm <- cbind(anno2_NA_Mm[,1:10], res_ent[,1:4])

t(colnames(anno2))
t(colnames(anno2_notNA))
t(colnames(anno2_NA_Mm))
t(colnames(anno2_NA_Other))

#rbind
anno3 <- rbind(anno2_notNA, anno2_NA_Mm, anno2_NA_Other)
anno3_NA <- anno3%>% filter(is.na(Protein.IDs)) #Check NA value
#Original order
t(colnames(anno3))
anno_final <- left_join(anno,anno3[,c(3,11:14)],by = "Protein.IDs")
################################################################################
################################################################################
################################################################################
#output xlsx
library(openxlsx) #入出力(write.xlsx)
smp <- list("anno_new"=anno_final,"anno"=anno)
write.xlsx(smp, file = "anno.xlsx", overwrite = T)
################################################################################


################################################################################
################################################################################
#Perseus_Like_Analysis
################################################################################
################################################################################
rm(list = ls(all = TRUE))
source("/home/rstudio/rproject/script/archive/functions.R")        # load functions
source("/home/rstudio/rproject/script/archive/functions_DEPpkg.R") # load DEPpkg functions
detach_all()
source("/home/rstudio/rproject/script/startup.R")                  # load packages
################################################################################
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/Other") # set directory
# setwd("/Users/user/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/Other")
anno <- read_excel("anno.xlsx", 1)                               # annotation
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/AMY")   # set directory
PLA2()                                                           # analysis
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/HIP")   # set directory
PLA2()                                                           # analysis
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/NAc")   # set directory
PLA2()                                                           # analysis
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")   # set directory
PLA2()                                                           # analysis
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/STR")   # set directory
PLA2()                                                           # analysis
setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/Heart") # set directory
# setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/Heart")
PLA2()                                                           # analysis
################################################################################