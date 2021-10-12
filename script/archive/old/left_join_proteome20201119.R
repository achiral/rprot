#Annotation table作成
#解析データの統合
##############################################################
library(readxl) #エクセル入力(read_excel)
library(tidyverse)
library(writexl) #xlsx出力
rm(list = ls(all.names = TRUE))
##############################################################
#Annotation table作成
##############################################################
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/AMY")
dat_a <- read_excel("SWATH.xlsx", 2)
dat_as <- read_excel("stat.xlsx", 1)
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/HIP")
dat_h <- read_excel("SWATH.xlsx", 2)
dat_hs <- read_excel("stat.xlsx", 1)
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/NAc")
dat_n <- read_excel("SWATH.xlsx", 2)
dat_ns <- read_excel("stat.xlsx", 1)
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/PFC")
dat_p <- read_excel("SWATH.xlsx", 2)
dat_ps <- read_excel("stat.xlsx", 1)
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/STR")
dat_s <- read_excel("SWATH.xlsx", 2)
dat_ss <- read_excel("stat.xlsx", 1)
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/Other")
getwd()
dir() 
t(colnames(dat_a))
x <- rbind(dat_a[,c(2,5)],dat_h[,c(2,5)],dat_n[,c(2,5)],dat_p[,c(2,5)],dat_s[,c(2,5)])
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
#Remove duplication
xxx <- xx %>% distinct(Protein.IDs,.keep_all=TRUE)
#Search Duplication
xxx$Protein.IDs %>% duplicated() %>% any()
#Duplication table
xxx %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
##############################################################
#Annotation table出力
##############################################################
write_xlsx(xxx, "anno.xlsx", format_headers = FALSE)
##############################################################
##############################################################
##############################################################
#解析データの統合
##############################################################
library(readxl) #エクセル入力(read_excel)
library(tidyverse)
library(writexl) #xlsx出力
rm(list = ls(all.names = TRUE))
##############################################################
#AMY
#setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/AMY")
#vdat_a <- read_excel("AMY.xlsx", 1)
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY") 
#dat_as <- read_excel("stat.xlsx", 1)
#HIP
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/HIP")
dat_h <- read_excel("HIP.xlsx", 1)
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP") 
dat_hs <- read_excel("stat.xlsx", 1)
#NAc
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/NAc")
dat_n <- read_excel("NAc.xlsx", 1)
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc") 
dat_ns <- read_excel("stat.xlsx", 1)
#PFC
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/PFC")
dat_p <- read_excel("PFC.xlsx", 1)
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC") 
dat_ps <- read_excel("stat.xlsx", 1)
#STR
#setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/STR")
#dat_s <- read_excel("STR.xlsx", 1)
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR") 
#dat_ss <- read_excel("stat.xlsx", 1)
#統合left_join
#dat_a2 <- left_join(dat_a, dat_as, by = c("Group" = "Description"))
dat_h2 <- left_join(dat_h, dat_hs, by = c("Group" = "Description"))
dat_n2 <- left_join(dat_n, dat_ns, by = c("Group" = "Description"))
dat_p2 <- left_join(dat_p, dat_ps, by = c("Group" = "Description"))
#dat_s2 <- left_join(dat_s, dat_ss, by = c("Group" = "Description"))
#統合right_join
#dat_a3 <- right_join(dat_a, dat_as, by = c("Group" = "Description"))
#dat_h3 <- right_join(dat_h, dat_hs, by = c("Group" = "Description"))
#dat_n3 <- right_join(dat_n, dat_ns, by = c("Group" = "Description"))
#dat_p3 <- right_join(dat_p, dat_ps, by = c("Group" = "Description"))
#dat_s3 <- right_join(dat_s, dat_ss, by = c("Group" = "Description"))
#統合left_join
#dat_a4 <- left_join(dat_as, dat_a, by = c("Description" ="Group"))
dat_h4 <- left_join(dat_hs, dat_h, by = c("Description" ="Group"))
dat_n4 <- left_join(dat_ns, dat_n, by = c("Description" ="Group"))
dat_p4 <- left_join(dat_ps, dat_p, by = c("Description" ="Group"))
#dat_s4 <- left_join(dat_ss, dat_s, by = c("Description" ="Group"))
#############################################################
#xlsx出力
library(writexl) #xlsx出力
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R/Other")
#sheets <- list("AMY" = dat_a2, "HIP" = dat_h2, "NAc" = dat_n2, "PFC" = dat_p2, "STR" = dat_s2)
sheets2 <- list("HIP" = dat_h2, "NAc" = dat_n2, "PFC" = dat_p2)
#sheets3 <- list("HIP" = dat_h3, "NAc" = dat_n3, "PFC" = dat_p3)
sheets4 <- list("HIP" = dat_h4, "NAc" = dat_n4, "PFC" = dat_p4)
write_xlsx(sheets2, "stat2.xlsx", format_headers = FALSE)
#write_xlsx(sheets3, "stat3.xlsx", format_headers = FALSE)
write_xlsx(sheets4, "stat4.xlsx", format_headers = FALSE)
#txt出力
#write.table (sdata2, file = "integ.txt", sep = "\t") #保存
#############################################################