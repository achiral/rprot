#Loading_Control_List
################################################################################
rm(list = ls(all = TRUE))
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", 
                 "package:utils", "package:datasets", "package:methods", "package:base")
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  pkg.list <- setdiff(pkg.list, basic.pkg)
  lapply(pkg.list, detach, character.only = TRUE)
}
detach_all()
library(tidyverse) #ggplot2,dplyr
library(readxl) #入力(read_excel)
library(xlsx) #入力
library(openxlsx) #入出力(write.xlsx)
library(writexl) #出力
library(multcomp)
library(tictoc)
tictoc::tic() #処理時間計測開始
################################################################################
#フルパスの確認
dir.choose <- function() {
  system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}
#dirname = dir.choose()
#dirname
#filename = file.choose()
################################################################################
################################################################################
################################################################################
#Annotation table作成####
################################################################################
################################################################################
################################################################################
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
anno <- read_excel("anno.xlsx", 1)
anno_LC <-read_excel("anno_LC.xlsx", 1)

setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY")
LC_prc_mean_a <- read_excel("LC_prc_mean.xlsx", 1)
colnames(LC_prc_mean_a) <- str_c("a", colnames(LC_prc_mean_a), sep="_") #rename
LC_prc_mean_a <- LC_prc_mean_a %>% rename(Protein.IDs = 1)
LC_prc_THSDp_a <- read_excel("LC_prc_THSDp.xlsx", 1)
colnames(LC_prc_THSDp_a) <- str_c("a", colnames(LC_prc_THSDp_a), sep="_") #rename
LC_prc_THSDp_a <- LC_prc_THSDp_a %>% rename(Protein.IDs = 1)

setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP")
LC_prc_mean_h <- read_excel("LC_prc_mean.xlsx", 1)
colnames(LC_prc_mean_h) <- str_c("h", colnames(LC_prc_mean_h), sep="_") #rename
LC_prc_mean_h <- LC_prc_mean_h %>% rename(Protein.IDs = 1)
LC_prc_THSDp_h <- read_excel("LC_prc_THSDp.xlsx", 1)
colnames(LC_prc_THSDp_h) <- str_c("h", colnames(LC_prc_THSDp_h), sep="_") #rename
LC_prc_THSDp_h <- LC_prc_THSDp_h %>% rename(Protein.IDs = 1)

setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc")
LC_prc_mean_n <- read_excel("LC_prc_mean.xlsx", 1)
colnames(LC_prc_mean_n) <- str_c("n", colnames(LC_prc_mean_n), sep="_") #rename
LC_prc_mean_n <- LC_prc_mean_n %>% rename(Protein.IDs = 1)
LC_prc_THSDp_n <- read_excel("LC_prc_THSDp.xlsx", 1)
colnames(LC_prc_THSDp_n) <- str_c("n", colnames(LC_prc_THSDp_n), sep="_") #rename
LC_prc_THSDp_n <- LC_prc_THSDp_n %>% rename(Protein.IDs = 1)

setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC")
LC_prc_mean_p <- read_excel("LC_prc_mean.xlsx", 1)
colnames(LC_prc_mean_p) <- str_c("p", colnames(LC_prc_mean_p), sep="_") #rename
LC_prc_mean_p <- LC_prc_mean_p %>% rename(Protein.IDs = 1)
LC_prc_THSDp_p <- read_excel("LC_prc_THSDp.xlsx", 1)
colnames(LC_prc_THSDp_p) <- str_c("p", colnames(LC_prc_THSDp_p), sep="_") #rename
LC_prc_THSDp_p <- LC_prc_THSDp_p %>% rename(Protein.IDs = 1)

setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/STR")
LC_prc_mean_s <- read_excel("LC_prc_mean.xlsx", 1)
colnames(LC_prc_mean_s) <- str_c("s", colnames(LC_prc_mean_s), sep="_") #rename
LC_prc_mean_s <- LC_prc_mean_s %>% rename(Protein.IDs = 1)
LC_prc_THSDp_s <- read_excel("LC_prc_THSDp.xlsx", 1)
colnames(LC_prc_THSDp_s) <- str_c("s", colnames(LC_prc_THSDp_s), sep="_") #rename
LC_prc_THSDp_s <- LC_prc_THSDp_s %>% rename(Protein.IDs = 1)

#共通
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}
LC_all <- as.data.frame(intersect_all(LC_prc_mean_a$Protein.IDs, LC_prc_mean_h$Protein.IDs, LC_prc_mean_n$Protein.IDs, LC_prc_mean_p$Protein.IDs, LC_prc_mean_s$Protein.IDs))
colnames(LC_all) <- "Protein.IDs"
LC_all <- left_join(LC_all, anno, by = "Protein.IDs")
LC_ahnps <- left_join(anno_LC, LC_prc_mean_a, by = "Protein.IDs")
LC_ahnps <- left_join(LC_ahnps, LC_prc_mean_h, by = "Protein.IDs")
LC_ahnps <- left_join(LC_ahnps, LC_prc_mean_n, by = "Protein.IDs")
LC_ahnps <- left_join(LC_ahnps, LC_prc_mean_p, by = "Protein.IDs")
LC_ahnps <- left_join(LC_ahnps, LC_prc_mean_s, by = "Protein.IDs")

LC_ahnps_p <- left_join(anno_LC, LC_prc_THSDp_a, by = "Protein.IDs")
LC_ahnps_p <- left_join(LC_ahnps_p, LC_prc_THSDp_h, by = "Protein.IDs")
LC_ahnps_p <- left_join(LC_ahnps_p, LC_prc_THSDp_n, by = "Protein.IDs")
LC_ahnps_p <- left_join(LC_ahnps_p, LC_prc_THSDp_p, by = "Protein.IDs")
LC_ahnps_p <- left_join(LC_ahnps_p, LC_prc_THSDp_s, by = "Protein.IDs")


## THSDp >= 0.05を999に修正
map_res <- LC_ahnps_p %>% select_if(is.numeric) %>% #型が数値の列のみ選択
  map_dbl(~ sum(. >= 0.05))  #mapで各列に0.05以上の要素をカウント
names(map_res)[map_res > 0]  #0.05以上の要素が存在する列名
res <- LC_ahnps_p %>% mutate_if(names(.) %in% names(map_res)[map_res > 0], ~ if_else(. >= 0.05, 999, .)) #999に置換


setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
sheets <- list("anno_LC" = anno_LC, "LC_all" = LC_all, "LC_ahnps" = LC_ahnps, "LC_ahnps_p" = LC_ahnps_p, "LC_ahnps_p_res" = res)
write_xlsx(sheets,"anno_LC.xlsx", format_headers = F)
