#Loading_Control
################################################################################
rm(list = ls(all = TRUE))
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", "package:utils", "package:datasets", "package:methods", "package:base")
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
  system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder", intern = FALSE, ignore.stderr = TRUE)
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
ls1 <- grep("60 kDa .eat shock protein|Actin|actin|.apdh|Histon|.inculin|Lamin.B1|Serotransferrin|Tubulin|tubulin|Vdac1|.yclophilin|.ytochrome",anno$Group)
ls2 <- grep("P68134|P68368|P05213|P60710|Q7TMM9|Q9ERD7|Q9D6F9|P99024|TUBB1|P18760|P19783|P16858|O09106|P68433|P63038|P14733|P14094|Q8VDN2|Q9WV27|P17918|P29037|Q60932|Q64727|Q00899",anno$Protein.IDs)
anno2 <- anno[ls1,]
anno3 <- anno[ls2,]
ls_rm1 <- c(grep("acting|binding|.ctinin|cuppling|deacetylase|folding|F.actin|.like|H2|H4|lysine|macro|NADPH|ractin|regulator|polymer|.related|.sociated|specific|subunit|tactin|.ytochrome.b|linking", anno2$Group),grep("P49070",anno2$`Peak Name`))
anno_LC <- anno2[-ls_rm1,]
anno_LC2 <- rbind(anno_LC, anno3)
anno_LC <- anno_LC2 %>% distinct(Protein.IDs, .keep_all = T) # anno_LCに戻す
ls_LC <- as.data.frame(anno_LC[,grep("Protein.IDs|GeneName",colnames(anno_LC))])
#output xlsx
sheets <- list("anno_LC" = anno_LC)
write_xlsx(sheets, "anno_LC.xlsx", format_headers = FALSE)
#anno_LC <- read_excel("anno_LC.xlsx", 1)
#ls_LC <- as.data.frame(anno_LC[,grep("Protein.IDs|GeneName",colnames(anno_LC))])
################################################################################
################################################################################
################################################################################
#SWATHの統計解析データ入力####←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←directory変更←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
################################################################################
################################################################################
################################################################################
#setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY") #←←←←←←directory変更←←←←←←
#setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP") #←←←←←←directory変更←←←←←←
#setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc") #←←←←←←directory変更←←←←←←
#setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC") #←←←←←←directory変更←←←←←←
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/STR") #←←←←←←directory変更←←←←←←
dat_raw <- read_excel("data.xlsx", 1)
dat_prc <- read_excel("data.xlsx", 4)
dat_raw2 <- dat_raw[,grep("SAL|PCP|Veh|CLZ",colnames(dat_raw))]
dat_prc2 <- dat_prc[,grep("SC|PC",colnames(dat_prc))]
rownames(dat_raw2) <- dat_prc$`rownames(data3)`
rownames(dat_prc2) <- dat_prc$`rownames(data3)`
colnames(dat_raw2) <- colnames(dat_prc2)
#rowname,colnameを統一させるためここまでは同時に実行
################################################################################
#raw_data解析####
fun1 <- function(dat_fun1){
  data_rm <- dat_fun1
  tdata_rm <- t(data_rm)
  tdata_rm <- cbind(as.data.frame(rownames(tdata_rm)),tdata_rm)
  colnames(tdata_rm)[1] <- "ID"
  #grouping
  group <- read_excel("SWATH.xlsx", 4) #シート4(G)入力
  PC <- factor(group$PC, levels = c("SC0", "SC10", "SC30", "PC0", "PC10", "PC30"))
  P <- factor(group$P, levels = c("S", "P"))
  C <- factor(group$C, levels = c("C0", "C10", "C30"))
  g <- cbind(PC,P,C)
  #annotation
  ganno <- group[,grep("condition|ID", colnames(group))]
  tdata_rm2 <- left_join(ganno, tdata_rm, by = "ID")
  tdata_rm3 <- tdata_rm2[,-grep("ID", colnames(tdata_rm2))]
  #statistic summary
  statv <- tdata_rm3 %>% gather(key = GeneName, value = expression, -condition) %>%
    group_by(condition, GeneName) %>%
    summarise_each(funs(N = length, mean = mean, sd = sd, se = sd/sqrt(n()), 
                        min = min, Q1 = quantile(.,0.25, na.rm=TRUE),
                        Q2 = quantile(.,0.5, na.rm=TRUE), #med = median, 
                        Q3 = quantile(., 0.75, na.rm=TRUE),
                        max = max, IQR = IQR))
  statSC0 <- statv %>% filter(condition == "SC0")
  statSC10 <- statv %>% filter(condition == "SC10")
  statSC30 <- statv %>% filter(condition == "SC30")
  statPC0 <- statv %>% filter(condition == "PC0")
  statPC10 <- statv %>% filter(condition == "PC10")
  statPC30 <- statv %>% filter(condition == "PC30")
  #colnames
  colnames(statSC0) <- str_c("SC0", colnames(statSC0), sep="_")
  colnames(statSC10) <- str_c("SC10", colnames(statSC10), sep="_")
  colnames(statSC30) <- str_c("SC30", colnames(statSC30), sep="_")
  colnames(statPC0) <- str_c("PC0", colnames(statPC0), sep="_")
  colnames(statPC10) <- str_c("PC10", colnames(statPC10), sep="_")
  colnames(statPC30) <- str_c("PC30", colnames(statPC30), sep="_")
  colnames(statSC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC30)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC30)[c(1,2)] <- c("condition","GeneName")
  #bind
  statSC0 <- statSC0[,-1]
  statSC10 <- statSC10[,-1]
  statSC30 <- statSC30[,-1]
  statPC0 <- statPC0[,-1]
  statPC10 <- statPC10[,-1]
  statPC30 <- statPC30[,-1]
  statv2 <- left_join(statSC0, statSC10, by = "GeneName")
  statv2 <- left_join(statv2, statSC30, by = "GeneName")
  statv2 <- left_join(statv2, statPC0, by = "GeneName")
  statv2 <- left_join(statv2, statPC10, by = "GeneName")
  statv2 <- left_join(statv2, statPC30, by = "GeneName")
  #############################################################
  #multcomp
  #1wANOVA function
  aof <- function(x) { 
    m <- data.frame(PC, x); 
    anova(aov(x ~ PC, m))
  }
  # apply analysis to the data and get the pvalues.
  onewayANOVA <- apply(data_rm, 1, aof)
  onewayANOVAp <- data.frame(lapply(onewayANOVA, function(x) { x["Pr(>F)"][1,] }))
  onewayANOVAp2 <- data.frame(t(onewayANOVAp))
  colnames(onewayANOVAp2) <- "p_PC" #rename
  #############################################################
  #2wANOVA function
  aof2 <- function(x) { 
    n <- data.frame(P,C, x); 
    anova(aov(x ~ P + C + P*C, n))
  }
  # apply analysis to the data and get the pvalues
  twowayANOVA <- apply(data_rm, 1, aof2)
  twowayANOVAp <- data.frame(lapply(twowayANOVA, function(x) { x["Pr(>F)"][1:3,] }))
  twowayANOVAp2 <- data.frame(t(twowayANOVAp))
  colnames(twowayANOVAp2) <- c("p_P","p_C","p_PxC") #rename
  sdata <- cbind(data_rm, onewayANOVAp2, twowayANOVAp2)
  #############################################################
  #2wANOVA BH-FDR
  #p値
  p_PC <- sdata$p_PC
  p_P  <- sdata$p_P 
  p_C <- sdata$p_C
  p_PxC <- sdata$p_PxC
  checkP <- data.frame(cbind(p_PC, p_P, p_C, p_PxC))
  rownames(checkP) <- rownames(data_rm)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "GeneName"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(data_rm)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "GeneName"
  sdata <- cbind(sdata, checkQ)
  #############################################################
  #TukeyHSD function
  THSD <- function(x) { 
    nn <- data.frame(P,C, x); 
    TukeyHSD(aov(x ~ P + C + P*C, nn))
  }
  THSDresults <- apply(data_rm, 1, THSD) 
  THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
  #############################################################
  #library(tidyverse) #ggplot2,dplyr
  #THSDp_PC <- select(THSD_PC, ends_with("p.adj")) #p値抽出
  THSDp_PC <- THSD_PC[,grep("p.adj$",colnames(THSD_PC))] #p値抽出
  #THSDd_PC <- select(THSD_PC, ends_with(".diff")) #diff値抽出
  THSDd_PC <- THSD_PC[,grep(".diff$",colnames(THSD_PC))] #diff値抽出
  #transpose
  THSDp_PC2 <- data.frame(t(THSDp_PC))
  THSDd_PC2 <- data.frame(t(THSDd_PC))
  #rename
  colnames(THSDp_PC2) <- str_c("THSDp", colnames(THSDp_PC2), sep="_")
  colnames(THSDd_PC2) <- str_c("diff", colnames(THSDd_PC2), sep="_")
  #bind
  THSDpd <- cbind(rownames(data_rm), THSDp_PC2, THSDd_PC2)
  names(THSDpd)[1] <- "GeneName"
  #############################################################
  #Annotation
  sdata2 <- cbind(rownames(sdata),sdata)
  names(sdata2)[1] <- "GeneName"
  sdata2 <- left_join(sdata2, statv2, by = "GeneName")
  sdata2 <- left_join(sdata2, THSDpd, by = "GeneName")
  sdata3 <- left_join(sdata2, anno, by = "GeneName")
  checkPr2 <- left_join(checkPr, anno, by = "GeneName")
  checkQr2 <- left_join(checkQr, anno, by = "GeneName")
  THSDpd2 <- left_join(THSDpd, anno, by = "GeneName")
  #############################################################
  #output xlsx
  sheets <- list("integ" = sdata3, "anovap" = checkPr2, "anovaq" = checkQr2, "THSDpd" = THSDpd2, "statvalue" = statv2)
  write_xlsx(sheets, "stat_raw.xlsx", format_headers = FALSE)
  #############################################################
  #DEP list
  twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
  twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
  twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
  sheets2 <- list("Pq005"=twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))],
                  "Cq005"=twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))],
                  "PxCq005"=twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
  write_xlsx(sheets2, "DEPtwANOVA_raw.xlsx", format_headers = FALSE)
}
################################################################################
#processed_data解析####
fun2 <- function(dat_fun2){
  data_rm <- dat_fun2
  tdata_rm <- t(data_rm)
  tdata_rm <- cbind(as.data.frame(rownames(tdata_rm)),tdata_rm)
  colnames(tdata_rm)[1] <- "ID"
  #grouping
  group <- read_excel("SWATH.xlsx", 4) #シート4(G)入力
  PC <- factor(group$PC, levels = c("SC0", "SC10", "SC30", "PC0", "PC10", "PC30"))
  P <- factor(group$P, levels = c("S", "P"))
  C <- factor(group$C, levels = c("C0", "C10", "C30"))
  g <- cbind(PC,P,C)
  #annotation
  ganno <- group[,grep("condition|ID", colnames(group))]
  tdata_rm2 <- left_join(ganno, tdata_rm, by = "ID")
  tdata_rm3 <- tdata_rm2[,-grep("ID", colnames(tdata_rm2))]
  #statistic summary
  statv <- tdata_rm3 %>% gather(key = GeneName, value = expression, -condition) %>%
    group_by(condition, GeneName) %>%
    summarise_each(funs(N = length, mean = mean, sd = sd, se = sd/sqrt(n()), 
                        min = min, Q1 = quantile(.,0.25, na.rm=TRUE),
                        Q2 = quantile(.,0.5, na.rm=TRUE), #med = median, 
                        Q3 = quantile(., 0.75, na.rm=TRUE),
                        max = max, IQR = IQR))
  statSC0 <- statv %>% filter(condition == "SC0")
  statSC10 <- statv %>% filter(condition == "SC10")
  statSC30 <- statv %>% filter(condition == "SC30")
  statPC0 <- statv %>% filter(condition == "PC0")
  statPC10 <- statv %>% filter(condition == "PC10")
  statPC30 <- statv %>% filter(condition == "PC30")
  #colnames
  colnames(statSC0) <- str_c("SC0", colnames(statSC0), sep="_")
  colnames(statSC10) <- str_c("SC10", colnames(statSC10), sep="_")
  colnames(statSC30) <- str_c("SC30", colnames(statSC30), sep="_")
  colnames(statPC0) <- str_c("PC0", colnames(statPC0), sep="_")
  colnames(statPC10) <- str_c("PC10", colnames(statPC10), sep="_")
  colnames(statPC30) <- str_c("PC30", colnames(statPC30), sep="_")
  colnames(statSC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC30)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC30)[c(1,2)] <- c("condition","GeneName")
  #bind
  statSC0 <- statSC0[,-1]
  statSC10 <- statSC10[,-1]
  statSC30 <- statSC30[,-1]
  statPC0 <- statPC0[,-1]
  statPC10 <- statPC10[,-1]
  statPC30 <- statPC30[,-1]
  statv2 <- left_join(statSC0, statSC10, by = "GeneName")
  statv2 <- left_join(statv2, statSC30, by = "GeneName")
  statv2 <- left_join(statv2, statPC0, by = "GeneName")
  statv2 <- left_join(statv2, statPC10, by = "GeneName")
  statv2 <- left_join(statv2, statPC30, by = "GeneName")
  #############################################################
  #multcomp
  #1wANOVA function
  aof <- function(x) { 
    m <- data.frame(PC, x); 
    anova(aov(x ~ PC, m))
  }
  # apply analysis to the data and get the pvalues.
  onewayANOVA <- apply(data_rm, 1, aof)
  onewayANOVAp <- data.frame(lapply(onewayANOVA, function(x) { x["Pr(>F)"][1,] }))
  onewayANOVAp2 <- data.frame(t(onewayANOVAp))
  colnames(onewayANOVAp2) <- "p_PC" #rename
  #############################################################
  #2wANOVA function
  aof2 <- function(x) { 
    n <- data.frame(P,C, x); 
    anova(aov(x ~ P + C + P*C, n))
  }
  # apply analysis to the data and get the pvalues
  twowayANOVA <- apply(data_rm, 1, aof2)
  twowayANOVAp <- data.frame(lapply(twowayANOVA, function(x) { x["Pr(>F)"][1:3,] }))
  twowayANOVAp2 <- data.frame(t(twowayANOVAp))
  colnames(twowayANOVAp2) <- c("p_P","p_C","p_PxC") #rename
  sdata <- cbind(data_rm, onewayANOVAp2, twowayANOVAp2)
  #############################################################
  #2wANOVA BH-FDR
  #p値
  p_PC <- sdata$p_PC
  p_P  <- sdata$p_P 
  p_C <- sdata$p_C
  p_PxC <- sdata$p_PxC
  checkP <- data.frame(cbind(p_PC, p_P, p_C, p_PxC))
  rownames(checkP) <- rownames(data_rm)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "GeneName"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(data_rm)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "GeneName"
  sdata <- cbind(sdata, checkQ)
  #############################################################
  #TukeyHSD function
  THSD <- function(x) { 
    nn <- data.frame(P,C, x); 
    TukeyHSD(aov(x ~ P + C + P*C, nn))
  }
  THSDresults <- apply(data_rm, 1, THSD) 
  THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
  #############################################################
  #library(tidyverse) #ggplot2,dplyr
  #THSDp_PC <- select(THSD_PC, ends_with("p.adj")) #p値抽出
  THSDp_PC <- THSD_PC[,grep("p.adj$",colnames(THSD_PC))] #p値抽出
  #THSDd_PC <- select(THSD_PC, ends_with(".diff")) #diff値抽出
  THSDd_PC <- THSD_PC[,grep(".diff$",colnames(THSD_PC))] #diff値抽出
  #transpose
  THSDp_PC2 <- data.frame(t(THSDp_PC))
  THSDd_PC2 <- data.frame(t(THSDd_PC))
  #rename
  colnames(THSDp_PC2) <- str_c("THSDp", colnames(THSDp_PC2), sep="_")
  colnames(THSDd_PC2) <- str_c("diff", colnames(THSDd_PC2), sep="_")
  #bind
  THSDpd <- cbind(rownames(data_rm), THSDp_PC2, THSDd_PC2)
  names(THSDpd)[1] <- "GeneName"
  #############################################################
  #Annotation
  sdata2 <- cbind(rownames(sdata),sdata)
  names(sdata2)[1] <- "GeneName"
  sdata2 <- left_join(sdata2, statv2, by = "GeneName")
  sdata2 <- left_join(sdata2, THSDpd, by = "GeneName")
  sdata3 <- left_join(sdata2, anno, by = "GeneName")
  checkPr2 <- left_join(checkPr, anno, by = "GeneName")
  checkQr2 <- left_join(checkQr, anno, by = "GeneName")
  THSDpd2 <- left_join(THSDpd, anno, by = "GeneName")
  #############################################################
  #output xlsx
  sheets <- list("integ" = sdata3, "anovap" = checkPr2, "anovaq" = checkQr2, "THSDpd" = THSDpd2, "statvalue" = statv2)
  write_xlsx(sheets, "stat_prc.xlsx", format_headers = FALSE)
  #############################################################
  #DEP list
  twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
  twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
  twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
  sheets2 <- list("Pq005"=twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))],
                  "Cq005"=twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))],
                  "PxCq005"=twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
  write_xlsx(sheets2, "DEPtwANOVA_prc.xlsx", format_headers = FALSE)
}
fun1(dat_raw2)
fun2(dat_prc2)
################################################################################
################################################################################
################################################################################
#データ加工####
################################################################################
################################################################################
################################################################################
dat_raw <- read_excel("stat_raw.xlsx", 1)
dat_prc <- read_excel("stat_prc.xlsx", 1)
NAval <- dat_prc %>% filter(is.na(Protein.IDs)) 
dupli <- dat_prc %>% group_by(GeneName) %>% filter(n()>1)
dat_prc <- dat_prc %>% distinct(GeneName, .keep_all = T) 
LC_raw <- left_join(ls_LC, dat_raw[,-1], by = "Protein.IDs")
LC_raw_mean <- LC_raw[,grep("^Protein|GeneName|mean$",colnames(LC_raw))] #mean値抽出
LC_raw_mean <- LC_raw_mean %>% filter(!is.na(SC0_mean)) #NA値除去
LC_prc <- left_join(ls_LC, dat_prc[,-1], by = "Protein.IDs")
LC_prc_mean <- LC_prc[,grep("^Protein|GeneName|mean$",colnames(LC_raw))] #mean値抽出
LC_prc_mean <- LC_prc_mean %>% filter(!is.na(SC0_mean)) #NA値除去
geomean <- function(dat_geomean){ #幾何平均関数定義
  exp(mean(log(abs(dat_geomean))))
} #幾何平均算出
mean_table <- function(dat_mean){
  #dat_mean<- LC_prc_mean #おまじない
  x2 <- dat_mean %>% mutate(all_samples = apply(dat_mean[,-(1:2)], 1, mean)) #全サンプル平均値算出
  #ACTB,G3P(GAPDH)
  ACTBxG3P <- x2[grep("P60710|P16858",dat_mean$Protein.IDs),]
  ACTBxG3P_m <- t(as.data.frame(c("P60710xP16858m","ACTBxG3Pm",apply(ACTBxG3P[,-(1:2)],2,mean)))) #ACTB,G3P(GAPDH)の算術平均値
  ACTBxG3P_gm <- t(as.data.frame(c("P60710xP16858gm","ACTBxG3Pgm",apply(ACTBxG3P[,-(1:2)],2,geomean)))) #ACTB,G3P(GAPDH)の幾何平均値
  colnames(ACTBxG3P_m) <- colnames(x2)
  rownames(ACTBxG3P_m) <- NULL
  colnames(ACTBxG3P_gm) <- colnames(x2)
  rownames(ACTBxG3P_gm) <- NULL
  #TBA4A,G3P(GAPDH)
  TBA4AxG3P <- x2[grep("P68368|P16858",dat_mean$Protein.IDs),]
  TBA4AxG3P_m <- t(as.data.frame(c("P68368xP16858m","TBA4AxG3Pm",apply(TBA4AxG3P[,-(1:2)],2,mean)))) #TBA4A,G3P(GAPDH)の算術平均値
  TBA4AxG3P_gm <- t(as.data.frame(c("P68368xP16858gm","TBA4AxG3Pgm",apply(TBA4AxG3P[,-(1:2)],2,geomean)))) #TBA4A,G3P(GAPDH)の幾何平均値
  colnames(TBA4AxG3P_m) <- colnames(x2)
  rownames(TBA4AxG3P_m) <- NULL
  colnames(TBA4AxG3P_gm) <- colnames(x2)
  rownames(TBA4AxG3P_gm) <- NULL
  #TBB3,G3P(GAPDH)
  TBB3xG3P <- x2[grep("Q9ERD7|P16858",dat_mean$Protein.IDs),]
  TBB3xG3P_m <- t(as.data.frame(c("Q9ERD7xP16858m","TBB3xG3Pm",apply(TBB3xG3P[,-(1:2)],2,mean)))) #TBB3,G3P(GAPDH)の算術平均値
  TBB3xG3P_gm <- t(as.data.frame(c("Q9ERD7xP16858gm","TBB3xG3Pgm",apply(TBB3xG3P[,-(1:2)],2,geomean)))) #TBB3,G3P(GAPDH)の幾何平均値
  colnames(TBB3xG3P_m) <- colnames(x2)
  rownames(TBB3xG3P_m) <- NULL
  colnames(TBB3xG3P_gm) <- colnames(x2)
  rownames(TBB3xG3P_gm) <- NULL 
  #TBB4A,G3P(GAPDH)
  TBB4AxG3P <- x2[grep("Q9D6F9|P16858",dat_mean$Protein.IDs),]
  TBB4AxG3P_m <- t(as.data.frame(c("Q9D6F9xP16858m","TBB4AxG3Pm",apply(TBB4AxG3P[,-(1:2)],2,mean)))) #TBB4A,G3P(GAPDH)の算術平均値
  TBB4AxG3P_gm <- t(as.data.frame(c("Q9D6F9xP16858gm","TBB4AxG3Pgm",apply(TBB4AxG3P[,-(1:2)],2,geomean)))) #TBB4A,G3P(GAPDH)の幾何平均値
  colnames(TBB4AxG3P_m) <- colnames(x2)
  rownames(TBB4AxG3P_m) <- NULL
  colnames(TBB4AxG3P_gm) <- colnames(x2)
  rownames(TBB4AxG3P_gm) <- NULL
  x3 <- rbind(x2,ACTBxG3P_m,ACTBxG3P_gm,TBA4AxG3P_m,TBA4AxG3P_gm,
              TBB3xG3P_m,TBB3xG3P_gm,TBB4AxG3P_m,TBB4AxG3P_gm)
  #return(invisible(x3))
} #LC_prc_meanフォーマットの全サンプル算術平均とACTBxG3Pの算術・幾何平均を追加
opt <- mean_table(LC_prc_mean) #output
LC_prc_mean <- opt #LC_prc_meanに戻す
dat_prc_value <- dat_prc[,grep("GeneName|_1$|_2$|_3$|_4$|_5$|_6$|_7$|_8$|_9$|_10$|_11$|_12$",colnames(dat_prc))] #数値データ抽出
#rownames(dat_prc_value) <- dat_prc$GeneName
#Loading_Control_Mean_List保存
write_xlsx(LC_prc_mean, "LC_prc_mean.xlsx", format_headers = FALSE)
################################################################################
#計算結果を保存した複数のオブジェクトをcsv出力
#https://ja.stackoverflow.com/questions/34975/
div <- NULL
div_all <- NULL
x <-NULL
x_all <- NULL
#list_name <- paste(1:nrow(LC_prc_mean), ".csv", sep ="")
for(i in 1:nrow(LC_prc_mean)){
  x <- dat_prc_value[,-1]/as.numeric(LC_prc_mean[i,grep("^SC0_mean$",colnames(LC_prc_mean))])
  x <- cbind(dat_prc_value[,1],x)
  x_all <- dat_prc_value[,-1]/as.numeric(LC_prc_mean[i,grep("^all_samples$",colnames(LC_prc_mean))])
  x_all <- cbind(dat_prc_value[,1],x_all)
  div[[i]] <- x
  div_all[[i]] <- x_all
  #write.csv(div[[i]], list_name[i]) #個別にcsv出力
}
write_xlsx(div, "div.xlsx", format_headers = FALSE) #まとめてexcel出力
write_xlsx(div_all, "div_all.xlsx", format_headers = FALSE) #まとめてexcel出力
################################################################################
################################################################################
################################################################################
#統計値####
################################################################################
################################################################################
################################################################################
#コントロール(Saline/Vehicle)サンプルの平均で補正####
################################################################################
ls_stat <- NULL
ls_DEPtwANOVA_Pq005 <-NULL
ls_DEPtwANOVA_Cq005 <-NULL
ls_DEPtwANOVA_PxCq005 <-NULL
for(j in 1:nrow(LC_prc_mean)){
#シート番号を入力する
dat_prc <- read_excel("div.xlsx", j)
dat_prc2 <- dat_prc[,-1]
rownames(dat_prc2) <- dat_prc$GeneName
################################################################################
#processed_data解析
data_rm <- dat_prc2
tdata_rm <- t(data_rm)
tdata_rm <- cbind(as.data.frame(rownames(tdata_rm)),tdata_rm)
colnames(tdata_rm)[1] <- "ID"
#grouping
group <- read_excel("SWATH.xlsx", 4) #シート4(G)入力
PC <- factor(group$PC, levels = c("SC0", "SC10", "SC30", "PC0", "PC10", "PC30"))
P <- factor(group$P, levels = c("S", "P"))
C <- factor(group$C, levels = c("C0", "C10", "C30"))
g <- cbind(PC,P,C)
#annotation
ganno <- group[,grep("condition|ID", colnames(group))]
tdata_rm2 <- left_join(ganno, tdata_rm, by = "ID")
tdata_rm3 <- tdata_rm2[,-grep("ID", colnames(tdata_rm2))]
#statistic summary
statv <- tdata_rm3 %>% gather(key = GeneName, value = expression, -condition) %>%
  group_by(condition, GeneName) %>%
  summarise_each(funs(N = length, mean = mean, sd = sd, se = sd/sqrt(n()), 
                      min = min, Q1 = quantile(.,0.25, na.rm=TRUE),
                      Q2 = quantile(.,0.5, na.rm=TRUE), #med = median, 
                      Q3 = quantile(., 0.75, na.rm=TRUE),
                      max = max, IQR = IQR))
statSC0 <- statv %>% filter(condition == "SC0")
statSC10 <- statv %>% filter(condition == "SC10")
statSC30 <- statv %>% filter(condition == "SC30")
statPC0 <- statv %>% filter(condition == "PC0")
statPC10 <- statv %>% filter(condition == "PC10")
statPC30 <- statv %>% filter(condition == "PC30")
#colnames
colnames(statSC0) <- str_c("SC0", colnames(statSC0), sep="_")
colnames(statSC10) <- str_c("SC10", colnames(statSC10), sep="_")
colnames(statSC30) <- str_c("SC30", colnames(statSC30), sep="_")
colnames(statPC0) <- str_c("PC0", colnames(statPC0), sep="_")
colnames(statPC10) <- str_c("PC10", colnames(statPC10), sep="_")
colnames(statPC30) <- str_c("PC30", colnames(statPC30), sep="_")
colnames(statSC0)[c(1,2)] <- c("condition","GeneName")
colnames(statSC10)[c(1,2)] <- c("condition","GeneName")
colnames(statSC30)[c(1,2)] <- c("condition","GeneName")
colnames(statPC0)[c(1,2)] <- c("condition","GeneName")
colnames(statPC10)[c(1,2)] <- c("condition","GeneName")
colnames(statPC30)[c(1,2)] <- c("condition","GeneName")
#bind
statSC0 <- statSC0[,-1]
statSC10 <- statSC10[,-1]
statSC30 <- statSC30[,-1]
statPC0 <- statPC0[,-1]
statPC10 <- statPC10[,-1]
statPC30 <- statPC30[,-1]
statv2 <- left_join(statSC0, statSC10, by = "GeneName")
statv2 <- left_join(statv2, statSC30, by = "GeneName")
statv2 <- left_join(statv2, statPC0, by = "GeneName")
statv2 <- left_join(statv2, statPC10, by = "GeneName")
statv2 <- left_join(statv2, statPC30, by = "GeneName")
#############################################################
#multcomp
#1wANOVA function
aof <- function(x) { 
  m <- data.frame(PC, x); 
  anova(aov(x ~ PC, m))
}
# apply analysis to the data and get the pvalues.
onewayANOVA <- apply(data_rm, 1, aof)
onewayANOVAp <- data.frame(lapply(onewayANOVA, function(x) { x["Pr(>F)"][1,] }))
onewayANOVAp2 <- data.frame(t(onewayANOVAp))
colnames(onewayANOVAp2) <- "p_PC" #rename
#############################################################
#2wANOVA function
aof2 <- function(x) { 
  n <- data.frame(P,C, x); 
  anova(aov(x ~ P + C + P*C, n))
}
# apply analysis to the data and get the pvalues
twowayANOVA <- apply(data_rm, 1, aof2)
twowayANOVAp <- data.frame(lapply(twowayANOVA, function(x) { x["Pr(>F)"][1:3,] }))
twowayANOVAp2 <- data.frame(t(twowayANOVAp))
colnames(twowayANOVAp2) <- c("p_P","p_C","p_PxC") #rename
sdata <- cbind(data_rm, onewayANOVAp2, twowayANOVAp2)
#############################################################
#BH-FDR
#p値
p_PC <- sdata$p_PC
p_P  <- sdata$p_P 
p_C <- sdata$p_C
p_PxC <- sdata$p_PxC
checkP <- data.frame(cbind(p_PC, p_P, p_C, p_PxC))
rownames(checkP) <- rownames(data_rm)
checkPr <- cbind(rownames(checkP),checkP)
names(checkPr)[1] <- "GeneName"
#q値
q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
q_P <- data.frame(p.adjust(p_P, method = "BH"))
q_C <- data.frame(p.adjust(p_C, method = "BH"))
q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
rownames(checkQ) <- rownames(data_rm)
checkQr <- cbind(rownames(checkQ),checkQ)
names(checkQr)[1] <- "GeneName"
sdata <- cbind(sdata, checkQ)
#############################################################
#TukeyHSD function
THSD <- function(x) { 
  nn <- data.frame(P,C, x); 
  TukeyHSD(aov(x ~ P + C + P*C, nn))
}
THSDresults <- apply(data_rm, 1, THSD) 
THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
#############################################################
#library(tidyverse) #ggplot2,dplyr
#THSDp_PC <- select(THSD_PC, ends_with("p.adj")) #p値抽出
THSDp_PC <- THSD_PC[,grep("p.adj$",colnames(THSD_PC))] #p値抽出
#THSDd_PC <- select(THSD_PC, ends_with(".diff")) #diff値抽出
THSDd_PC <- THSD_PC[,grep(".diff$",colnames(THSD_PC))] #diff値抽出
#transpose
THSDp_PC2 <- data.frame(t(THSDp_PC))
THSDd_PC2 <- data.frame(t(THSDd_PC))
#rename
colnames(THSDp_PC2) <- str_c("THSDp", colnames(THSDp_PC2), sep="_")
colnames(THSDd_PC2) <- str_c("diff", colnames(THSDd_PC2), sep="_")
#bind
THSDpd <- cbind(rownames(data_rm), THSDp_PC2, THSDd_PC2)
names(THSDpd)[1] <- "GeneName"
#############################################################
#Annotation
sdata2 <- cbind(rownames(sdata),sdata)
names(sdata2)[1] <- "GeneName"
sdata2 <- left_join(sdata2, statv2, by = "GeneName")
sdata2 <- left_join(sdata2, THSDpd, by = "GeneName")
sdata3 <- left_join(sdata2, anno, by = "GeneName")
checkPr2 <- left_join(checkPr, anno, by = "GeneName")
checkQr2 <- left_join(checkQr, anno, by = "GeneName")
THSDpd2 <- left_join(THSDpd, anno, by = "GeneName")
#############################################################
#output xlsx
ls_stat[[j]] <- as.data.frame(sdata3)
#write_xlsx(ls_stat, "stat_prc.xlsx", format_headers = FALSE)
#############################################################
#DEP list
twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
ls_DEPtwANOVA_Pq005[[j]] <- as.data.frame(twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))])
ls_DEPtwANOVA_Cq005[[j]] <- as.data.frame(twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))])
ls_DEPtwANOVA_PxCq005[[j]] <- as.data.frame(twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
}
write_xlsx(ls_stat, "stat_prc_div.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_Pq005, "DEPtwANOVA_Pq005_prc_div.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_Cq005, "DEPtwANOVA_Cq005_prc_div.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_PxCq005, "DEPtwANOVA_PxCq005_prc_div.xlsx", format_headers = FALSE)
################################################################################
#全サンプルの平均で補正####
################################################################################
ls_stat_all <- NULL
ls_DEPtwANOVA_Pq005_all <-NULL
ls_DEPtwANOVA_Cq005_all <-NULL
ls_DEPtwANOVA_PxCq005_all <-NULL
for(k in 1:nrow(LC_prc_mean)){
  #シート番号を入力する
  dat_prc <- read_excel("div_all.xlsx", k)
  dat_prc2 <- dat_prc[,-1]
  rownames(dat_prc2) <- dat_prc$GeneName
  ################################################################################
  #processed_data解析
  data_rm <- dat_prc2
  tdata_rm <- t(data_rm)
  tdata_rm <- cbind(as.data.frame(rownames(tdata_rm)),tdata_rm)
  colnames(tdata_rm)[1] <- "ID"
  #grouping
  group <- read_excel("SWATH.xlsx", 4) #シート4(G)入力
  PC <- factor(group$PC, levels = c("SC0", "SC10", "SC30", "PC0", "PC10", "PC30"))
  P <- factor(group$P, levels = c("S", "P"))
  C <- factor(group$C, levels = c("C0", "C10", "C30"))
  g <- cbind(PC,P,C)
  #annotation
  ganno <- group[,grep("condition|ID", colnames(group))]
  tdata_rm2 <- left_join(ganno, tdata_rm, by = "ID")
  tdata_rm3 <- tdata_rm2[,-grep("ID", colnames(tdata_rm2))]
  #statistic summary
  statv <- tdata_rm3 %>% gather(key = GeneName, value = expression, -condition) %>%
    group_by(condition, GeneName) %>%
    summarise_each(funs(N = length, mean = mean, sd = sd, se = sd/sqrt(n()), 
                        min = min, Q1 = quantile(.,0.25, na.rm=TRUE),
                        Q2 = quantile(.,0.5, na.rm=TRUE), #med = median, 
                        Q3 = quantile(., 0.75, na.rm=TRUE),
                        max = max, IQR = IQR))
  statSC0 <- statv %>% filter(condition == "SC0")
  statSC10 <- statv %>% filter(condition == "SC10")
  statSC30 <- statv %>% filter(condition == "SC30")
  statPC0 <- statv %>% filter(condition == "PC0")
  statPC10 <- statv %>% filter(condition == "PC10")
  statPC30 <- statv %>% filter(condition == "PC30")
  #colnames
  colnames(statSC0) <- str_c("SC0", colnames(statSC0), sep="_")
  colnames(statSC10) <- str_c("SC10", colnames(statSC10), sep="_")
  colnames(statSC30) <- str_c("SC30", colnames(statSC30), sep="_")
  colnames(statPC0) <- str_c("PC0", colnames(statPC0), sep="_")
  colnames(statPC10) <- str_c("PC10", colnames(statPC10), sep="_")
  colnames(statPC30) <- str_c("PC30", colnames(statPC30), sep="_")
  colnames(statSC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statSC30)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC0)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC10)[c(1,2)] <- c("condition","GeneName")
  colnames(statPC30)[c(1,2)] <- c("condition","GeneName")
  #bind
  statSC0 <- statSC0[,-1]
  statSC10 <- statSC10[,-1]
  statSC30 <- statSC30[,-1]
  statPC0 <- statPC0[,-1]
  statPC10 <- statPC10[,-1]
  statPC30 <- statPC30[,-1]
  statv2 <- left_join(statSC0, statSC10, by = "GeneName")
  statv2 <- left_join(statv2, statSC30, by = "GeneName")
  statv2 <- left_join(statv2, statPC0, by = "GeneName")
  statv2 <- left_join(statv2, statPC10, by = "GeneName")
  statv2 <- left_join(statv2, statPC30, by = "GeneName")
  #############################################################
  #multcomp
  #1wANOVA function
  aof <- function(x) { 
    m <- data.frame(PC, x); 
    anova(aov(x ~ PC, m))
  }
  # apply analysis to the data and get the pvalues.
  onewayANOVA <- apply(data_rm, 1, aof)
  onewayANOVAp <- data.frame(lapply(onewayANOVA, function(x) { x["Pr(>F)"][1,] }))
  onewayANOVAp2 <- data.frame(t(onewayANOVAp))
  colnames(onewayANOVAp2) <- "p_PC" #rename
  #############################################################
  #2wANOVA function
  aof2 <- function(x) { 
    n <- data.frame(P,C, x); 
    anova(aov(x ~ P + C + P*C, n))
  }
  # apply analysis to the data and get the pvalues
  twowayANOVA <- apply(data_rm, 1, aof2)
  twowayANOVAp <- data.frame(lapply(twowayANOVA, function(x) { x["Pr(>F)"][1:3,] }))
  twowayANOVAp2 <- data.frame(t(twowayANOVAp))
  colnames(twowayANOVAp2) <- c("p_P","p_C","p_PxC") #rename
  sdata <- cbind(data_rm, onewayANOVAp2, twowayANOVAp2)
  #############################################################
  #BH-FDR
  #p値
  p_PC <- sdata$p_PC
  p_P  <- sdata$p_P 
  p_C <- sdata$p_C
  p_PxC <- sdata$p_PxC
  checkP <- data.frame(cbind(p_PC, p_P, p_C, p_PxC))
  rownames(checkP) <- rownames(data_rm)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "GeneName"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(data_rm)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "GeneName"
  sdata <- cbind(sdata, checkQ)
  #############################################################
  #TukeyHSD function
  THSD <- function(x) { 
    nn <- data.frame(P,C, x); 
    TukeyHSD(aov(x ~ P + C + P*C, nn))
  }
  THSDresults <- apply(data_rm, 1, THSD) 
  THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
  #############################################################
  #library(tidyverse) #ggplot2,dplyr
  #THSDp_PC <- select(THSD_PC, ends_with("p.adj")) #p値抽出
  THSDp_PC <- THSD_PC[,grep("p.adj$",colnames(THSD_PC))] #p値抽出
  #THSDd_PC <- select(THSD_PC, ends_with(".diff")) #diff値抽出
  THSDd_PC <- THSD_PC[,grep(".diff$",colnames(THSD_PC))] #diff値抽出
  #transpose
  THSDp_PC2 <- data.frame(t(THSDp_PC))
  THSDd_PC2 <- data.frame(t(THSDd_PC))
  #rename
  colnames(THSDp_PC2) <- str_c("THSDp", colnames(THSDp_PC2), sep="_")
  colnames(THSDd_PC2) <- str_c("diff", colnames(THSDd_PC2), sep="_")
  #bind
  THSDpd <- cbind(rownames(data_rm), THSDp_PC2, THSDd_PC2)
  names(THSDpd)[1] <- "GeneName"
  #############################################################
  #Annotation
  sdata2 <- cbind(rownames(sdata),sdata)
  names(sdata2)[1] <- "GeneName"
  sdata2 <- left_join(sdata2, statv2, by = "GeneName")
  sdata2 <- left_join(sdata2, THSDpd, by = "GeneName")
  sdata3 <- left_join(sdata2, anno, by = "GeneName")
  checkPr2 <- left_join(checkPr, anno, by = "GeneName")
  checkQr2 <- left_join(checkQr, anno, by = "GeneName")
  THSDpd2 <- left_join(THSDpd, anno, by = "GeneName")
  #############################################################
  #output xlsx
  ls_stat_all[[k]] <- as.data.frame(sdata3)
  #write_xlsx(ls_stat, "stat_prc.xlsx", format_headers = FALSE)
  #############################################################
  #DEP list
  twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
  twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
  twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
  ls_DEPtwANOVA_Pq005_all[[k]] <- as.data.frame(twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))])
  ls_DEPtwANOVA_Cq005_all[[k]] <- as.data.frame(twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))])
  ls_DEPtwANOVA_PxCq005_all[[k]] <- as.data.frame(twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
}
write_xlsx(ls_stat_all, "stat_prc_div_all.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_Pq005_all, "DEPtwANOVA_Pq005_prc_div_all.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_Cq005_all, "DEPtwANOVA_Cq005_prc_div_all.xlsx", format_headers = FALSE)
write_xlsx(ls_DEPtwANOVA_PxCq005_all, "DEPtwANOVA_PxCq005_prc_div_all.xlsx", format_headers = FALSE)
################################################################################
################################################################################
################################################################################
tictoc::toc() #処理時間計測終了
