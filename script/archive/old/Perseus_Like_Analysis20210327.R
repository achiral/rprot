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
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", 
                 "package:utils", "package:datasets", "package:methods", "package:base")
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  pkg.list <- setdiff(pkg.list, basic.pkg)
  lapply(pkg.list, detach, character.only = TRUE)
}
detach_all()
library(DEP)
library(tidyverse) #ggplot2,dplyr
library(dplyr)
library(readxl) #入力(read_excel)
library(xlsx) #入力
library(openxlsx) #入出力(write.xlsx)
library(writexl) #出力
library(multcomp)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(GO.db)
#library(PANTHER.db)
#library(biomaRt)
################################################################################
#フルパスの確認(https://qiita.com/h398qy988q5/items/7e0052b29ec876407f5d)
dir.choose <- function() {
  system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}
#dirname = dir.choose()
#filename = file.choose()
################################################################################
#Annotation table作成
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY")
dat_a <- read_excel("SWATH.xlsx", 2)
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP")
dat_h <- read_excel("SWATH.xlsx", 2)
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc")
dat_n <- read_excel("SWATH.xlsx", 2)
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC")
dat_p <- read_excel("SWATH.xlsx", 2)
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/STR")
dat_s <- read_excel("SWATH.xlsx", 2)
setwd("~/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
getwd()
dir() 
t(colnames(dat_a))
num <- grep("(Peak Name|Group)",colnames(dat_a))
x <- rbind(dat_a[,num],dat_h[,num],dat_n[,num],dat_p[,num],dat_s[,num])
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
#Annotation table出力
write_xlsx(xxx, "anno.xlsx", format_headers = FALSE)
################################################################################
#SWATHのAnnotation情報にEntrezIDなど追加
anno <- xxx
#生物種レベルのアノテーション（OrgDb）
id <- anno$`Protein.IDs`
GN <- anno$GN
#GeneName <- anno$GeneName
res_id <- select(org.Mm.eg.db, keys = id, keytype = "UNIPROT",
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN <- select(org.Mm.eg.db, keys = GN, keytype = "SYMBOL",
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN <- res_GN[,c(6,2,3,4,5,1)]
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
ent <- c("18563", "234695", "14467", "14070")
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

#output xlsx
library(openxlsx) #入出力(write.xlsx)
smp <- list("anno_new"=anno_final,"anno"=anno)
write.xlsx(smp, "anno.xlsx")
################################################################################
################################################################################
#Perseus_Like_Analysis
################################################################################
################################################################################
rm(list = ls(all.names = TRUE))
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", 
                 "package:utils", "package:datasets", "package:methods", "package:base")
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  pkg.list <- setdiff(pkg.list, basic.pkg)
  lapply(pkg.list, detach, character.only = TRUE)
}
detach_all()
library(DEP)
library(tidyverse) #ggplot2,dplyr
library(dplyr)
library(readxl) #入力(read_excel)
library(xlsx) #入力
library(openxlsx) #入出力(write.xlsx)
library(writexl) #出力
library(multcomp)
################################################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
anno <- read_excel("anno.xlsx", 1) #シート1入力
################################################################################
#統計解析関数(引数2)
#Log2transform,Imputation(MNAR),Subtraction(Median),1wANOVA,2wANOVA,THSD
fun2 <- function(x,y){
  data <- x
  ExpDesign <- y
  #split
  split <- str_split(data$`Peak Name`, pattern = "\\|", simplify = TRUE)
  colnames(split) <- c("sp", "Protein.IDs", "GeneName") #列名変更
  class(split)
  x <- data.frame(split)
  #extract
  Protein.IDs <- str_sub(x$`Protein.IDs`, start = 1, end = 6) #`Peak Name`列の1-6文字目(Protein.IDs)抽出
  Gene.names <- str_sub(x$`GeneName`, start = 1, end = -7) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
  Species <- str_sub(x$`GeneName`, start = -5, end = -1) #`GeneName`列の-5〜-1文字目(Species)抽出
  #bind
  data <- cbind(data, Protein.IDs, Gene.names, Species) #data, Protein.IDs, Gene.names, Speciesを列ベクトル単位で結合
  #Search Duplication
  data$Protein.IDs %>% duplicated() %>% any()
  data$Gene.names %>% duplicated() %>% any()
  data$Species %>% duplicated() %>% any()
  #Duplication table
  data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  data %>% group_by(Species) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  #Unique Uniprot ID
  data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
  data_unique$Protein.IDs %>% duplicated() %>% any() # Are there any duplicated names?
  #SummarizedExperiment
  Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
  experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
  ###############################################################################
  #Log2-transform
  data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数, #Log2-transformation
  data1 <- data.frame(data_se@assays@data) #log2
  #Impute:left-shifted Gaussian distribution (for MNAR)
  data_imp_man <- impute(data_se, fun = "man", shift = 1.8, scale = 0.3) #Perseus,imputation
  data2 <- data.frame(data_imp_man@assays@data) #Subtract前log2imp
  #Subtract(Median):Perseus
  standardize <- function(z) {
    colmed <- apply(z, 2, median) #Median of Each Sample's Protein Expression level
    colmad <- apply(z, 2, mad)  # median absolute deviation
    rv <- sweep(z, 2, colmed,"-")  #subtracting median expression
    #rv <- sweep(rv, 2, colmad, "/")  # dividing by median absolute deviation
    return(rv)
  }
  data3 <- data2 #Subtract前log2impをコピー
  Sample_columns <- grep("(SC|PC)", colnames(data3)) # get Sample column numbers
  data3[Sample_columns] <- standardize(data3[Sample_columns]) #Subtract(Median),log2impsub
  #############################################################
  dat1 <- cbind(rownames(data1),data1) #log2
  dat2 <- cbind(rownames(data2),data2) #log2imp
  dat3 <- cbind(rownames(data3),data3) #log2impsub
  #integration
  dat <- cbind(data$Gene.names,data) #行名追加
  dat4 <- left_join(dat, dat1, by = c("Gene.names" = "rownames(data1)")) #raw+log2
  dat4 <- left_join(dat4, dat2, by = c("Gene.names" = "rownames(data2)")) #raw+log2+log2imp
  dat4 <- left_join(dat4, dat3, by = c("Gene.names" = "rownames(data3)")) #raw+log2+log2imp+log2impsub
  #output xlsx
  smp <- list("raw"=dat,"log2"=dat1,"log2imp"=dat2,"log2impsub"=dat3,"integ"=dat4,"anno"=anno) #リスト作成,rawdata,log2,imputation,subtract,integration
  write.xlsx(smp, "data.xlsx")
  #############################################################
  #statistic summary
  data_rm  <- data3
  data_rm[,1:2] <- NULL #列削除
  #transpose
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
  rownames(checkP) <- rownames(data3)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "GeneName"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(data3)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "GeneName"
  sdata <- cbind(sdata, checkQ)
  #############################################################
  #TukeyHSD function
  #diff群間の平均値の差(例)B-Aが-127.3であればデータBの平均がデータAの平均より-127.3大きい
  #lwr,upr=下方信頼限界,情報信頼限界:信頼区間の下限値 (lower) と上限値 (upper)
  #0を含まない場合 (例)B-A は含まず D-A は含む=2群間差は0ではないので有意差あり
  #p.adj < 0.05=2群間に有意差あり(信頼区間内に0を含まない)
  #############################################################
  THSD <- function(x) { 
    nn <- data.frame(P,C, x); 
    TukeyHSD(aov(x ~ P + C + P*C, nn))
  }
  THSDresults <- apply(data_rm, 1, THSD) 
  THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
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
  THSDpd <- cbind(rownames(data3), THSDp_PC2, THSDd_PC2)
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
  sheets <- list("integ" = sdata3, "anovap" = checkPr2, 
                 "anovaq" = checkQr2, "THSDpd" = THSDpd2, 
                 "statvalue" = statv2) #assume sheet1-4 are data frames
  write_xlsx(sheets, "stat.xlsx", format_headers = FALSE)
  #############################################################
  #DEP list
  twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
  twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
  twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
  sheets2 <- list("Pq005"=twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))],
                  "Cq005"=twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))],
                  "PxCq005"=twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
  write_xlsx(sheets2, "DEPtwANOVA.xlsx", format_headers = FALSE)
}
################################################################################
#Amygdala
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY")
data <- read_excel("SWATH.xlsx", 2) #swath data
ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
fun2(data, ExpDesign)
#Hippocampus
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP")
data <- read_excel("SWATH.xlsx", 2) #swath data
ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
fun2(data, ExpDesign)
#Nucleus accumbens
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc")
data <- read_excel("SWATH.xlsx", 2) #swath data
ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
fun2(data, ExpDesign)
#Prefrontal cortex
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC")
data <- read_excel("SWATH.xlsx", 2) #swath data
ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
fun2(data, ExpDesign)
#Striatum
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR")
data <- read_excel("SWATH.xlsx", 2) #swath data
ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
fun2(data, ExpDesign)
################################################################################
################################################################################
#引数1つのパターン
################################################################################
################################################################################
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
#anno <- read_excel("anno.xlsx", 1) #シート1入力
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY")
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP")
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc")
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC")
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR")
#data <- read_excel("SWATH.xlsx", 2) #swath data
#ExpDesign <- read_excel("SWATH.xlsx", 3) #DEP.packcage SE file
#fun1 <- function(data){
  #split
  split <- str_split(data$`Peak Name`, pattern = "\\|", simplify = TRUE)
  colnames(split) <- c("sp", "Protein.IDs", "GeneName") #列名変更
  class(split)
  x <- data.frame(split)
  #extract
  Protein.IDs <- str_sub(x$`Protein.IDs`, start = 1, end = 6) #`Peak Name`列の1-6文字目(Protein.IDs)抽出
  Gene.names <- str_sub(x$`GeneName`, start = 1, end = -7) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
  Species <- str_sub(x$`GeneName`, start = -5, end = -1) #`GeneName`列の-5〜-1文字目(Species)抽出
  #bind
  data <- cbind(data, Protein.IDs, Gene.names, Species) #data, Protein.IDs, Gene.names, Speciesを列ベクトル単位で結合
  #Search Duplication
  data$Protein.IDs %>% duplicated() %>% any()
  data$Gene.names %>% duplicated() %>% any()
  data$Species %>% duplicated() %>% any()
  #Duplication table
  data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  data %>% group_by(Species) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
  #Unique Uniprot ID
  data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
  data_unique$Protein.IDs %>% duplicated() %>% any() # Are there any duplicated names?
  #SummarizedExperiment
  Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
  experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
  ###############################################################################
  #Log2-transform
  data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数, #Log2-transformation
  data1 <- data.frame(data_se@assays@data) #log2
  #Impute:left-shifted Gaussian distribution (for MNAR)
  data_imp_man <- impute(data_se, fun = "man", shift = 1.8, scale = 0.3) #Perseus,imputation
  data2 <- data.frame(data_imp_man@assays@data) #Subtract前log2imp
  #Subtract(Median):Perseus
  standardize <- function(z) {
    colmed <- apply(z, 2, median) #Median of Each Sample's Protein Expression level
    colmad <- apply(z, 2, mad)  # median absolute deviation
    rv <- sweep(z, 2, colmed,"-")  #subtracting median expression
    #rv <- sweep(rv, 2, colmad, "/")  # dividing by median absolute deviation
    return(rv)
  }
  data3 <- data2 #Subtract前log2impをコピー
  Sample_columns <- grep("(SC|PC)", colnames(data3)) # get Sample column numbers
  data3[Sample_columns] <- standardize(data3[Sample_columns]) #Subtract(Median),log2impsub
  #############################################################
  dat1 <- cbind(rownames(data1),data1) #log2
  dat2 <- cbind(rownames(data2),data2) #log2imp
  dat3 <- cbind(rownames(data3),data3) #log2impsub
  #integration
  dat <- cbind(data$Gene.names,data) #行名追加
  dat4 <- left_join(dat, dat1, by = c("Gene.names" = "rownames(data1)")) #raw+log2
  dat4 <- left_join(dat4, dat2, by = c("Gene.names" = "rownames(data2)")) #raw+log2+log2imp
  dat4 <- left_join(dat4, dat3, by = c("Gene.names" = "rownames(data3)")) #raw+log2+log2imp+log2impsub
  #output xlsx
  smp <- list("raw"=dat,"log2"=dat1,"log2imp"=dat2,"log2impsub"=dat3,"integ"=dat4,"anno"=anno) #リスト作成,rawdata,log2,imputation,subtract,integration
  write.xlsx(smp, "data.xlsx")
  #############################################################
  #statistic summary
  data_rm  <- data3
  data_rm[,1:2] <- NULL #列削除
  #transpose
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
  rownames(checkP) <- rownames(data3)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "GeneName"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(data3)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "GeneName"
  sdata <- cbind(sdata, checkQ)
  #############################################################
  #TukeyHSD function
  #diff群間の平均値の差(例)B-Aが-127.3であればデータBの平均がデータAの平均より-127.3大きい
  #lwr,upr=下方信頼限界,情報信頼限界:信頼区間の下限値 (lower) と上限値 (upper)
  #0を含まない場合 (例)B-A は含まず D-A は含む=2群間差は0ではないので有意差あり
  #p.adj < 0.05=2群間に有意差あり(信頼区間内に0を含まない)
  #############################################################
  THSD <- function(x) { 
    nn <- data.frame(P,C, x); 
    TukeyHSD(aov(x ~ P + C + P*C, nn))
  }
  THSDresults <- apply(data_rm, 1, THSD) 
  THSD_PC <- data.frame(lapply(THSDresults, function(x) {x["P:C"]}))
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
  THSDpd <- cbind(rownames(data3), THSDp_PC2, THSDd_PC2)
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
  sheets <- list("integ" = sdata3, "anovap" = checkPr2, 
                 "anovaq" = checkQr2, "THSDpd" = THSDpd2, 
                 "statvalue" = statv2) #assume sheet1-4 are data frames
  write_xlsx(sheets, "stat.xlsx", format_headers = FALSE)
  #############################################################
  #DEP list
  twANOVA_Pq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_P < 0.05)
  twANOVA_Cq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_C < 0.05)
  twANOVA_PxCq005 <- sdata3 %>% filter(Species == "MOUSE") %>% filter(q_PxC < 0.05)
  sheets2 <- list("Pq005"=twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Pq005))],
                  "Cq005"=twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_Cq005))],
                  "PxCq005"=twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN)", colnames(twANOVA_PxCq005))])
  write_xlsx(sheets2, "DEPtwANOVA.xlsx", format_headers = FALSE)
}
#fun1(data)
################################################################################
################################################################################