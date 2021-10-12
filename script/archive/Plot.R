#Plot
#############################################################
#ライブラリ読み込み
#library(DEP)
library(tidyverse)
library(readxl) #エクセル入力(read_excel)
#library(xlsx) #エクセル入力
library(openxlsx) #エクセル入出力(write.xlsx)
#library(writexl) #xlsx出力
#library(multcomp) #多重比較検定
library(devtools)
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")
library(reshape2)
#############################################################
#rm(list = ls(all.names = TRUE))
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/Other") #作業ディレクトリ設定
#anno <- read_excel("anno.xlsx", 1) #シート1入力
#anno2 <- anno %>% select(GeneName, `Peak Name`, Description, OS, GN, PE, SV)
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc") #作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
##############################################################
#xlsx入力
dat1 <- read_excel("stat.xlsx", 1)
#PCPxCLZ
dat1_f <- dat1 %>% 
  filter(Species == "MOUSE") %>%
  filter(q_PxC < 0.05) %>%
  filter(THSDp_P.C0.S.C0 < 0.05) %>%
  filter(THSDp_P.C10.P.C0 < 0.05|THSDp_P.C30.P.C0 < 0.05) #%>%
  #filter(abs(diff_P.C0.S.C0) > 0.67) 
#PCP
dat1_f <- dat1 %>% 
  filter(Species == "MOUSE") %>%
  filter(q_P < 0.05) %>%
  filter(abs(diff_P.C0.S.C0) > 0.67) 
#CLZ
dat1_f <- dat1 %>% 
  filter(Species == "MOUSE") %>%
  filter(q_C < 0.05) %>%
  filter(abs(diff_S.C10.S.C0) > 1|abs(diff_S.C30.S.C0) > 1)
#intersect
dat1_f <- dat1 %>% 
  filter(Species == "MOUSE") %>%
  filter(q_P < 0.05) %>%
  filter(q_C < 0.05) %>%
  filter(q_PxC < 0.05)

#top絞り込み
top_10 <- dat1_f %>% arrange(desc(THSDp_P.C0.S.C0)) %>% group_by(Species) %>% slice(1:10)
top_50 <- dat1_f %>% arrange(desc(THSDp_P.C0.S.C0)) %>% group_by(Species) %>% slice(1:50)
dat1_f <- top_10
dat1_f <- top_50

#列削除
dat1_f <- dat1_f[,-grep("^p_|^q_|^THSD|^diff_|^Peak|^Group|^Protein|^Species|^Description|^OS|^GN|^PE|^SV", colnames(dat1_f))]
t(colnames(dat1_f))
#行列入れ替え
tdat1 <- t(dat1_f)
#列名入力
colnames(tdat1) <- (tdat1[1,])
tdat1 <- as.data.frame(tdat1[-1,])
tdat1 <- rownames_to_column(tdat1, var = "ID")
#csv出力
write_csv(tdat1, "tstat.csv")
#xlsx出力
#smp <- list("tstat"=tdat1)
#write.xlsx(smp, "tstat.xlsx", row.names=T)
#csv入力
tdat1 <- read_csv("tstat.csv", col_types = NULL)
#xlsx入力
#tdat1 <- as.data.frame(read_excel("tstat.xlsx", 1))
#colnames(tdat1)[1] <- "ID"
#アノテーション
group <- read_excel("SWATH.xlsx", 4)
cond <- group[,grep("ID|condition", colnames(group))]
#left_join
tdat2 <- left_join(tdat1, cond, by = c("ID"="ID"))
tdat2 <- tdat2[,-1]
##############################################################
#データ形式確認
#typeof(tdat2)
#mode(tdat2)
#class(tdat2)
#class(tdat2$AT1A2)
#if(mode(tdat2) == "numeric") {print("tdat_m$value is numeric.")} else {print("tdat_m$value is not numeric.")}
#if (is.numeric(tdat2)) {print("tdat_m$value is numeric.")} else {print("tdat_m$value is not numeric.")}
tdat_m = reshape2::melt(tdat2,
                        id.vars="condition",
                        variable.name="Protein", na.rm=T)
#tdat_m <- drop_na(tdat_m, everything())
#無意味？as.factor(tdat_m$condition)
#無意味？levels(tdat_m$condition)
tdat_m2 <- transform(tdat_m, condition = factor(condition, levels = c("SC0", "SC10", "SC30", "PC0", "PC10", "PC30")))
levels(tdat_m2$condition)
x <- arrange(as.data.frame(levels(tdat_m2$Protein)),`levels(tdat_m2$Protein)`)
y <- unlist(x)
tdat_m2 <- transform(tdat_m2, Protein = factor(Protein, levels = y))
levels(tdat_m2$Protein)
##############################################################
#群選択
SC0 <- tdat_m %>% filter(condition == "SC0")
SC10 <- tdat_m %>% filter(condition == "SC10")
SC30 <- tdat_m %>% filter(condition == "SC30")
PC0 <- tdat_m %>% filter(condition == "PC0")
PC10 <- tdat_m %>% filter(condition == "PC10")
PC30 <- tdat_m %>% filter(condition == "PC30")
tdat_m2 <- rbind(SC0, PC0)
tdat_m2 <- rbind(SC0, PC0, PC10, PC30)
tdat_m2 <- rbind(SC0, SC10, SC30)
tdat_m2 <- rbind(PC0, PC10, PC30)
#tdat_m2 <- tdat_m %>% filter(condition == c("SC0", "PC0", "PC10", "PC30"))
#無意味？tdat_m3 <- arrange(tdat_m2, as.character(Protein))
tdat_m2 <- transform(tdat_m2, condition = factor(condition, levels = c("SC0", "PC0", "PC10", "PC30")))
levels(tdat_m2$condition)
x <- arrange(as.data.frame(levels(tdat_m2$Protein)),`levels(tdat_m2$Protein)`)
y <- unlist(x)
tdat_m2 <- transform(tdat_m2, Protein = factor(Protein, levels = y))
levels(tdat_m2$Protein)
##############################################################
t<-proc.time()
plot1 <- ggplot(data = tdat_m2, mapping = aes(x = condition, y = value, fill = condition)) +
  geom_flat_violin(scale = "count", trim = F) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               position = position_nudge(0.05),
               colour = "black",
               size = 0.1) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.5, 
               stackdir = "down", 
               binwidth = 0.1, 
               position = position_nudge(-0.025))+ 
  facet_wrap(~Protein, ncol = 5, scales = "free")+
  theme(legend.position = "none") + #theme_bw()
  labs(x = "condition", y = "value")
proc.time()-t
# 描画して確認
#t<-proc.time()
print(plot1)
#proc.time()-t
# PNG画像として保存
t<-proc.time()
ggsave(file = "plot1.png", plot = plot1, dpi = 300, width = 20, height = 20)
proc.time()-t
# PDF画像として保存
t<-proc.time()
ggsave(file = "plot1.pdf", plot = plot1, width = 28, height = 21.2)
proc.time()-t
##############################################################
q()
