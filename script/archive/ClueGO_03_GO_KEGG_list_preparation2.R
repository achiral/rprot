#ClueGO_03_GO_KEGG_list_preparation2.R
#R dplyrパッケージで複数の列を文字列として指定し結合された列を追加する方法
#https://www.trifields.jp/how-to-add-multiple-joined-columns-by-specifying-multiple-columns-as-strings-in-dplyr-package-in-r-2812
#https://ja.stackoverflow.com/questions/46965/rでデータフレーム内の二列を結合し新しい列を作りたい
#文字列結合
#https://rpubs.com/hoxo_m/5731
#列の変換
#https://kazutan.github.io/kazutanR/hands_on_170730/mutate.html#複数の条件で場合分け
#https://ando-roid.hatenablog.com/entry/2019/04/17/124842
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/SWATH")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#install.packages("xlsx") #エクセル出力
#install.packages("xlsx2") #エクセル出力,動作しない
library(tidyverse) #ライブラリ
#library(dplyr)
library(readxl) #エクセル入力
library(xlsx) #エクセル出力
library(openxlsx) #write.xlsxによるエクセル出力
#library(xlsx2) #エクセル出力,動作しない
library(rJava)

#Excelファイル読み込み
GO_top2 <- read_excel("GO_KEGG.xlsx", 1) #シート1の読み込み
KEGG_top2 <- read_excel("GO_KEGG.xlsx", sheet = 2) #シート2の読み込み

#GO_top2 <- read_excel("GO_top2.xlsx", 1) #シート1の読み込み
#KEGG_top2 <- read_excel("KEGG_top2.xlsx", sheet = 1) #シート1の読み込み

#列の作成と結合
GO_top2 <- GO_top2 %>% as.data.frame() %>%
  mutate(DBsourse = case_when(
    `Ontology Source.x` == "GO_BiologicalProcess-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "BP",
    `Ontology Source.x` == "GO_CellularComponent-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "CC",
    `Ontology Source.x` == "GO_MolecularFunction-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "MP",
    `Ontology Source.x` == "KEGG_08.05.2020" ~ "KEGG", TRUE ~ as.character(`Ontology Source.x`)
    ))%>%
  mutate(GOTerm.n = paste(!!!rlang::syms(c("GOTerm.x", "DBsourse")), sep="_")) #シンボル化：map(as.symbol)と同じ：paste(!!!rlang::syms(c("列名", "列名", "列名")), sep="_")
KEGG_top2 <- KEGG_top2 %>% as.data.frame() %>%
  mutate(DBsourse = case_when(
    `Ontology Source.x` == "GO_BiologicalProcess-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "BP",
    `Ontology Source.x` == "GO_CellularComponent-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "CC",
    `Ontology Source.x` == "GO_MolecularFunction-EBI-UniProt-GOA-ACAP-ARAP_08.05.2020_00h00" ~ "MP",
    `Ontology Source.x` == "KEGG_08.05.2020" ~ "KEGG", TRUE ~ as.character(`Ontology Source.x`)
  ))%>%
  mutate(GOTerm.n = paste(!!!rlang::syms(c("GOTerm.x", "DBsourse")), sep="_")) #シンボル化：map(as.symbol)と同じ：paste(!!!rlang::syms(c("列名", "列名", "列名")), sep="_")

#出力
#install.packages("openxlsx") 
library(openxlsx) #write.xlsxによるエクセル出力
smp <- list("GO"=GO_top2, "KEGG"=KEGG_top2) #GOリスト,KEGGリスト作成
write.xlsx(smp, "GO_KEGG.xlsx") #GOシート,KEGGシート出力

#エクセル出力
write.xlsx(GO_top2, "GO_top2.xlsx", sheetName="sheet1") #出力
write.xlsx(KEGG_top2, "KEGG_top2.xlsx", sheetName="sheet2") #出力

quit() #Rを閉じてシートを再読込してグラフを描写する