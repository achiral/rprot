#Annotation
#biomaRt（https://ukulele.nig.ac.jp/SurfWiki/biomaRt.html）（https://qiita.com/yuifu/items/a757629506c1cd98156b）
###########################################################################
rm(list = ls(all = TRUE))
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", 
                 "package:utils", "package:datasets", "package:methods", "package:base")
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  pkg.list <- setdiff(pkg.list, basic.pkg)
  lapply(pkg.list, detach, character.only = TRUE)
}
detach_all()
###########################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db", "mouse4302.db","GO.db", "PANTHER.db", "biomaRt"))
###########################################################################
#改善の余地が大いにあり
###########################################################################
#library(tidyverse)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(GO.db)
#library(PANTHER.db)
#library(biomaRt)

#library(EnhancedVolcano)
#library(magrittr)
#library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
#library(dplyr)
#library(scales) #muted()関数使用のため
#library(gridExtra) #svg出力のため
#library(rJava)
#library(readxl) #エクセル読み込み
#library(openxlsx) #JAVA不使用で大きなデータも読み込める
#library(cowplot)
###########################################################################
#CSV入力
library(readxl) #エクセル読み込み
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/Other")
anno <- read_excel("anno.xlsx", 1) #シート1入力
#library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
id <- anno$`Protein.IDs`
GN <- anno$GN
#GeneName <- anno$GeneName
###########################################################################
###########################################################################
#Bioconductor
#columns()：列名
#t(columns(org.Hs.eg.db))
#t(columns(org.Mm.eg.db))
#col <- columns(org.Mm.eg.db)
#t(columns(GO.db))
#columns(PANTHER.db)
#keytypes()：検索条件の指定に用いる項目
#t(keytypes(org.Hs.eg.db))
#t(keytypes(org.Mm.eg.db))
#t(keytypes(GO.db))
###########################################################################
#生物種レベルのアノテーション（OrgDb）
#ENSid <- c("ENSG00000004660","ENSG00000162946")
#res <- select(org.Hs.eg.db,
#              keys = c("ENSG00000004660","ENSG00000162946"),
#              keytype = "ENSEMBL",
#              columns = c("ENSEMBL", "ENTREZID", "GO", "PATH"))
res_id <- select(org.Mm.eg.db, keys = id, keytype = "UNIPROT",
              columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN <- select(org.Mm.eg.db, keys = GN, keytype = "SYMBOL",
              columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))


###########################################################################
#left_join, remove duplicates
library(tidyverse)
#res <- inner_join(res_id, res_GN, by = c ("UNIPROT"))
#anno2 <- left_join(anno, res, by = c("Protein.IDs" = "UNIPROT", "GN" = "SYMBOL"))
#anno2 <- left_join(anno, res, by = c(("Protein.IDs" = "UNIPROT", "GN" = "SYMBOL")) %>% distinct(Protein.IDs, .keep_all = T))

#search Duplicates
#res_id$UNIPROT %>% duplicated() %>% any()
#res_GN$UNIPROT %>% duplicated() %>% any()
#anno$Protein.IDs %>% duplicated() %>% any()
#anno2$Protein.IDs %>% duplicated() %>% any()

#duplication table
#res_id_dup <- res_id %>% group_by(UNIPROT) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
#res_GN_dup <- res_GN %>% group_by(UNIPROT) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
#anno2_dup <- anno2 %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)

#semi_join
semi1 <- rbind(semi_join(res_id, res_GN, by = "UNIPROT"), semi_join(res_id, res_GN, by = "SYMBOL"))
semi2 <- rbind(semi_join(res_GN, res_id, by = "UNIPROT"), semi_join(res_GN, res_id, by = "SYMBOL"))
semi2 <- semi2[,c(6,2,3,4,5,1)]

#anti_join
anti1 <- rbind(anti_join(res_id, res_GN, by = "UNIPROT"), anti_join(res_id, res_GN, by = "SYMBOL"))
anti2 <- rbind(anti_join(res_GN, res_id, by = "UNIPROT"), anti_join(res_GN, res_id, by = "SYMBOL"))
anti2 <- anti2[,c(6,2,3,4,5,1)]

#rbind
semi_anti1 <- rbind(semi1, anti1)
semi_anti2 <- rbind(semi2, anti2)
semi_anti3 <- rbind(semi_anti1, semi_anti2)

#remove duplicates
#ex1_id <- semi_anti1 %>% distinct(UNIPROT, .keep_all = T)
#ex2_id <- semi_anti2 %>% distinct(UNIPROT, .keep_all = T)
ex3_id <- semi_anti3 %>% distinct(UNIPROT, .keep_all = T)
#ex1_GN <- semi_anti1 %>% distinct(SYMBOL, .keep_all = T)
#ex2_GN <- semi_anti2 %>% distinct(SYMBOL, .keep_all = T)
ex3_GN <- semi_anti3 %>% distinct(SYMBOL, .keep_all = T)

#left_join
anno_id <- left_join(anno, ex3_id, by = c("Protein.IDs" = "UNIPROT")) 
anno_id2 <- anno_id %>% filter(!is.na(ENTREZID))
#anno_id2 <- anno_id2 %>% distinct(ENTREZID, .keep_all = T)
anno_GN <- left_join(anno, ex3_GN, by = c("GN" = "SYMBOL"))
anno_GN2 <- anno_GN %>% filter(!is.na(ENTREZID))
semi4 <- semi_join(anno_id2, anno_GN2, by = "Protein.IDs")
semi5 <- semi_join(anno_GN2, anno_id2, by = "Protein.IDs")
semi6 <- rbind(semi4[,1:14], semi5[,-11])
semi6 <- semi6 %>% distinct(ENTREZID, .keep_all = T)
t(colnames(semi6))
semi7 <- semi6[,c(3,11:14)]
semi7 <- semi7 %>% distinct(Protein.IDs, .keep_all = T)
anno_n <- left_join(anno, semi7, by = "Protein.IDs")
anno_n <- anno_n %>% distinct(Protein.IDs, .keep_all = T)
anno_n1 <- anno_n %>% filter(!is.na(ENTREZID))

#NA value
sub <- anno_n %>% subset(is.na(anno_n$ENTREZID))
t(colnames(sub))
t(colnames(ex3_GN))
anno_n2 <- left_join(sub[,1:10], ex3_GN, by = c("GN" = "SYMBOL"))
anno_n2 <- anno_n2 %>% distinct(Protein.IDs, .keep_all = T)

#NA value
sub2 <- anno_n2 %>% subset(is.na(anno_n2$ENTREZID))
sub2 <- sub2[,1:10]
sub2_f <- sub2 %>% filter(Species == "MOUSE")
detach_all()
library(org.Mm.eg.db)
id2 <- sub2_f$`Protein.IDs`
GN2 <- sub2_f$`GN`
res_id2 <- select(org.Mm.eg.db, keys = id2, keytype = "UNIPROT",
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN2 <- select(org.Mm.eg.db, keys = GN2, keytype = "SYMBOL",
                 columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
library(tidyverse)
res_id2 <- res_id2 %>% filter(!is.na(ENTREZID)) %>% distinct(UNIPROT, .keep_all = T)
anno_n3 <- left_join(sub2_f, res_id2, by = c("Protein.IDs" = "UNIPROT"))
anno_n3 <- anno_n3 %>% distinct(Protein.IDs, .keep_all = T)

#subset including NA
sub3 <- subset(anno_n3, is.na(anno_n3$ENTREZID))
t(colnames(sub))
sub3 <- sub3[,1:10]

#remove libraries
detach_all()
library(org.Mm.eg.db)
#entrezID
ent <- c("18563", "234695", "14467", "14070")
res_ent <- select(org.Mm.eg.db, keys = ent, keytype = "ENTREZID",
                  columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
#remove duplicates
library(tidyverse)
res_ent <- res_ent %>% filter(!is.na(ENTREZID)) %>% distinct(ENTREZID, .keep_all = T)

#cbind
res_ent[1,]
res_ent <- res_ent[,c(2,1,3,4,5)]
res_ent[1,]
t(colnames(sub3))
anno_n3 <- cbind(sub3, res_ent)

t(colnames(anno_n))
t(colnames(anno_n2))
anno_n2 <- anno_n2[,-11]
t(colnames(anno_n3))
anno_n3 <- anno_n3[,-15]

#remove 4 rows includes above 4 proteins 
grep("(Q05920|Q3V3V9|O55126|Q00558)", anno_n2$Protein.IDs)
anno_n2 <- anno_n2[-grep("(Q05920|Q3V3V9|O55126|Q00558)", anno_n2$Protein.IDs),]
grep("(Q05920|Q3V3V9|O55126|Q00558)", anno_n2$Protein.IDs)
anno_nn <- rbind(anno_n1, anno_n2, anno_n3)

#list
list_ok <- anno_nn %>% filter(!is.na(ENTREZID))
list_ng_mouse <- anno_nn %>% subset(is.na(anno_nn$ENTREZID)) %>% filter(Species == "MOUSE")
list_ng_other <- anno_nn %>% subset(is.na(anno_nn$ENTREZID)) %>% filter(!Species == "MOUSE")

id3 <- list_ng_mouse$`Protein.IDs`
GN3 <- list_ng_mouse$GN

#remove libraries
detach_all()
library(org.Mm.eg.db)
res_id3 <- select(org.Mm.eg.db, keys = id3, keytype = "UNIPROT",
                  columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))
res_GN3 <- select(org.Mm.eg.db, keys = GN3, keytype = "SYMBOL",
                  columns = c("ENSEMBL", "ENTREZID", "GENENAME", "MGI", "SYMBOL", "UNIPROT"))





#remove duplicates
library(tidyverse)
res_id3 <- res_id3 %>% filter(!is.na(ENTREZID)) %>% distinct(ENTREZID, .keep_all = T) %>% distinct(UNIPROT, .keep_all = T)
list_ng_mouse <- left_join(list_ng_mouse[,1:10], res_id3, by = c("Protein.IDs" = "UNIPROT"))

t(colnames(list_ok))
t(colnames(list_ng_mouse))
list_ng_mouse <- list_ng_mouse[,-15]
t(colnames(list_ng_other))
anno_nnn <- rbind(list_ok, list_ng_mouse, list_ng_other)

#output xlsx
library(openxlsx) #入出力(write.xlsx)
smp <- list("anno"=anno,"anno_new"=anno_nnn)
write.xlsx(smp, "anno2.xlsx")
###########################################################################
###########################################################################
###########################################################################

#システム生物学レベルのアノテーション（GO.db）
goid <- c("GO:0004683", "GO:0005524", "GO:0007268")
res <- select(GO.db,
              keys = goid,
              keytype = "GOID",
              columns = c("GOID", "TERM", "ONTOLOGY"))
###########################################################################
db <- useMart("ensembl")
head(listMarts())
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
head(listFilters(hg))
head(listAttributes(hg))
#染色体名、遺伝子の開始位置と終了位置を検索
ensid <- c("ENSG00000000003", "ENSG00000000457", "ENSG00000000938",
           "ENSG00000006327", "ENSG00000011405", "ENSG00000001497")
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position"),
             filters = "ensembl_gene_id", 
             values = ensid,
             mart = hg)
res
str(ensid)
#Ensembl ID から GO term を検索する
ensid <- c("ENSG00000000003", "ENSG00000000457", "ENSG00000000938",
           "ENSG00000006327", "ENSG00000011405", "ENSG00000001497")
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", "go_id"),
             filters = "ensembl_gene_id", 
             values = ensid,
             mart = hg)
head(res)
tail(res)
#Ensembl ID から Entrez 遺伝子 ID を検索する
ensid <- c("ENSG00000000003", "ENSG00000000457", "ENSG00000000938",
           "ENSG00000006327", "ENSG00000011405", "ENSG00000001497")
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "go_id"),
             filters = "ensembl_gene_id", 
             values = ensid,
             mart = hg)
head(res)
###########################################################################
#Ensembl
#マートとデータセットの決定
library(biomaRt)
ensembl = useMart("ensembl", "scerevisiae_gene_ensembl")
listMarts()
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("scerevisiae_gene_ensembl", ensembl)
#どんなデータやフィルタが利用可能か調べる
attributePages(ensembl)
listAttributes(ensembl, "feature_page")
subset(ensembl@filters, select=c(name, description, type, operation))
#フィルタの選択肢を調べる
subset(ensembl@filters, name=="biotype")$options
keys(ensembl, "biotype")
keys(ensembl, "go_evidence_code")
#近いミラーを使う
listMarts(host="asia.ensembl.org")
ensembl = useMart("ENSEMBL_MART_SNP", "scerevisiae_snp", host="asia.ensembl.org")
###########################################################################
#UniProt
unimart = useMart("unimart", "uniprot")
unimart = useMart("unimart")
unimart = useMart("uniprot")
listDatasets(uniprot)
subset(unimart@filters, select=c(name, description, type, operation))
#フィルタの選択肢
subset(unimart@filters, 3 < nchar(options) & nchar(options) < 120, select=c(name, options))
#“Complete proteome” の選択肢(すげえ長い)を抜き出す
proteome_name = biomaRt::keys(unimart, "proteome_name")
grep("Sac.* cer.*", proteome_name, value=TRUE)
#アトリビュート列挙
listAttributes(unimart)



###########################################################################
#biomaRtバージョン情報
packageVersion("biomaRt")
#接続できるBioMartデータベースの種類
list <- listMarts()
#Ensemblデータベースにアクセス
ensembl <- useMart("ensembl")
#table of dataset
dataset <- as.data.frame(listDatasets(ensembl))
#select dataset
#ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
#table of filters
filters <- as.data.frame(listFilters(ensembl))
filters1 <- c("ensembl_gene_id", "entrezgene_id")


#table of attributes
attributes <- as.data.frame(listAttributes(ensembl))
attributes1 <- c("ensembl_gene_id", "entrezgene", "go_id")
###########################################################################
#biomaRt実践
db <- useMart("ensembl") #ensembl
mg <- useDataset("mmusculus_gene_ensembl", mart = db) #mouse
#id <- c() #protein_id




res <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "go_id"),
             filters = "with_uniprot_gn", 
             values = "Dlg4",
             mart = mg)
head(res)

id <- iidd[[1]][["Protein.IDs"]]
str(id)
idd <- as.data.frame(id)
iddd <- as.vector(idd)

str(iddd)
iidd <- list(anno[,"Protein.IDs"])

?getBM()
