#perseus like analysis (3 arguments)
#Log2transform,Imputation(MNAR),Subtraction(Median),1wANOVA,2wANOVA,THSD
PLA2 <- function(x = read_excel("SWATH.xlsx", 2), 
                 y = read_excel("SWATH.xlsx", 3), 
                 z = read_excel("SWATH.xlsx", 4)){
  data <- x
  ExpDesign <- y
  group <- z
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
  #annotation
  dat <- cbind(data$Gene.names,data) #行名追加
  # annotation_columns <- dat[,c(1:7,80:82)]
  annotation_columns <- grep("(data$Gene.names|Index|Peak Name|m/z|Ret. Time|Group|Use|Protein.IDs|Gene.names|Species)", colnames(dat)) # get annotation column numbers
  annotation_columns <- dat[annotation_columns]
  dat1 <- cbind(annotation_columns,data1) #log2
  dat2 <- cbind(annotation_columns,data2) #log2imp
  dat3 <- cbind(annotation_columns,data3) #log2impsub
  #integration
  dat4 <- left_join(dat, dat1[,-c(1:7,9:12)], by = c("Protein.IDs" = "Protein.IDs")) #raw+log2
  dat4 <- left_join(dat4, dat2[,-c(1:7,9:12)], by = c("Protein.IDs" = "Protein.IDs")) #raw+log2+log2imp
  dat4 <- left_join(dat4, dat3[,-c(1:7,9:12)], by = c("Protein.IDs" = "Protein.IDs")) #raw+log2+log2imp+log2impsub
  #annotation
  dat <- left_join(dat, anno[,-c(1,2,4,5,7,9,10)], by = "Protein.IDs")
  dat1 <- left_join(dat1, anno[,-c(1,2,4,5,7,9,10)], by = "Protein.IDs")
  dat2 <- left_join(dat2, anno[,-c(1,2,4,5,7,9,10)], by = "Protein.IDs")
  dat3 <- left_join(dat3, anno[,-c(1,2,4,5,7,9,10)], by = "Protein.IDs")
  dat4 <- left_join(dat4, anno[,-c(1,2,4,5,7,9,10)], by = "Protein.IDs")
  #output xlsx
  #write.table(dat, "data_raw.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  #write.table(dat1, "data_log2.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  #write.table(dat2, "data_log2imp.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  #write.table(dat3, "data_log2impsub.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  #write.table(dat4, "data_integ.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  #write.table(anno, "data_anno.csv", quote=F, sep="," ,row.names=F, col.names=T, dec=".", append=F)
  smp <- list("raw"=dat,"log2"=dat1,"log2imp"=dat2,"log2impsub"=dat3,"integ"=dat4,"anno"=anno) #リスト作成,rawdata,log2,imputation,subtract,integration
  write_xlsx(smp, "data.xlsx", format_headers = FALSE)
  #############################################################
  #statistic summary
  data_rm  <- data3 #log2impsub
  rownames(data_rm) <- annotation_columns$Protein.IDs
  data_rm[,1:2] <- NULL #列削除
  #transpose
  tdata_rm <- t(data_rm)
  tdata_rm <- cbind(as.data.frame(rownames(tdata_rm)),tdata_rm)
  colnames(tdata_rm)[1] <- "ID"
  #grouping
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
  colnames(statSC0)[c(1,2)] <- c("condition","Protein.IDs")
  colnames(statSC10)[c(1,2)] <- c("condition","Protein.IDs")
  colnames(statSC30)[c(1,2)] <- c("condition","Protein.IDs")
  colnames(statPC0)[c(1,2)] <- c("condition","Protein.IDs")
  colnames(statPC10)[c(1,2)] <- c("condition","Protein.IDs")
  colnames(statPC30)[c(1,2)] <- c("condition","Protein.IDs")
  #bind
  statSC0 <- statSC0[,-1]
  statSC10 <- statSC10[,-1]
  statSC30 <- statSC30[,-1]
  statPC0 <- statPC0[,-1]
  statPC10 <- statPC10[,-1]
  statPC30 <- statPC30[,-1]
  statv2 <- left_join(statSC0, statSC10, by = "Protein.IDs")
  statv2 <- left_join(statv2, statSC30, by = "Protein.IDs")
  statv2 <- left_join(statv2, statPC0, by = "Protein.IDs")
  statv2 <- left_join(statv2, statPC10, by = "Protein.IDs")
  statv2 <- left_join(statv2, statPC30, by = "Protein.IDs")
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
  rownames(checkP) <- rownames(sdata)
  checkPr <- cbind(rownames(checkP),checkP)
  names(checkPr)[1] <- "Protein.IDs"
  #q値
  q_PC <- data.frame(p.adjust(p_PC, method = "BH"))
  q_P <- data.frame(p.adjust(p_P, method = "BH"))
  q_C <- data.frame(p.adjust(p_C, method = "BH"))
  q_PxC <- data.frame(p.adjust(p_PxC, method = "BH"))
  checkQ <- data.frame(cbind(q_PC, q_P, q_C, q_PxC))
  colnames(checkQ) <- c("q_PC", "q_P", "q_C","q_PxC") #rename
  rownames(checkQ) <- rownames(sdata)
  checkQr <- cbind(rownames(checkQ),checkQ)
  names(checkQr)[1] <- "Protein.IDs"
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
  THSDpd <- cbind(rownames(sdata), THSDp_PC2, THSDd_PC2)
  names(THSDpd)[1] <- "Protein.IDs"
  #############################################################
  #Annotation
  sdata2 <- cbind(rownames(sdata),sdata)
  names(sdata2)[1] <- "Protein.IDs"
  sdata2 <- left_join(sdata2, statv2, by = "Protein.IDs")
  sdata2 <- left_join(sdata2, THSDpd, by = "Protein.IDs")
  sdata3 <- left_join(sdata2, anno, by = "Protein.IDs")
  checkPr2 <- left_join(checkPr, anno, by = "Protein.IDs")
  checkQr2 <- left_join(checkQr, anno, by = "Protein.IDs")
  THSDpd2 <- left_join(THSDpd, anno, by = "Protein.IDs")
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
  sheets2 <- list("Pq005"=twANOVA_Pq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN|ENSEMBL|ENTREZID|GENENAME|MGI)", colnames(twANOVA_Pq005))],
                  "Cq005"=twANOVA_Cq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN|ENSEMBL|ENTREZID|GENENAME|MGI)", colnames(twANOVA_Cq005))],
                  "PxCq005"=twANOVA_PxCq005[,grep("(GeneName|p_P$|p_C$|p_PxC$|q_P$|q_C$|q_PxC$|Protein.IDs|Description|GN|ENSEMBL|ENTREZID|GENENAME|MGI)", colnames(twANOVA_PxCq005))])
  write_xlsx(sheets2, "DEPtwANOVA.xlsx", format_headers = FALSE)
}
