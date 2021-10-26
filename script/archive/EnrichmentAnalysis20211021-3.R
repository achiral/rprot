# enrichment analysis
# clusterProfiler(https://yulab-smu.top/biomedical-knowledge-mining-book/index.html)
# Repidemiology(https://jojoshin.hatenablog.com/entry/2016/09/27/161744)
# GO, KEGG(https://shiokoji11235.com/go-analysis-for-transcriptome-analysis)
# GSEA(https://www.ncc.go.jp/jp/ri/division/rare_cancer_research/labo/20181011150525.html)
# GSEA(https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#abstract)
# GSEA(https://www.biostars.org/p/445059/)
# reactome pathway(https://shiokoji11235.com/pathway-enrichment-analysis-in-r)
# enrichment map trouble shooting(https://github.com/YuLab-SMU/enrichplot/issues/79)
################################################################################
## 1. instal/load pkg #####
################################################################################
# BiocManager::install("clusterProfiler")
# install.packages("ggnewscale")
# require(AnnotationHub)           # BioC 3.14 (Nov. 2021, with R-4.2.0)
require(clusterProfiler)
require(DOSE)
require(europepmc)
require(ggnewscale)           # enrichment map plot
require(ggplot2)
require(ggupset)
require(Hmisc)
# require(MeSHDbi)                 # BioC 3.14 (Nov. 2021, with R-4.2.0)
require(meshes)
# require(MeSH.Hsa.eg.db)          # BioC 2.14-3.13
require(MeSH.Mmu.eg.db)            # BioC 2.14-3.13
require(msigdbr)
# require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(pathview)
require(ReactomePA)
require(stringr)
################################################################################
## 2. setup #####
################################################################################
rm(list = ls(all = TRUE))
setwd("~/Dropbox/GitHub/local/Docker/R/rprot")
# source("/home/rstudio/rproject/script/archive/functions.R")        # load functions
# source("/home/rstudio/rproject/script/archive/functions_DEPpkg.R") # load DEPpkg functions
source("script/startup.R")                  # load packages
source("script/archive/functions.R")
detach_all()
source("script/startup.R")                  # load packages
source("script/archive/functions.R")
################################################################################
## 3. set working directory #####
################################################################################
getwd()
# setwd("/home/rstudio/rproject")
# setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")
# setwd("~/Dropbox/GitHub/local/Docker/R/rprot")

setwd("~/Dropbox/GitHub/local/Docker/R/rprot/data/Perseus_Like_Analysis/PFC")
DIR_HOST <- getwd()  # host working directory
# setwd(DIR_HOST)

## input data set
dat <- read_excel("stat.xlsx", 1)
################################################################################
## 4. prepare geneList #####
################################################################################
# t(colnames(dat))                                            # check colnames
# d <- dat %>% dplyr::select(ENTREZID,diff_P.C0.S.C0)         # column1=ID, column2=FC
# d <- dat[,grep("ENTREZID|diff_P.C0.S.C0|diff_*", colnames(dat))] # x column1=ID, column2=FC
# dat_rm <- dat %>% drop_na(ENTREZID)
# dat_rm$ENTREZID
## twANOVAq < 0.05
## PCP
dat_p <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_P < 0.05) %>% 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## CLZ
dat_c <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_C < 0.05) %>% 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## PCPxCLZ(interaction)
dat_pc <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## twANOVAq_PxC(interaction) < 0.05 & THSDp < 0.05
## PCP
dat_ps <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05) %>%
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## attenuate(improve) PCP reaction
## PCP up/CLZ down
dat_ud <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 > 0) %>%                                                    # PCP -> up-regulate
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 < 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 < 0) %>%  # CLZ -> down-regulate 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## PCP down/CLZ up
dat_du <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 < 0) %>%                                                    # PCP -> down-regulate
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 > 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 > 0) %>%  # CLZ -> up-regulate 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## amplify(exacerbate) PCP reaction
## PCP up/CLZ up
dat_uu <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 > 0) %>%                                                    # PCP -> up-regulate
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 > 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 > 0) %>%  # CLZ -> up-regulate
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()
## PCP down/CLZ down
dat_dd <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 < 0) %>%                                                    # PCP -> down-regulate
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 < 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 < 0) %>%  # CLZ -> down-regulate 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()

## geneLiset #####
geneList_p <- prepare_geneList(dat_p)
geneList_c <- prepare_geneList(dat_c)
geneList_pc <- prepare_geneList(dat_pc)
geneList_ps <- prepare_geneList(dat_ps)
geneList_ud <- prepare_geneList(dat_ud)
geneList_du <- prepare_geneList(dat_du)
geneList_uu <- prepare_geneList(dat_uu)
geneList_dd <- prepare_geneList(dat_dd)
GENELIST <- list(geneList_p, geneList_c, geneList_pc, geneList_ps,
                 geneList_ud, geneList_du, geneList_uu, geneList_dd)

gene_p <- names(geneList_p)
gene_c <- names(geneList_c)
gene_pc <- names(geneList_pc)
gene_ps <- names(geneList_ps)
gene_ud <- names(geneList_ud)
gene_du <- names(geneList_du)
gene_uu <- names(geneList_uu)
gene_dd <- names(geneList_dd)
GENE <- Hmisc::llist(gene_p, gene_c, gene_pc, gene_ps,
                     gene_ud, gene_du, gene_uu, gene_dd)
################################################################################
## set parameters #####
################################################################################
## expression data #####
DAT_GENE <- dat_p
# DAT_GENE <- dat_c
# DAT_GENE <- dat_pc
# DAT_GENE <- dat_ps
# DAT_GENE <- dat_ud
# DAT_GENE <- dat_du
# DAT_GENE <- dat_uu
# DAT_GENE <- dat_dd

## visualization parameters #####
geneList <- prepare_geneList(DAT_GENE)
gene <- names(geneList)[abs(geneList) > 0]
ont <- "BP"                                     # "BP"(default ),"CC","MF","ALL"
showCategory <- 20                              # default = 20
layout = "kk"                                   # "kk", "nicely"
labeln = "all"                                  # "category", "gene", "all", "none"
db <- org.Mm.eg.db
d <- GOSemSim::godata(OrgDb = db, ont = ont)    # emap plot
MeSHDb <- "MeSH.Mmu.eg.db"                      # MeSHDb <- MeSH.Hsa.eg.db    # BioC 2.14-3.13
################################################################################
## 5. GO enrichment analysis #####
################################################################################
ego <- run_enrichGO(DAT_GENE, db = db, ont = ont)  # modify GO term: "BP","CC","MF","ALL"
egosim <- clusterProfiler::simplify(ego)      # remove redundancy
################################################################################
## 6. KEGG pathway analysis #####
################################################################################
# KEGG pathway over-representation analysis #####
ek <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',       # 'hsa'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
## ID conversion(http://guangchuangyu.github.io/2016/05/convert-biological-id-with-kegg-api-using-clusterprofiler/)
eg2np <- bitr_kegg(gene, fromType='kegg', toType='ncbi-proteinid', organism='mmu') # Path, Module, ncbi-proteinid, ncbi-geneid, uniprot, kegg
eg2up <- bitr_kegg(eg2np[,2], fromType='ncbi-proteinid', toType='uniprot', organism='mmu') # Path, Module, ncbi-proteinid, ncbi-geneid, uniprot, kegg
## KEGG pathway enrichment analysis
ek2 <- enrichKEGG(gene = eg2up[,2],   # gene = eg2np[,2]
                  organism = 'mmu',
                  keyType='uniprot',  # keyType='ncbi-proteinid'
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1)
## ID conversion #####
ek3 <- setReadable(ek2, OrgDb = db, keyType="UNIPROT")
# head(summary(ek3))

## KEGG module over-representation analysis
mek <- enrichMKEGG(gene = gene,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
mek2 <- enrichMKEGG(gene = eg2up[,2],   # gene = eg2np[,2]
                    organism = 'mmu',
                    keyType='uniprot',  # keyType='ncbi-proteinid'
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
## ID conversion
mek3 <- setReadable(mek2, OrgDb = db, keyType="UNIPROT")
################################################################################
## 7. WikiPathways analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html) #####
################################################################################
# get_wp_organisms()
## WikiPathways over-representation analysis
ewp <- enrichWP(gene, organism = "Mus musculus")
ewp2 <- setReadable(ewp, OrgDb = db, keyType="ENTREZID")
################################################################################
## 8. reactome pathway analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html) #####
################################################################################
er <- enrichPathway(gene = gene,                   # entrez gene id
                    organism = "mouse",            # "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
                    pvalueCutoff = 0.05,
                    # pAdjustMethod = "BH",        # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                    # qvalueCutoff = 0.05,
                    # universe = dat_rm$ENTREZID,  # background genes
                    # minGSSize = 3,               # 10
                    # maxGSSize = 500,
                    readable = TRUE)                # gene ID -> gene Name
################################################################################
## 9. Disease enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html) #####
################################################################################
# only for human(see bottom)
################################################################################
## 10. MeSH enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/meshes-enrichment.html) #####
################################################################################
## MeSH over-representation analysis
# str(geneList)
# de <- names(geneList)[1:100]
de <- names(geneList)
emesh <- enrichMeSH(de, MeSHDb = MeSHDb, 
                    database = 'gendoo',      # 'gendoo', 'gene2pubmed' or 'RBBH'
                    category = 'C',           # "A", "B", "C"(Diseases), "D", "E", "F", "G"(Phenomena and Processes), "H", "I", "J", "K", "L","M", "N", "V", "Z"
                    # universe,
                    # minGSSize = 10, maxGSSize = 500,
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH",     # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                    qvalueCutoff = 1)
emesh2 <- setReadable(emesh, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
################################################################################
## 11. Universal enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html) #####
################################################################################
## 11. Cell Marker Data #####
################################################################################
## cell_marker_data #####
cmd <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt')
# cells <- cell_marker_data %>%                          # error
#   dplyr::select(cellName, geneID) %>%
#   dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
#   tidyr::unnest()
cmd2 <- dplyr::select(cmd, cellName, geneID)
cmd3 <- dplyr::mutate(cmd2, geneID = strsplit(geneID, ', '))
cells <- tidyr::unnest(cmd3)
## Cell Marker over-presentaton analysis
ecm <- enricher(gene, TERM2GENE = cells)
ecm2 <- setReadable(ecm, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
################################################################################
## 12. MSigDb analysis #####
################################################################################
## check species
# msigdbr_species()     # msigdbr_show_species()
# m_df <- msigdbr(species = "Mus musculus")
# head(m_df, 2) %>% as.data.frame

## MSigDb gene sets
H_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, entrez_gene)    # H(hallmark gene sets)
C1_t2g <- msigdbr(species = "Mus musculus", category = "C1") %>% dplyr::select(gs_name, entrez_gene)  # C1(positional gene sets)
C2_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::select(gs_name, entrez_gene)  # C2(curated gene sets)
C3_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% dplyr::select(gs_name, entrez_gene)  # C3(motif gene sets)
C4_t2g <- msigdbr(species = "Mus musculus", category = "C4") %>% dplyr::select(gs_name, entrez_gene)  # C4(computational gene sets)
C5_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% dplyr::select(gs_name, entrez_gene)  # C5(GO gene sets)
C6_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% dplyr::select(gs_name, entrez_gene)  # C6(oncogenic gene sets)
C7_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% dplyr::select(gs_name, entrez_gene)  # C7(immunologic gene sets)
## MSigDb over-presentaton analysis
emsig_H <- enricher(gene, TERM2GENE=H_t2g)
emsig_C1 <- enricher(gene, TERM2GENE=C1_t2g)
emsig_C2 <- enricher(gene, TERM2GENE=C2_t2g)
emsig_C3 <- enricher(gene, TERM2GENE=C3_t2g)
emsig_C4 <- enricher(gene, TERM2GENE=C4_t2g)
emsig_C5 <- enricher(gene, TERM2GENE=C5_t2g)
emsig_C6 <- enricher(gene, TERM2GENE=C6_t2g)
emsig_C7 <- enricher(gene, TERM2GENE=C7_t2g)
emsig_H_2 <- setReadable(emsig_H, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C1_2 <- setReadable(emsig_C1, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C2_2 <- setReadable(emsig_C2, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C3_2 <- setReadable(emsig_C3, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C4_2 <- setReadable(emsig_C4, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C5_2 <- setReadable(emsig_C5, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C6_2 <- setReadable(emsig_C6, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
emsig_C7_2 <- setReadable(emsig_C7, OrgDb = db, keyType="ENTREZID")    # ID -> GeneName
################################################################################
## 13. Biological theme comparison(https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html) #####
################################################################################
## Comparing multiple gene lists #####
GENE <- Hmisc::llist(gene_p, gene_c)
## compare KEGG pathway
ck <- compareCluster(geneCluster = GENE, 
                     fun = enrichKEGG,   # "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
                     organism="mmu",
                     # OrgDb = db,
                     # ont = ont,
                     # organism="mmu", 
                     # keyType = "ENTREZID",
                     pvalueCutoff=0.05)
ck <- setReadable(ck, OrgDb = db, keyType="ENTREZID")

## compare reactome pathway
cr <- compareCluster(geneCluster = GENE, 
                     fun = enrichPathway   # "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
                     # organism="mmu",
                     # OrgDb = db,
                     # ont = ont,
                     # organism="mmu", 
                     # keyType = "ENTREZID",
                     # pvalueCutoff=0.05
                     )
################################################################################
## Formula interface of compareCluster
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 0,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "Log2FClessequal1"
mydf$othergroup[abs(mydf$FC) > 1] <- "Log2FCmorethan1"
fr <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG", organism="mmu")
################################################################################
################################################################################
################################################################################
## 99. Visualize enriched pathways #####
################################################################################
################################################################################
################################################################################
## Visualize enriched GO term #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("go")){dir.create("go")}
setwd("./go")

plot_set_go(ego, 
            showCategory = showCategory,             # modify showCategory
            filegp  = "gp_ego.svg",                  # filegp = "gp.svg"
            filebp  = "bp_ego.svg",                  # filebp = "bp.svg"
            filedp  = "dp_ego.svg",                  # filedp = "dp.svg"
            filehp  = "hp_ego.svg",                  # filehp = "hp.svg"
            filetp  = "tp_ego.svg",                  # filetp = "tp.svg"
            filemp  = "mp_ego.svg",                  # filemp = "mp.svg"
            fileusp = "usp_ego.svg",                 # fileusp = "usp.svg"
            filecp  = "cp_ego.svg")                  # filecp = "cp.svg"

plot_set_go(egosim, 
            showCategory = showCategory,             # modify showCategory
            filegp  = "gp_egosim.svg",               # filegp = "gp.svg"
            filebp  = "bp_egosim.svg",               # filebp = "bp.svg"
            filedp  = "dp_egosim.svg",               # filedp = "dp.svg"
            filehp  = "hp_egosim.svg",               # filehp = "hp.svg"
            filetp  = "tp_egosim.svg",               # filetp = "tp.svg"
            filemp  = "mp_egosim.svg",               # filemp = "mp.svg"
            fileusp = "usp_egosim.svg",              # fileusp = "usp.svg"
            filecp  = "cp_egosim.svg")               # filecp = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
## pubmed trend of enriched terms #####
# terms <- ego$Description[1:5]
# p <- enrichplot::pmcplot(terms, 2015:2020)
# p2 <- enrichplot::pmcplot(terms, 2015:2020, proportion=FALSE)
# plot_grid(p, p2, ncol=2)
################################################################################
## Visualize enriched KEGG pathways #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("kegg")){dir.create("kegg")}
setwd("./kegg")

plot_set(ek, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ek.svg",            # filegp  = "gp.svg"
         filebp       = "bp_ek.svg",            # filebp  = "bp.svg"
         filedp       = "dp_ek.svg",            # filedp  = "dp.svg"
         filehp       = "hp_ek.svg",            # filehp  = "hp.svg"
         filetp       = "tp_ek.svg",            # filetp  = "tp.svg"
         filemp       = "mp_ek.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_ek.svg",           # fileusp = "usp.svg"
         filecp       = "cp_ek.svg")            # filecp  = "cp.svg"

plot_set(ek2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ek2.svg",           # filegp  = "gp.svg"
         filebp       = "bp_ek2.svg",           # filebp  = "bp.svg"
         filedp       = "dp_ek2.svg",           # filedp  = "dp.svg"
         filehp       = "hp_ek2.svg",           # filehp  = "hp.svg"
         filetp       = "tp_ek2.svg",           # filetp  = "tp.svg"
         filemp       = "mp_ek2.svg",           # filemp  = "mp.svg"
         fileusp      = "usp_ek2.svg",          # fileusp = "usp.svg"
         filecp       = "cp_ek2.svg")           # filecp  = "cp.svg"

plot_set(ek3,
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ek3.svg",           # filegp  = "gp.svg"
         filebp       = "bp_ek3.svg",           # filebp  = "bp.svg"
         filedp       = "dp_ek3.svg",           # filedp  = "dp.svg"
         filehp       = "hp_ek3.svg",           # filehp  = "hp.svg"
         filetp       = "tp_ek3.svg",           # filetp  = "tp.svg"
         filemp       = "mp_ek3.svg",           # filemp  = "mp.svg"
         fileusp      = "usp_ek3.svg",          # fileusp = "usp.svg"
         filecp       = "cp_ek3.svg")           # filecp  = "cp.svg"

plot_set(mek, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_mek.svg",            # filegp  = "gp.svg"
         filebp       = "bp_mek.svg",            # filebp  = "bp.svg"
         filedp       = "dp_mek.svg",            # filedp  = "dp.svg"
         filehp       = "hp_mek.svg",            # filehp  = "hp.svg"
         filetp       = "tp_mek.svg",            # filetp  = "tp.svg"
         filemp       = "mp_mek.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_mek.svg",           # fileusp = "usp.svg"
         filecp       = "cp_mek.svg")            # filecp  = "cp.svg"

plot_set(mek2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_mek2.svg",           # filegp  = "gp.svg"
         filebp       = "bp_mek2.svg",           # filebp  = "bp.svg"
         filedp       = "dp_mek2.svg",           # filedp  = "dp.svg"
         filehp       = "hp_mek2.svg",           # filehp  = "hp.svg"
         filetp       = "tp_mek2.svg",           # filetp  = "tp.svg"
         filemp       = "mp_mek2.svg",           # filemp  = "mp.svg"
         fileusp      = "usp_mek2.svg",          # fileusp = "usp.svg"
         filecp       = "cp_mek2.svg")           # filecp  = "cp.svg"

plot_set(mek3,
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_mek3.svg",           # filegp  = "gp.svg"
         filebp       = "bp_mek3.svg",           # filebp  = "bp.svg"
         filedp       = "dp_mek3.svg",           # filedp  = "dp.svg"
         filehp       = "hp_mek3.svg",           # filehp  = "hp.svg"
         filetp       = "tp_mek3.svg",           # filetp  = "tp.svg"
         filemp       = "mp_mek3.svg",           # filemp  = "mp.svg"
         fileusp      = "usp_mek3.svg",          # fileusp = "usp.svg"
         filecp       = "cp_mek3.svg")           # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
## KEGG pathway over-representation analysis #####
# browseKEGG(ek,"mmu04721")
# mmu04721 <- pathview(gene.data  = geneList,
#                      pathway.id = "mmu04721",
#                      species    = "mmu",
#                      # limit      = list(gene=max(abs(geneList)), cpd=1)  # limit = list(gene=2, cpd=1)
#                      low = list(gene = "green", cpd = "blue"), 
#                      mid = list(gene = "yellow", cpd = "white"), 
#                      high = list(gene = "red", cpd = "orange"))
################################################################################
## WikiPathways over-representation analysis
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("wiki")){dir.create("wiki")}
setwd("./wiki")

plot_set(ewp2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ewp.svg",            # filegp  = "gp.svg"
         filebp       = "bp_ewp.svg",            # filebp  = "bp.svg"
         filedp       = "dp_ewp.svg",            # filedp  = "dp.svg"
         filehp       = "hp_ewp.svg",            # filehp  = "hp.svg"
         filetp       = "tp_ewp.svg",            # filetp  = "tp.svg"
         filemp       = "mp_ewp.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_ewp.svg",           # fileusp = "usp.svg"
         filecp       = "cp_ewp.svg")            # filecp  = "cp.svg"
## unknown to visuaize the results of wiki pathway analysis
setwd(DIR_HOST)        # move to host working directory
################################################################################
## Visualize enriched reactome pathways
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("reac")){dir.create("reac")}
setwd("./reac")
getwd()

plot_set(er, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_er.svg",            # filegp  = "gp.svg"
         filebp       = "bp_er.svg",            # filebp  = "bp.svg"
         filedp       = "dp_er.svg",            # filedp  = "dp.svg"
         filehp       = "hp_er.svg",            # filehp  = "hp.svg"
         filetp       = "tp_er.svg",            # filetp  = "tp.svg"
         filemp       = "mp_er.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_er.svg",           # fileusp = "usp.svg"
         filecp       = "cp_er.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
## Pathway Visualization #####
# viewPathway("Transmission across Chemical Synapses",
#             organism = "mouse",                          # "human"
#             readable = TRUE, 
#             foldChange = geneList)

## reactome pathway map #####
# pathwayID <- er[1]$ID
# URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png")
# output_filename <- paste0(pathwayID, ".png")
# download.file(URL, output_filename)
## highlighted map #####
# genelist <- str_split(er[1]$geneID, "/")[[1]]
# genelist.str <- str_c(genelist, collapse = ",")
# URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png", "?flg=", genelist.str)
# output_filename <- paste0(pathwayID, "_decorated",".png")
# download.file(URL, output_filename)
################################################################################
## Visualize enriched MeSH terms #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("mesh")){dir.create("mesh")}
setwd("./mesh")
getwd()

plot_set(emesh2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emesh.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emesh.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emesh.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emesh.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emesh.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emesh.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emesh.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emesh.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
################################################################################
## Visualize enriched Cell marker #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("cell")){dir.create("cell")}
setwd("./cell")
getwd()

plot_set(ecm2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ecm.svg",            # filegp  = "gp.svg"
         filebp       = "bp_ecm.svg",            # filebp  = "bp.svg"
         filedp       = "dp_ecm.svg",            # filedp  = "dp.svg"
         filehp       = "hp_ecm.svg",            # filehp  = "hp.svg"
         filetp       = "tp_ecm.svg",            # filetp  = "tp.svg"
         filemp       = "mp_ecm.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_ecm.svg",           # fileusp = "usp.svg"
         filecp       = "cp_ecm.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
################################################################################
## Visualize enriched MSig term #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("msig")){dir.create("msig")}
setwd("./msig")
getwd()

plot_set(emsig_H_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_H.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_H.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_H.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_H.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_H.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_H.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_H.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_H.svg")            # filecp  = "cp.svg"

plot_set(emsig_C1_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C1.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C1.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C1.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C1.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C1.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C1.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C1.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C1.svg")            # filecp  = "cp.svg"

plot_set(emsig_C2_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C2.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C2.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C2.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C2.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C2.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C2.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C2.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C2.svg")            # filecp  = "cp.svg"

plot_set(emsig_C3_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C3.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C3.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C3.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C3.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C3.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C3.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C3.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C3.svg")            # filecp  = "cp.svg"

plot_set(emsig_C4_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C4.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C4.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C4.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C4.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C4.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C4.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C4.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C4.svg")            # filecp  = "cp.svg"

plot_set(emsig_C5_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C5.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C5.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C5.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C5.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C5.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C5.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C5.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C5.svg")            # filecp  = "cp.svg"

plot_set(emsig_C6_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C6.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C6.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C6.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C6.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C6.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C6.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C6.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C6.svg")            # filecp  = "cp.svg"

plot_set(emsig_C7_2, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_emsig_C7.svg",            # filegp  = "gp.svg"
         filebp       = "bp_emsig_C7.svg",            # filebp  = "bp.svg"
         filedp       = "dp_emsig_C7.svg",            # filedp  = "dp.svg"
         filehp       = "hp_emsig_C7.svg",            # filehp  = "hp.svg"
         filetp       = "tp_emsig_C7.svg",            # filetp  = "tp.svg"
         filemp       = "mp_emsig_C7.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_emsig_C7.svg",           # fileusp = "usp.svg"
         filecp       = "cp_emsig_C7.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
################################################################################
## Visualize enriched compare gene list #####
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("cgl")){dir.create("cgl")}
setwd("./cgl")
getwd()

plot_set(ck, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_ck.svg",            # filegp  = "gp.svg"
         filebp       = "bp_ck.svg",            # filebp  = "bp.svg"
         filedp       = "dp_ck.svg",            # filedp  = "dp.svg"
         filehp       = "hp_ck.svg",            # filehp  = "hp.svg"
         filetp       = "tp_ck.svg",            # filetp  = "tp.svg"
         filemp       = "mp_ck.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_ck.svg",           # fileusp = "usp.svg"
         filecp       = "cp_ck.svg")            # filecp  = "cp.svg"

plot_set(cr, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_cr.svg",            # filegp  = "gp.svg"
         filebp       = "bp_cr.svg",            # filebp  = "bp.svg"
         filedp       = "dp_cr.svg",            # filedp  = "dp.svg"
         filehp       = "hp_cr.svg",            # filehp  = "hp.svg"
         filetp       = "tp_cr.svg",            # filetp  = "tp.svg"
         filemp       = "mp_cr.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_cr.svg",           # fileusp = "usp.svg"
         filecp       = "cp_cr.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
################################################################################
## Visualization of functional profile comparison ###
## file save directory #####
## make directory to save images of enrichment analysis
setwd(DIR_HOST)        # move to host working directory
getwd()
if(!dir.exists("p")){dir.create("p")}
setwd("./p")
if(!dir.exists("fr")){dir.create("fr")}
setwd("./fr")
getwd()

plot_set(fr, 
         showCategory = showCategory,             # modify showCategory
         filegp       = "gp_fr.svg",              # filegp  = "gp.svg"
         filebp       = "bp_fr.svg",              # filebp  = "bp.svg"
         filedp       = "dp_fr.svg",            # filedp  = "dp.svg"
         filedp2      = "dp2_fr.svg",           # filedp2  = "dp2.svg"
         filehp       = "hp_fr.svg",            # filehp  = "hp.svg"
         filetp       = "tp_fr.svg",            # filetp  = "tp.svg"
         filemp       = "mp_fr.svg",            # filemp  = "mp.svg"
         fileusp      = "usp_fr.svg",           # fileusp = "usp.svg"
         filecp       = "cp_fr.svg")            # filecp  = "cp.svg"

setwd(DIR_HOST)        # move to host working directory
################################################################################
# ここから






################################################################################
## 999. gene set enrichment analysis(GSEA) # for more big data set(e.g., microarray, sequencing) #####
################################################################################
##  GO analysis using other gene ID #####
# gene.df <- bitr(gene, fromType = "ENTREZID",      # Gene ID -> ENSEMBL, SYMBOL
#                 toType = c("ENSEMBL", "SYMBOL"),  # add annotation
#                 OrgDb = org.Mm.eg.db              # org.Hs.eg.db
#                 )
# ego2 <- enrichGO(gene          = gene.df$ENSEMBL,
#                  OrgDb         = org.Mm.eg.db,    # org.Hs.eg.db
#                  keyType       = 'ENSEMBL',       # analysed by ENSEMBLE
#                  ont           = "BP", 　　　　   # "BP","CC","MF","ALL"
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.05,
#                  qvalueCutoff  = 0.05,
#                  readable      = TRUE 　　　　　   # Gene ID -> gene name
#                  )
# head(ego2, 3)

## GO GSEA  # error #####
# ego3 <- gseGO(geneList     = geneList,   # log2FC
#               OrgDb        = db,         # org.Hs.eg.db
#               ont          = ont, 　　   # "BP","CC","MF","ALL"
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               qvalueCutoff = 1,
#               pAdjustMethod = "BH",
#               verbose      = FALSE)
# head(as.data.frame(ego3))
# ridgeplot(ego3, showCategory = showCategory)
# gseaplot(ego3, geneSetID = "GO:0048598")

## KEGG pathway GSEA  # error #####
# ek2 <- gseKEGG(geneList       = geneList,
#                organism       = 'mmu',       # 'hsa'
#                # minGSSize      = 120,
#                pvalueCutoff   = 0.5,
#                verbose        = FALSE)
# head(ek2)

## KEGG module GSEA  # error #####
# mek2 <- gseMKEGG(geneList = geneList,
#                  organism = 'mmu',
#                  pvalueCutoff = 1)
# head(mkk2)

## WikiPathways GSEA  # error #####
# ewp2 <- gseWP(geneList, organism = "Mus musculus")

## Reactome pathway gene set enrichment analysis #####
# er2 <- gsePathway(geneList, 
#                   pvalueCutoff = 0.2,
#                   pAdjustMethod = "BH", 
#                   verbose = FALSE)
# head(er2)

## human only #####
## Disease enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html) #####
## Disease over-representation analysis
# edo <- enrichDO(gene          = gene,
#                 ont           = "DO",
#                 pvalueCutoff  = 0.05,
#                 pAdjustMethod = "BH",
#                 universe      = names(geneList),
#                 minGSSize     = 5,
#                 maxGSSize     = 500,
#                 qvalueCutoff  = 0.05,
#                 readable      = FALSE)
# head(edo)

## Over-representation analysis for the network of cancer gene #####
# ncg <- enrichNCG(gene) 
# head(ncg)

## Over-representation analysis for the disease gene network #####
# dgn <- enrichDGN(gene) 
# head(dgn)

# snp <- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392",
#          "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396",
#          "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646",
#          "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221")
# dgnv <- enrichDGNv(snp)
# head(dgnv)

## Disease GSEA #####
# edo2 <- gseDO(geneList,
#               minGSSize     = 120,
#               pvalueCutoff  = 0.2,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# head(edo2, 3)

## Disease GSEA for the network of cancer gene #####
# ncg <- gseNCG(geneList,
#               pvalueCutoff  = 0.5,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# ncg <- setReadable(ncg, 'org.Hs.eg.db')
# head(ncg, 3) 

## gseDGN fuction #####
# dgn <- gseDGN(geneList,
#               pvalueCutoff  = 0.2,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# dgn <- setReadable(dgn, 'org.Hs.eg.db')    # ID conversion
# head(dgn, 3) 

## MeSH GSEA #####
# emesh2 <- gseMeSH(geneList, MeSHDb = MeSHDb, database = 'gene2pubmed', category = "G")

## Cell Marker GSEA #####
# ecm2 <- GSEA(geneList, TERM2GENE = cells)
# head(ecm2)

## MSigDb GSEA #####
# emsig2_C3 <- GSEA(geneList, TERM2GENE = C3_t2g)
# head(emsig2_C3)

## visualization of GSEA(https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) #####
## ridgeline plot for expression distribution of GSEA result
# ridgeplot(edo2)
# enrichplot::gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])
# gsearank(edo2, 1, title = edo2[1, "Description"])
# pp <- lapply(1:3, function(i) {
#   anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
#   lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
#   
#   gsearank(edo2, i, edo2[i, 2]) + xlab(NULL) +ylab(NULL) +
#     annotate("text", 10000, edo2[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
# })
# plot_grid(plotlist=pp, ncol=1)
################################################################################
## 9999. note #####
################################################################################
