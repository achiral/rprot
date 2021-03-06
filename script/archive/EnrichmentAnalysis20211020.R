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
require(ggnewscale)           # enrichment map plot
require(ggplot2)
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
setwd("data/Perseus_Like_Analysis/PFC")
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
## 5. GO enrichment analysis #####
################################################################################
ont <- "BP"                                     # "BP"(default ),"CC","MF","ALL"
showCategory <- 20                              # default = 20
layout = "kk"                                   # "kk", "nicely"
labeln = "all"                                  # "category", "gene", "all", "none"
db <- org.Mm.eg.db
d <- GOSemSim::godata(OrgDb = db, ont = ont)    # emap plot
MeSHDb <- "MeSH.Mmu.eg.db"                      # MeSHDb <- MeSH.Hsa.eg.db    # BioC 2.14-3.13

geneList <- prepare_geneList(dat_p)
gene <- names(geneList)[abs(geneList) > 0]

ego_p <- run_enrichGO(dat_p, db = db,
                      ont = ont)                # modify GO term: "BP","CC","MF","ALL"
ego_psim <- clusterProfiler::simplify(ego_p)    # remove redundancy
plot_set(ego_p, 
         showCategory = showCategory,           # modify showCategory
         filegp = "gp_ego_p.svg",               # filegp = "gp.svg"
         filebp = "bp_ego_p.svg",               # filebp = "bp.svg"
         filedp = "dp_ego_p.svg",               # filedp = "dp.svg"
         filehp = "hp_ego_p.svg",               # filehp = "hp.svg"
         filemp = "mp_ego_p.svg",               # filemp = "mp.svg"
         filecp = "cp_ego_p.svg")               # filecp = "cp.svg"
plot_set(ego_psim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_psim.svg",            # filegp = "gp.svg"
         filebp = "bp_ego_psim.svg",            # filebp = "bp.svg"
         filedp = "dp_ego_psim.svg",            # filedp = "dp.svg"
         filehp = "hp_ego_psim.svg",            # filehp = "hp.svg"
         filemp = "mp_ego_psim.svg",            # filemp = "mp.svg"
         filecp = "cp_ego_psim.svg")            # filecp = "cp.svg"

## 保留 #####
geneList <- prepare_geneList(dat_c)
ego_c <- run_enrichGO(dat_c, 
                      ont = ont)                # modify GO term: "BP","CC","MF","ALL"
ego_csim <- clusterProfiler::simplify(ego_c)    # remove redundancy
plot_set(ego_c, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_c.svg",               # filegp = "gp.svg"
         filebp = "bp_ego_c.svg",               # filebp = "bp.svg"
         filedp = "dp_ego_c.svg",               # filedp = "dp.svg"
         filehp = "hp_ego_c.svg")               # filehp = "hp.svg"
plot_set(ego_csim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_csim.svg",            # filegp = "gp.svg"
         filebp = "bp_ego_csim.svg",            # filebp = "bp.svg"
         filedp = "dp_ego_csim.svg",            # filedp = "dp.svg"
         filehp = "hp_ego_cs.svg")              # filehp = "hp.svg"   

geneList <- prepare_geneList(dat_pc)
ego_pc <- run_enrichGO(dat_pc, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_pcsim <- clusterProfiler::simplify(ego_pc)  # remove redundancy
plot_set(ego_pc, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_pc.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_pc.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_pc.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_pc.svg")              # filehp = "hp.svg"
plot_set(ego_pcsim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_pcsim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_pcsim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_pcsim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_pcsim.svg")           # filehp = "hp.svg"

geneList <- prepare_geneList(dat_ps)
ego_ps <- run_enrichGO(dat_ps, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_pssim <- clusterProfiler::simplify(ego_ps)  # remove redundancy
plot_set(ego_ps, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_ps.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_ps.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_ps.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_ps.svg")              # filehp = "hp.svg"
plot_set(ego_pssim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_pssim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_pssim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_pssim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_pssim.svg")           # filehp = "hp.svg"

geneList <- prepare_geneList(dat_ud)
ego_ud <- run_enrichGO(dat_ud, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_udsim <- clusterProfiler::simplify(ego_ud)  # remove redundancy
plot_set(ego_ud, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_ud.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_ud.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_ud.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_ud.svg")              # filehp = "hp.svg"
plot_set(ego_udsim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_udsim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_udsim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_udsim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_udsim.svg")           # filehp = "hp.svg"

geneList <- prepare_geneList(dat_du)
ego_du <- run_enrichGO(dat_du, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_dusim <- clusterProfiler::simplify(ego_du)  # remove redundancy
plot_set(ego_du, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_du.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_du.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_du.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_du.svg")              # filehp = "hp.svg"
plot_set(ego_dusim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_dusim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_dusim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_dusim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_dusim.svg")           # filehp = "hp.svg"

geneList <- prepare_geneList(dat_uu)
ego_uu <- run_enrichGO(dat_uu, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_uusim <- clusterProfiler::simplify(ego_uu)  # remove redundancy
plot_set(ego_uu, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_uu.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_uu.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_uu.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_uu.svg")              # filehp = "hp.svg"
plot_set(ego_uusim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_uusim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_uusim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_uusim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_uusim.svg")           # filehp = "hp.svg"

geneList <- prepare_geneList(dat_dd)
ego_dd <- run_enrichGO(dat_dd, 
                       ont = ont)               # modify GO term: "BP","CC","MF","ALL"
ego_ddsim <- clusterProfiler::simplify(ego_dd)  # remove redundancy
plot_set(ego_dd, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_dd.svg",              # filegp = "gp.svg"
         filebp = "bp_ego_dd.svg",              # filebp = "bp.svg"
         filedp = "dp_ego_dd.svg",              # filedp = "dp.svg"
         filehp = "hp_ego_dd.svg")              # filehp = "hp.svg"
plot_set(ego_ddsim, 
         showCategory = showCategory,           # modify showCategory 
         filegp = "gp_ego_ddsim.svg",           # filegp = "gp.svg"
         filebp = "bp_ego_ddsim.svg",           # filebp = "bp.svg"
         filedp = "dp_ego_ddsim.svg",           # filedp = "dp.svg"
         filehp = "hp_ego_ddsim.svg")           # filehp = "hp.svg"
################################################################################
## 99. visualize enriched GO terms function modification #####
################################################################################
## output plot set as svg #####
plot_set <- function(e, showCategory = 20, layout = "nicely", filegp = "gp.svg", filebp = "bp.svg", filedp = "dp.svg", filehp = "hp.svg", filemp = "mp.svg", filecp = "cp.svg", 
                     labeln = "all", # "category", "gene", "all", "none"
                     label = 0.7, labelg = 0.5, pie = "count"){
  ## goplot
  gp <- goplot(e)
  ## bar plot
  bp <- barplot(e, drop=TRUE, showCategory=showCategory)
  ## dot plot
  dp <- clusterProfiler::dotplot(e, showCategory=showCategory)
  ## heat plot
  hp <- heatplot(e, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()
  ## emapplot
  # d <- GOSemSim::godata(OrgDb = db, ont = ont)
  em <- enrichplot::pairwise_termsim(e, semData = d,  method="Wang")
  mp <- clusterProfiler::emapplot(em, showCategory = showCategory, layout = layout, cex_label_category = label, pie = pie)
  ## cnet plot
  cp <- clusterProfiler::cnetplot(e, showCategory = showCategory, foldChange = geneList, layout = layout, colorEdge = TRUE, node_label = labeln, cex_label_category = label, cex_label_gene = labelg, color_category='firebrick', categorySize = "pvalue")
  ## output svg
  dev_svg(plot = gp, file = filegp)
  dev_svg(plot = bp, file = filebp)
  dev_svg(plot = dp, file = filedp)
  dev_svg(plot = hp, file = filehp)
  dev_svg(plot = mp, file = filemp)
  dev_svg(plot = cp, file = filecp)
}
################################################################################
## 6. KEGG pathway analysis #####
################################################################################
# KEGG pathway over-representation analysis #####
ek <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',       # 'hsa'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
head(ek)
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
head(summary(ek3))
## KEGG module over-representation analysis
mek <- enrichMKEGG(gene = gene,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mek)     
mek2 <- enrichMKEGG(gene = eg2up[,2],   # gene = eg2np[,2]
                    organism = 'mmu',
                    keyType='uniprot',  # keyType='ncbi-proteinid'
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
## ID conversion
mek3 <- setReadable(mek2, OrgDb = db, keyType="UNIPROT")
head(summary(mek3))
################################################################################
## 99. Visualize enriched KEGG pathways #####
################################################################################
## Visualize enriched KEGG pathways
# KEGG pathway over-representation analysis
browseKEGG(ek,"mmu04721")
mmu04721 <- pathview(gene.data  = geneList,
                     pathway.id = "mmu04721",
                     species    = "mmu",
                     # limit      = list(gene=max(abs(geneList)), cpd=1)  # limit = list(gene=2, cpd=1)
                     low = list(gene = "green", cpd = "blue"), 
                     mid = list(gene = "yellow", cpd = "white"), 
                     high = list(gene = "red", cpd = "orange"))
barplot(ek, drop=FALSE, showCategory = showCategory)
dotplot(ek)
heatplot(ek, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()  # DENTREZID, color
heatplot(ek2, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # uniprot, gray
heatplot(ek3, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # gene name, gray
cnetplot(ek, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", DENTREZID, color
cnetplot(ek3, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", gene name, gray
# KEGG module over-representation analysis
barplot(mek, drop=FALSE, showCategory = showCategory)
dotplot(mek)
heatplot(mek, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()  # DENTREZID, color
heatplot(mek2, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # uniprot, gray
heatplot(mek3, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()    # gene name, gray
cnetplot(mek, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", DENTREZID, color
cnetplot(mek3, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", gene name, gray
################################################################################
## 7. WikiPathways analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html) #####
################################################################################
get_wp_organisms()
## WikiPathways over-representation analysis
ewp <- enrichWP(gene, organism = "Mus musculus")
head(ewp@result)
## visualization???
# ここ保留
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
                    readable = TRUE                # gene ID -> gene Name
                    )
head(er)

## Pathway Visualization
viewPathway("Transmission across Chemical Synapses",
            organism = "mouse",                          # "human"
            readable = TRUE, 
            foldChange = geneList)
## visualization
barplot(er, showCategory = showCategory, x = "Count")
dotplot(er, showCategory = showCategory)
heatplot(er, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()
cnetplot(er, foldChange = geneList, showCategory = showCategory, categorySize = "pvalue")
## reactome pathway map
pathwayID <- er[1]$ID
URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png")
output_filename <- paste0(pathwayID, ".png")
download.file(URL, output_filename)
## highlighted map
genelist <- str_split(er[1]$geneID, "/")[[1]]
genelist.str <- str_c(genelist, collapse = ",")
URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png", "?flg=", genelist.str)
output_filename <- paste0(pathwayID, "_decorated",".png")
download.file(URL, output_filename)
################################################################################
## 9. Disease enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html) #####
################################################################################
# only for human(see bottom)
################################################################################
## 10. MeSH enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/meshes-enrichment.html) #####
################################################################################
## MeSH over-representation analysis
de <- names(geneList)[1:100]
emesh <- enrichMeSH(de, MeSHDb = MeSHDb, 
                    database = 'gendoo',        # 'gendoo', 'gene2pubmed' or 'RBBH'
                    category = 'C',           # "A", "B", "C"(Diseases), "D", "E", "F", "G"(Phenomena and Processes), "H", "I", "J", "K", "L","M", "N", "V", "Z"
                    # universe,
                    # minGSSize = 10, maxGSSize = 500,
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH",     # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                    qvalueCutoff = 1)
head(emesh)
################################################################################
## 11. Universal enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html) #####
################################################################################
## 11. Cell Marker Data #####
################################################################################
cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt')
# cells <- cell_marker_data %>%                          # error
#   dplyr::select(cellName, geneID) %>%
#   dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
#   tidyr::unnest()
cmd2 <- dplyr::select(cell_marker_data, cellName, geneID)
cmd3 <- dplyr::mutate(cmd2, geneID = strsplit(geneID, ', '))
cells <- tidyr::unnest(cmd3)
## Cell Marker over-presentaton analysis
ecm <- enricher(gene, TERM2GENE = cells)
head(ecm)
################################################################################
## 12. MSigDb analysis #####
################################################################################
## check species
msigdbr_species()     # msigdbr_show_species()
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame

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
head(emsig_H)
head(emsig_C1)
head(emsig_C2)
head(emsig_C3)
head(emsig_C4)
head(emsig_C5)
head(emsig_C6)
head(emsig_C7)
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
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck)
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
head(cr)
## Visualization of functional profile comparison
# dotplot(ck)  # error
## Gene-Concept Network
cnetplot(ck)
cnetplot(cr)

# head(ck)
## Visualization of functional profile comparison
# dotplot(ck)
## Gene-Concept Network
# cnetplot(ck)




## Formula interface of compareCluster
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 0,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"
formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG", organism="mmu")
head(formula_res)

## Visualization of functional profile comparison
dotplot(formula_res)
dotplot(formula_res, x="group") + facet_grid(~othergroup)



# ここから

################################################################################
## 999. gene set enrichment analysis(GSEA) # for more big data set(e.g., microarray, sequencing)
################################################################################
##  GO analysis using other gene ID
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

## GO GSEA  # error
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

## KEGG pathway GSEA  # error
# ek2 <- gseKEGG(geneList       = geneList,
#                organism       = 'mmu',       # 'hsa'
#                # minGSSize      = 120,
#                pvalueCutoff   = 0.5,
#                verbose        = FALSE)
# head(ek2)

## KEGG module GSEA  # error
# mek2 <- gseMKEGG(geneList = geneList,
#                  organism = 'mmu',
#                  pvalueCutoff = 1)
# head(mkk2)

## WikiPathways GSEA  # error
# ewp2 <- gseWP(geneList, organism = "Mus musculus")

## Reactome pathway gene set enrichment analysis
# er2 <- gsePathway(geneList, 
#                   pvalueCutoff = 0.2,
#                   pAdjustMethod = "BH", 
#                   verbose = FALSE)
# head(er2)

## human only
## Disease enrichment analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html)
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

## Over-representation analysis for the network of cancer gene
# ncg <- enrichNCG(gene) 
# head(ncg)

## Over-representation analysis for the disease gene network
# dgn <- enrichDGN(gene) 
# head(dgn)

# snp <- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392",
#          "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396",
#          "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646",
#          "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221")
# dgnv <- enrichDGNv(snp)
# head(dgnv)

## Disease GSEA
# edo2 <- gseDO(geneList,
#               minGSSize     = 120,
#               pvalueCutoff  = 0.2,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# head(edo2, 3)

## Disease GSEA for the network of cancer gene
# ncg <- gseNCG(geneList,
#               pvalueCutoff  = 0.5,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# ncg <- setReadable(ncg, 'org.Hs.eg.db')
# head(ncg, 3) 

## gseDGN fuction
# dgn <- gseDGN(geneList,
#               pvalueCutoff  = 0.2,
#               pAdjustMethod = "BH",
#               verbose       = FALSE)
# dgn <- setReadable(dgn, 'org.Hs.eg.db')    # ID conversion
# head(dgn, 3) 

## MeSH GSEA
# emesh2 <- gseMeSH(geneList, MeSHDb = MeSHDb, database = 'gene2pubmed', category = "G")

## Cell Marker GSEA
# ecm2 <- GSEA(geneList, TERM2GENE = cells)
# head(ecm2)

## MSigDb GSEA
# emsig2_C3 <- GSEA(geneList, TERM2GENE = C3_t2g)
# head(emsig2_C3)
################################################################################
## 9999. note #####
################################################################################
ego_p
ego_psim
ego_c
ego_csim
ego_pc
ego_pcsim
ego_ps
ego_pssim
ego_ud
ego_udsim
ego_du
ego_uusim
ego_dd
ego_ddsim