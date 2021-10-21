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
# BiocManager::install("clusterProfiler")
# install.packages("ggnewscale")
require(clusterProfiler)
require(DOSE)
require(ggnewscale)           # enrichment map plot
require(ggplot2)
# require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(pathview)
require(ReactomePA)
require(stringr)
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
## set working directory
getwd()
# setwd("/home/rstudio/rproject")
# setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")
# setwd("~/Dropbox/GitHub/local/Docker/R/rprot")
setwd("data/Perseus_Like_Analysis/PFC")
## input data set
dat <- read_excel("stat.xlsx", 1)
## prepare geneList
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
################################################################################
## run GO enrichment analysis
ont <- "BP"                                     # "BP"(default ),"CC","MF","ALL"
showCategory <- 20                              # default = 20
layout = "kk"                                   # "kk", "nicely"
labeln = "all"                                  # "category", "gene", "all", "none"
db <- org.Mm.eg.db
d <- GOSemSim::godata(OrgDb = db, ont = ont)    # emap plot

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
         filehp = "hp_ego_psim.svg")            # filehp = "hp.svg"   

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
## visualize enriched GO terms
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





## output plot set as svg
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
# KEGG pathway over-representation analysis
kk <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',       # 'hsa'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
head(kk)
## ID conversion(http://guangchuangyu.github.io/2016/05/convert-biological-id-with-kegg-api-using-clusterprofiler/)
eg2np <- bitr_kegg(gene, fromType='kegg', toType='ncbi-proteinid', organism='mmu') # Path, Module, ncbi-proteinid, ncbi-geneid, uniprot, kegg
eg2up <- bitr_kegg(eg2np[,2], fromType='ncbi-proteinid', toType='uniprot', organism='mmu') # Path, Module, ncbi-proteinid, ncbi-geneid, uniprot, kegg
## KEGG pathway enrichment analysis
kk2 <- enrichKEGG(gene = eg2up[,2],   # gene = eg2np[,2]
                  organism = 'mmu',
                  keyType='uniprot',  # keyType='ncbi-proteinid'
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1)
## ID conversion
kk3 <- setReadable(kk2, OrgDb = db, keyType="UNIPROT")
head(summary(kk3))
## KEGG module over-representation analysis
mkk <- enrichMKEGG(gene = gene,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)     
mkk2 <- enrichMKEGG(gene = eg2up[,2],   # gene = eg2np[,2]
                    organism = 'mmu',
                    keyType='uniprot',  # keyType='ncbi-proteinid'
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
## ID conversion
mkk3 <- setReadable(mkk2, OrgDb = db, keyType="UNIPROT")
head(summary(mkk3))
################################################################################
## Visualize enriched KEGG pathways
# KEGG pathway over-representation analysis
browseKEGG(kk,"mmu04721")
mmu04721 <- pathview(gene.data  = geneList,
                     pathway.id = "mmu04721",
                     species    = "mmu",
                     # limit      = list(gene=max(abs(geneList)), cpd=1)  # limit = list(gene=2, cpd=1)
                     low = list(gene = "green", cpd = "blue"), 
                     mid = list(gene = "yellow", cpd = "white"), 
                     high = list(gene = "red", cpd = "orange"))
barplot(kk, drop=FALSE, showCategory = showCategory)
dotplot(kk)
heatplot(kk, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()  # DENTREZID, color
heatplot(kk2, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # uniprot, gray
heatplot(kk3, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # gene name, gray
cnetplot(kk, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", DENTREZID, color
cnetplot(kk3, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", gene name, gray
# KEGG module over-representation analysis
barplot(mkk, drop=FALSE, showCategory = showCategory)
dotplot(mkk)
heatplot(mkk, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()  # DENTREZID, color
heatplot(mkk2, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip() # uniprot, gray
heatplot(mkk3, foldChange = geneList, showCategory = showCategory) + ggplot2::coord_flip()    # gene name, gray
cnetplot(mkk, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", DENTREZID, color
cnetplot(mkk3, categorySize="pvalue", foldChange=geneList) # categorySize="geneNum", gene name, gray
################################################################################
# 8 WikiPathways analysis(https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html)
################################################################################

# ここから

################################################################################
# reactome pathway
################################################################################
Reactome_enrichment_result <- enrichPathway(gene = gene,                 # entrez gene id
                                            organism = "mouse",          # "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
                                            pvalueCutoff = 0.05,
                                            # pAdjustMethod = "BH",        # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                            # qvalueCutoff = 0.05,
                                            # universe = dat_rm$ENTREZID,  # background genes
                                            # minGSSize = 3,               # 10
                                            # maxGSSize = 500,
                                            readable = TRUE              # gene ID -> gene Name
                                            )
head(as.data.frame(Reactome_enrichment_result))
# head(summary(Reactome_enrichment_result))
# plot(Reactome_enrichment_result)               # error
barplot(Reactome_enrichment_result, showCategory=100, x = "Count")
dotplot(Reactome_enrichment_result, showCategory=100)
heatplot(Reactome_enrichment_result, foldChange = gene_ls, showCategory = 100) + ggplot2::coord_flip()
# emapplot(Reactome_enrichment_result)  # error
cnetplot(Reactome_enrichment_result, 
         categorySize = "pvalue", 
         foldChange = gene_ls,
         showCategory = 100)

pathwayID <- Reactome_enrichment_result[1]$ID
URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png")
output_filename <- paste0(pathwayID, ".png")
download.file(URL, output_filename)

genelist <- str_split(Reactome_enrichment_result[1]$geneID, "/")[[1]]
genelist.str <- str_c(genelist, collapse = ",")
URL <- paste0("https://reactome.org/ContentService/exporter/diagram/", pathwayID, ".png", "?flg=", genelist.str)
output_filename <- paste0(pathwayID, "_decorated",".png")
download.file(URL, output_filename)


################################################################################
## run GO enrichment analysis
################################################################################
##  other gene ID
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
# kk2 <- gseKEGG(geneList       = geneList,
#                organism       = 'mmu',       # 'hsa'
#                # minGSSize      = 120,
#                pvalueCutoff   = 0.5,
#                verbose        = FALSE)
# head(kk2)

## KEGG module GSEA  # error
# mkk2 <- gseMKEGG(geneList = geneList,
#                  organism = 'mmu',
#                  pvalueCutoff = 1)
# head(mkk2)
################################################################################