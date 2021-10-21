# enrichment analysis
# clusterProfiler(https://yulab-smu.top/biomedical-knowledge-mining-book/index.html)
# Repidemiology(https://jojoshin.hatenablog.com/entry/2016/09/27/161744)
# GO, KEGG(https://shiokoji11235.com/go-analysis-for-transcriptome-analysis)
# GSEA(https://www.ncc.go.jp/jp/ri/division/rare_cancer_research/labo/20181011150525.html)
# GSEA(https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#abstract)
# GSEA(https://www.biostars.org/p/445059/)
# reactome pathway(https://shiokoji11235.com/pathway-enrichment-analysis-in-r)
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
source("/home/rstudio/rproject/script/archive/functions.R")        # load functions
source("/home/rstudio/rproject/script/archive/functions_DEPpkg.R") # load DEPpkg functions
detach_all()
source("/home/rstudio/rproject/script/startup.R")                  # load packages
################################################################################
getwd()
# setwd("/home/rstudio/rproject")
# setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")
setwd("data/Perseus_Like_Analysis/PFC")
dat <- read_excel("stat.xlsx", 1)

# geneList Ref(https://www.biostars.org/p/192968/)
t(colnames(dat))
dat_rm <- dat %>% drop_na(ENTREZID)

dat_ps <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05) %>%
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()

dat_ud <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 > 0) %>%
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 < 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 < 0) %>% 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()

dat_du <- dat %>%
  drop_na(ENTREZID) %>% 
  filter(q_PxC < 0.05) %>% 
  filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 < 0) %>%
  filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 > 0 | THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 > 0) %>% 
  dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
  dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
  t()

# dat_du2 <- dat %>%
#   drop_na(ENTREZID) %>% 
#   filter(q_PxC < 0.05) %>% 
#   filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 < 0) %>%
#   # filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 > 0) %>%
#   filter(THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 > 0) %>%
#   dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
#   dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
#   t()

# dat_du3 <- dat %>%
#   drop_na(ENTREZID) %>% 
#   filter(q_PxC < 0.05) %>% 
#   filter(THSDp_P.C0.S.C0 < 0.05 & diff_P.C0.S.C0 < 0) %>%
#   # filter(THSDp_P.C10.P.C0 < 0.05 & diff_P.C10.P.C0 > 0) %>%
#   filter(THSDp_P.C30.P.C0 < 0.05 & diff_P.C30.P.C0 > 0) %>%
#   dplyr::select(ENTREZID, diff_P.C0.S.C0) %>%
#   dplyr::arrange(desc(diff_P.C0.S.C0)) %>% 
#   t()

# intersect(dat_du2, dat_du3)

# define object as "dat_gse"(substitution)
dat_gse <- dat_ps
dat_gse <- dat_ud
dat_gse <- dat_du


gene_ls <- dat_gse[-1,]                        # log2FC
gene_ls <- as.numeric(gene_ls)                 # chr -> num
names(gene_ls) <- dat_gse[1,]                  # colnames
str(gene_ls)
names(gene_ls)
dat_rm$ENTREZID

# write.table(gene_ls, file = "gse.txt", row.names = TRUE, col.names = TRUE)
# df <- read.table(file = "gse.txt", header = TRUE, dec = ".")
# colnames(df) <- "log2FC"
################################################################################
# Ref(https://shiokoji11235.com/go-analysis-for-transcriptome-analysis)
################################################################################
# load geneList from DOSE pkg
# data(geneList, package = "DOSE")
# gene <- names(geneList)[abs(geneList) > 2]

#run GO enrichment analysis
gene <- names(gene_ls)[abs(gene_ls) > 0]
ego_result <- enrichGO(gene          = gene,                         # Entrez Gene ID
                       # universe      = names(gene_ls),
                       # universe      = dat_rm$ENTREZID,              # Entrez Gene ID detected in analysis
                       # OrgDb         = org.Hs.eg.db,
                       OrgDb         = org.Mm.eg.db,                 # mouse
                       ont           = "BP", 　　　　　              # "BP","CC","MF","ALL"
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE 　　　　　              # Gene ID -> gene name
                       )

head(as.data.frame(ego_result))
ego_result.simple <- clusterProfiler::simplify(ego_result)           # remove redundancy
head(as.data.frame(ego_result.simple))
################################################################################
# bar plot
barplot(ego_result, drop=TRUE, showCategory=20)
barplot(ego_result.simple, drop=TRUE, showCategory=20)

# dot plot
clusterProfiler::dotplot(ego_result, showCategory=20)
clusterProfiler::dotplot(ego_result.simple, showCategory=20)
# heat plot
heatplot(ego_result, foldChange = gene_ls, showCategory = 20) + ggplot2::coord_flip()
heatplot(ego_result.simple, foldChange = gene_ls, showCategory = 20) + ggplot2::coord_flip()
# enrichment map plot
# trouble shooting: Ref(https://github.com/YuLab-SMU/enrichplot/issues/79)
# d <- GOSemSim::godata("org.Hs.eg.db", ont = "CC")
d <- GOSemSim::godata("org.Mm.eg.db", ont = "BP")    
compare_cluster_GO_emap <- enrichplot::pairwise_termsim(ego_result, semData = d,  method="Wang")
compare_cluster_GO_emap2 <- enrichplot::pairwise_termsim(ego_result.simple, semData = d,  method="Wang")
clusterProfiler::emapplot(compare_cluster_GO_emap,     # enrichment result
                          showCategory = 30,           # number of term
                          color = "p.adjust",          # pvalue, p.adjust, qvalue
                          layout = "nicely",           # "kk", "nicely"
                          node_scale = NULL,           # scale of node, changed to cex_category
                          line_scale = NULL,           # scale of line width, changed to cex_line
                          min_edge = 0.2,              # minimum percentage of overlap genes 0-1(0.2)
                          node_label_size = NULL,      # size of node label,changed to cex_label_category
                          cex_label_category = 0.7,    # scale of category node(pathway) label size
                          # cex_category = 1,          # number indicating the amount by which plotting category nodes should be scaled relative to the default, NULL=error, 1 =<
                          # cex_line = 1,              # scale of line width, NULL=error
                          # split = 100,                 # separate result by 'category' variable 
                          pie = 'equal',               # proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
                          # legend_n = 100,            # number of circle in legend
                          pie_scale = NULL             # scale of pie chart or point, this parameter has been changed to "node_scale"
                          )
clusterProfiler::emapplot(compare_cluster_GO_emap2)

# enrichMap(ego_result.simple)                            # error

#cnet plot
clusterProfiler::cnetplot(ego_result, 
                          categorySize = "pvalue", 
                          # foldChange=geneList,
                          foldChange = gene_ls
                          )
clusterProfiler::cnetplot(ego_result.simple, 
                          categorySize = "pvalue", 
                          foldChange = gene_ls
                          )

#GO DAG graph
goplot(ego_result)
goplot(ego_result.simple)
################################################################################
# gene set enrichment analysis
gse_result <- gseGO(# geneList     = geneList,
                    geneList     = gene_ls,
                    # OrgDb        = org.Hs.eg.db,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "ALL",　                　#"BP","CC","MF","ALL"
                    # nPerm        = 1000,
                    minGSSize    = 3,
                    maxGSSize    = 800, 
                    pvalueCutoff = 0.5,
                    # qvalueCutoff = 1,
                    verbose      = TRUE,
                    pAdjustMethod = "BH")
head(as.data.frame(gse_result))


ridgeplot(gse_result,showCategory = 8)
gseaplot(gse_result,geneSetID = "GO:0048598")


mkk <- enrichMKEGG(gene = gene,
                   organism = 'mmu',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1
                   )
head(mkk)
barplot(mkk, drop=FALSE, showCategory=12)
dotplot(mkk)
heatplot(mkk, foldChange = gene_ls, showCategory = 20000) + ggplot2::coord_flip()
# DOSE::enrichMap(mkk)  # error
cnetplot(mkk, categorySize="pvalue", foldChange=gene_ls)
cnetplot(mkk, categorySize="geneNum", foldChange=gene_ls)
browseKEGG(mkk,"M00143")
pathview(gene.data = gene_ls, pathway.id = "M00143", limit = list(gene=2, cpd=1))

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
