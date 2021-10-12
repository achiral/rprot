# Docker base image
# FROM bioconductor/bioconductor_docker
# FROM bioconductor/bioconductor_docker:RELEASE_3_13
FROM achiral/rbioc:latest

# information
# LABEL version="2021.10"
# LABEL description="Docker RStudio bioconductor"
# LABEL maintainer="https://hub.docker.com/u/bioconductor"

# # OS env JP locale (Character, DateTime)
# ENV LANG ja_JP.UTF-8\
#   && LC_ALL ja_JP.UTF-8\
#   && R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
  
# # change directory
# WORKDIR /home/rstudio/rproject

# # set global options
# COPY --chown=rstudio:rstudio . /home/rstudio/

# # update/install apt package and Japanese font (IPAex)
# RUN apt-get update \
#   && apt-get install -y --no-install-recommends \
#     ssh \
#     libxt-dev \
#     libxml2-dev \
#     libgit2-dev \
#     fonts-ipaexfont \
#     fonts-noto-cjk\
#   && apt-get clean\
#   && sed -i '$d' /etc/locale.gen \
#   && echo "ja_JP.UTF-8 UTF-8" >> /etc/locale.gen \
#     && locale-gen ja_JP.UTF-8 \
#     && /usr/sbin/update-locale LANG=ja_JP.UTF-8 LANGUAGE="ja_JP:ja"\
#   && /bin/bash -c "source /etc/default/locale"\
#   && ln -sf  /usr/share/zoneinfo/Asia/Tokyo /etc/localtime

# # Install R packages
# # RUN Rscript -e "install.packages('vector_of_package_names')"
# # RUN install2.r pkgname1 pkgname2
# # RUN Rscript -e 'BiocManager::install("vector_of_package_names")'
# # example code (below):
# RUN Rscript -e "install.packages(c('agricolae', 'BH', 'BiocManager', 'cowplot', 'devtools', 'doSNOW', 'ggcorrplot',\
#   'ggThemeAssist', 'githubinstall', 'imputeMissings', 'mice', 'minpack.lm', 'missForest', 'mlbench', 'multcomp',\
#   'openxlsx', 'pacman', 'ranger', 'remotes', 'reshape', 'reshape2', 'randomForest', 'rJava', 'rsvg',\
#   'rvg', 'sets', 'sgof', 'tidyverse', 'VIM', 'writexl', 'xlsx', 'XLConnect'))"
# # RUN install2.r --error --skipinstalled\
# #     agricolae BH BiocManager cowplot devtools doSNOW ggcorrplot ggThemeAssist\
# #     githubinstall imputeMissings mice minpack.lm missForest mlbench multcomp \
# #     openxlsx pacman ranger remotes reshape reshape2 randomForest rJava rsvg\
# #     rvg sets sgof tidyverse VIM writexl xlsx XLConnect
# RUN Rscript -e "BiocManager::install(update = TRUE, ask=FALSE, Ncpus=2)"\
#   && Rscript -e 'BiocManager::install("airway")'\
#   && Rscript -e 'BiocManager::install("aLFQ")'\
#   && Rscript -e 'BiocManager::install("BaylorEdPsych")'\
#   && Rscript -e 'BiocManager::install("biomaRt")'\
#   && Rscript -e 'BiocManager::install("BiocStyle")'\
#   && Rscript -e 'BiocManager::install("ComplexHeatmap")'\
#   && Rscript -e 'BiocManager::install("DEP")'\
#   && Rscript -e 'BiocManager::install("DESeq")'\
#   && Rscript -e 'BiocManager::install("DESeq2")'\
#   && Rscript -e 'BiocManager::install("EnhancedVolcano")'\
#   && Rscript -e 'BiocManager::install("enrichR")'\
#   && Rscript -e 'BiocManager::install("genefilter")'\
#   && Rscript -e 'BiocManager::install("gmm")'\
#   && Rscript -e 'BiocManager::install("GO.db")'\
#   && Rscript -e 'BiocManager::install("gplots")'\
#   && Rscript -e 'BiocManager::install("limma")'\
#   && Rscript -e 'BiocManager::install("imsbInfer")'\
#   && Rscript -e 'BiocManager::install("loadTransitonsMSExperiment")'\
#   && Rscript -e 'BiocManager::install("mouse4302.db")'\
#   && Rscript -e 'BiocManager::install("MSstats")'\
#   && Rscript -e 'BiocManager::install("mzR")'\
#   && Rscript -e 'BiocManager::install("org.Hs.eg.db")'\
#   && Rscript -e 'BiocManager::install("org.Mm.eg.db")'\
#   && Rscript -e 'BiocManager::install("PANTHER.db")'\
#   && Rscript -e 'BiocManager::install("PECA")'\
#   && Rscript -e 'BiocManager::install("proDA")'\
#   && Rscript -e 'BiocManager::install("RCyjs")'\
#   && Rscript -e 'BiocManager::install("RforProteomics")'\
#   && Rscript -e 'BiocManager::install("readinteger_binary")'\
#   && Rscript -e 'BiocManager::install("SummarizedExperiment")'\
#   && Rscript -e 'BiocManager::install("SWATH2stats")'