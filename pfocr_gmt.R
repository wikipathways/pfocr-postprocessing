library(magrittr)
library(dplyr)
library(rWikiPathways)

genes <- readRDS("../pfocr_genes_draft.rds")
figs <- readRDS("../pfocr_figures_draft.rds")

##################################################################
## Default ##
##############
# Starting with all PFOCR, then applying 3+ filter

genes.cnt<- genes %>%
  distinct(figid,entrez) %>%
  group_by(figid) %>%
  summarise(cnt=n())

figid_list3<- as.list(genes.cnt %>%
                        filter(cnt >=3) %>%
                        select(figid))[[1]]

gmt.df<-genes %>%
  distinct(figid,entrez) %>%
  dplyr::filter(figid %in% figid_list3) %>%
  dplyr::left_join(figs, by="figid") %>%
  dplyr::select(figid, figtitle, entrez)

# This take ~7 minutes to complete
rWikiPathways::writeGMT(gmt.df, "../pfocr_genes.gmt")

##################################################################
## hgnc3 ##
##############
# Removing bioentities, previous and alias matches, then
# filtering for 3+ genes

genes.nobe0 <- genes %>%
  dplyr::filter(source == "hgnc_symbol")

genes.cnt<- genes.nobe0 %>%
  distinct(figid,entrez) %>%
  group_by(figid) %>%
  summarise(cnt=n())

figid_list3<- as.list(genes.cnt %>%
                        filter(cnt >=3) %>%
                        select(figid))[[1]]

gmt.df<-genes.nobe0 %>%
  distinct(figid,entrez) %>%
  dplyr::filter(figid %in% figid_list3) %>%
  dplyr::left_join(figs, by="figid") %>%
  dplyr::select(figid, figtitle, entrez)

# This take ~1 minutes to complete
rWikiPathways::writeGMT(gmt.df, "../pfocr_genes_hgnc3.gmt")


