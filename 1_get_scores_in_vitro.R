library(Seurat)
library(tidyverse)
library(escape)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

pp = "Standard"
datasets = c("hayflick","eto1", "hca", "huvec", "notch", "oskm")
# gs = read.csv("./data/coreGS_v2.csv", row.names = 1)
gs = readRDS("./data/coreGS_v2.rds")
gs = list(
  "trans" = rownames(gs)[gs$trans.d=="up"],
  "network" = rownames(gs)[gs$network],
  "sensig" = rownames(gs)[gs$SenSig.d=="up"],
  "Senmayo" = rownames(gs)[gs$Senmayo],
  "geneAge" = rownames(gs)[gs$geneAge],
  "cellAge" = rownames(gs)[gs$cellAge.d=="up"],
  "CSgene" = rownames(gs)[gs$CSgene],
  "SASP" = rownames(gs)[gs$SASP],
  "Quest" = rownames(gs)[gs$quest.d=="up"],
  "CoreScence" = rownames(gs)[gs$n >= 5]
)

for (d in datasets) {
  adata = readRDS(paste0("./data/VADLIATION_DATA/IN_VITRO/", pp, "/rds/", d, ".rds"))
  adata@active.assay = "RNA"
  adata = NormalizeData(adata)
  
  # ssgsea
  adata = runEscape(
    adata,
    gene.sets = gs,
    method = "ssGSEA",
    normalize = T,
    new.assay.name = "escape_ssGESA"
  )
  
  # AUCell
  adata = runEscape(
    adata,
    gene.sets = gs,
    method = "AUCell",
    normalize = T,
    new.assay.name = "escape_AUCell"
  )
  
  # single-gene
  single_exp = t(as.matrix(GetAssayData(adata, layer = "data")[c("CDKN1A", "CDKN2A"),]))
  
  # put together and save
  s1 = t(adata[["escape_ssGESA"]]$data)
  colnames(s1) = paste0(colnames(s1), "_ssGSEA")
  s2 = t(adata[["escape_AUCell"]]$data)
  colnames(s2) = paste0(colnames(s2), "_AUCell")
  s = cbind(s1, s2, single_exp)
  
  write.csv(s, paste0("./data/VADLIATION_DATA/IN_VITRO/", pp, "/metas/", d, "_scores_part1.csv"), row.names = T)
}

