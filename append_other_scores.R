library(tidyverse)
library(Seurat)
library(reshape2)
library(GSVA)
library(BiocParallel)
library(AUCell)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

# mtx = readRDS("./data/in_vivo/exp_mtx(converted)/mouse_injury.rds")
mtx = readRDS("./data/in_vitro/exp_mtx/mtx_oskm.rds")

# get gene sets
geneset = readRDS("./data/coreGS_v2.rds")
CoreScence = rownames(geneset)[geneset$n >= 5 & geneset$direction=="up"]
senmayo = rownames(geneset)[geneset$Senmayo]
cellAge = rownames(geneset)[geneset$cellAge & geneset$cellAge.d=="up"]
geneAge = rownames(geneset)[geneset$geneAge]
csgene = rownames(geneset)[geneset$CSgene]
cellsig = rownames(geneset)[geneset$sensig & geneset$SenSig.d == "up"]
quest = rownames(geneset)[geneset$Quest & geneset$quest.d == "up"]
sasp = rownames(geneset)[geneset$SASP]
trans = rownames(geneset)[geneset$trans & geneset$trans.d == "up"]

# gene sets
gs = list(
  CoreScence = CoreScence,
  senmayo=senmayo,
  cellAge = cellAge,
  geneAge = geneAge,
  csgene = csgene,
  cellsig = cellsig,
  quest = quest,
  sasp = sasp,
  trans = trans
)

# run ssGSEA
p = ssgseaParam(mtx, gs)
ssGSEA_scores = gsva(p, BPPARAM = MulticoreParam(workers = 4, progressbar = T))
rownames(ssGSEA_scores) = paste0("ssGSEA_", rownames(ssGSEA_scores))
ssGSEA_scores = as.data.frame(t(as.matrix(ssGSEA_scores)))

# run AUCell
aucell_scores = AUCell_run(mtx, gs)
aucell_scores = getAUC(aucell_scores)
aucell_scores = as.data.frame(t(aucell_scores))
colnames(aucell_scores) = paste0("AUCell_", colnames(aucell_scores))

# single marker
marker_scores = data.frame(
  CDKN1A = mtx["CDKN1A",],
  CDKN2A = mtx["CDKN2A",],
  row.names = colnames(mtx)
)

all_scores = cbind(ssGSEA_scores, aucell_scores, marker_scores)

# binary
all_scores$senmayo_binary = ifelse(all_scores$ssGSEA_senmayo >= quantile(all_scores$ssGSEA_senmayo, 0.9), 1, 0)

write.csv(all_scores, "./data/in_vitro/scores/oskm_other_scores.csv", row.names = T)
