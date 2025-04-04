library(Seurat)
library(tidyverse)
library(escape)
rm(list = ls())
options(Seurat.object.assay.version = "v3")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

# get stereo-seq rds first, then get "y" for all data
# counts = read.csv("./data/additional_sc/stereo-seq/counts.csv", row.names = 1)
# m = read.csv("./data/additional_sc/stereo-seq/meta.csv", row.names = 1)
# counts <- t(as.matrix(counts))
# adata <- CreateSeuratObject(counts = counts, meta.data = m)
# adata$y = ifelse(adata$condition == "Severe", 1,0)
# saveRDS(adata, "./data/additional_sc/stereo-seq/micro.rds")

# # ctx data
# adata = readRDS("./data/in_vivo/sp_mouse_injury/ctx_scored_final.rds")
# adata@meta.data = adata@meta.data[,1:8]
# adata = adata[,adata$condition != "Newly repaired"]
# adata$y = ifelse(adata$condition == "Injured", 1,0)
# adata[["RNA"]] = adata[["Spatial"]]; adata@active.assay = "RNA"
# adata[["Spatial"]] = NULL
# saveRDS(adata, "./data/in_vivo/sp_mouse_injury/ctx.rds")
# 
# # notexin data
# adata = readRDS("./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_notexin_d2.rds")
# adata$y = ifelse(adata$is.injury, 1, 0)
# adata[["RNA"]] = adata[["Spatial"]]; adata@active.assay = "RNA"
# adata[["Spatial"]] = NULL
# saveRDS(adata, "./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_notexin_d2.rds")
# adata = readRDS("./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_notexin_d5.rds")
# adata$y = ifelse(adata$is.injury, 1, 0)
# adata[["RNA"]] = adata[["Spatial"]]; adata@active.assay = "RNA"
# adata[["Spatial"]] = NULL
# saveRDS(adata , "./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_notexin_d5.rds")

# # mouse aging data
# adata = readRDS("./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_aging.rds")
# adata$y = ifelse(startsWith(adata$sample, "18M"), 1,0)
# saveRDS(adata, "./data/VADLIATION_DATA/SPATIAL/Current/rds/mouse_aging.rds")

pp = "Current"
datasets = c("mouse_notexin_d2", "mouse_notexin_d5", "human_micro", "mouse_ctx", "mouse_aging")
gs = read.csv("./data/coreGS_v2.csv", row.names = 1)

gs_human = list(
  "trans" = rownames(gs)[gs$trans],
  "network" = rownames(gs)[gs$network],
  "sensig" = rownames(gs)[gs$sensig],
  "Senmayo" = rownames(gs)[gs$Senmayo],
  "geneAge" = rownames(gs)[gs$geneAge],
  "cellAge" = rownames(gs)[gs$cellAge],
  "CSgene" = rownames(gs)[gs$CSgene],
  "SASP" = rownames(gs)[gs$SASP],
  "Quest" = rownames(gs)[gs$Quest],
  "CoreScence" = rownames(gs)[gs$n >= 5]
)
gs_mouse = list(
  "trans" = gs$mouse_gene[gs$trans],
  "network" = gs$mouse_gene[gs$network],
  "sensig" = gs$mouse_gene[gs$sensig],
  "Senmayo" = gs$mouse_gene[gs$Senmayo],
  "geneAge" = gs$mouse_gene[gs$geneAge],
  "cellAge" = gs$mouse_gene[gs$cellAge],
  "CSgene" = gs$mouse_gene[gs$CSgene],
  "SASP" = gs$mouse_gene[gs$SASP],
  "Quest" = gs$mouse_gene[gs$Quest],
  "CoreScence" = gs$mouse_gene[gs$n >= 5]
)

for (d in datasets) {
  print(paste0("processing ", d, "..."))
  adata = readRDS(paste0("./data/VADLIATION_DATA/SPATIAL/", pp, "/rds/", d, ".rds"))
  adata = NormalizeData(adata)
  
  if (startsWith(d, "mouse_")) {
    gs.use = gs_mouse
    single_exp = t(as.matrix(GetAssayData(adata, layer = "data")[c("Cdkn1a", "Cdkn2a"),]))
  } else {
    gs.use = gs_human
    single_exp = t(as.matrix(GetAssayData(adata, layer = "data")[c("CDKN1A", "CDKN2A"),]))}
  
  # ssgsea
  adata = runEscape(
    adata,
    gene.sets = gs.use,
    method = "ssGSEA",
    normalize = T,
    new.assay.name = "escape_ssGESA"
  )
  
  # AUCell
  adata = runEscape(
    adata,
    gene.sets = gs.use,
    method = "AUCell",
    normalize = T,
    new.assay.name = "escape_AUCell"
  )
  
  
  # put together and save
  s1 = t(adata[["escape_ssGESA"]]$data)
  colnames(s1) = paste0(colnames(s1), "_ssGSEA")
  s2 = t(adata[["escape_AUCell"]]$data)
  colnames(s2) = paste0(colnames(s2), "_AUCell")
  s = cbind(s1, s2, single_exp, adata@meta.data["y"])
  
  write.csv(s, paste0("./data/VADLIATION_DATA/SPATIAL/Current/metas/", d, "_scores_part1.csv"), row.names = T)
  
}

