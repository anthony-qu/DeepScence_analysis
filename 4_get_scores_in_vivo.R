library(Seurat)
library(tidyverse)
library(escape)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

# # get "y" for each dataset, overwrite data
# mouse_heart_ECs = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_heart_ECs.rds")
# mouse_heart_CMs = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_heart_CMs.rds")
# mouse_tauopathy = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_tauopathy.rds")
# mouse_lung_cancer = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_lung_cancer.rds")
# human_IPF = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/human_IPF.rds")
# human_oral = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/human_oral.rds")
# mouse_muscle = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_muscle.rds")
# mouse_testes = readRDS("./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_testes.rds")
# 
# mouse_heart_ECs$y = ifelse(mouse_heart_ECs$condition=="TAC2w", 1, 0)
# mouse_heart_CMs = mouse_heart_CMs[,mouse_heart_CMs$condition!="TAC2w+CDKi"]
# mouse_heart_CMs$y = ifelse(mouse_heart_CMs$condition=="TAC2w", 1, 0)
# mouse_tauopathy$y = ifelse(mouse_tauopathy$condition=="tau", 1,0)
# mouse_lung_cancer$y = ifelse(mouse_lung_cancer$mCherry.Count>0, 1,0)
# human_IPF$y = ifelse(human_IPF$health_status=="IPF", 1, 0)
# human_oral$y = ifelse(human_oral$Condition=="periodontitis", 1, 0)
# mouse_muscle$y = ifelse(mouse_muscle$stim == "CTX", 1, 0)
# mouse_testes$y = ifelse(mouse_testes$stim == "O50t", 1,0)
# 
# saveRDS(mouse_heart_ECs, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_heart_ECs.rds")
# saveRDS(mouse_heart_CMs, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_heart_CMs.rds")
# saveRDS(mouse_tauopathy, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_tauopathy.rds")
# saveRDS(mouse_lung_cancer, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_lung_cancer.rds")
# saveRDS(human_IPF, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/human_IPF.rds")
# saveRDS(human_oral, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/human_oral.rds")
# saveRDS(mouse_muscle, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_muscle.rds")
# saveRDS(mouse_testes, "./data/VADLIATION_DATA/IN_VIVO/Current/rds/mouse_testes.rds")

# # add mouse genes to gs
# gs = read.csv("./data/coreGS_v2.csv", row.names = 1)
# c = read.csv("./data/in_vivo/gene_convert.csv")
# gs$mouse_gene = plyr::mapvalues(rownames(gs), c$converted, c$original)
# write.csv(gs, "./data/coreGS_v2.csv", row.names = T)

pp = "Standard"
datasets = c("human_IPF", "human_oral", "mouse_muscle", "mouse_testes", "mouse_lung_cancer","mouse_tauopathy", "mouse_heart_CMs", "mouse_heart_ECs")
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
  adata = readRDS(paste0("./data/VADLIATION_DATA/IN_VIVO/", pp, "/rds/", d, ".rds"))
  adata@active.assay = "RNA"
  adata = NormalizeData(adata)
  adata$y = adata$SnC
  
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
  
  write.csv(s, paste0("./data/VADLIATION_DATA/IN_VIVO/", pp, "/metas/", d, "_scores_part1.csv"), row.names = T)
  
}
