library(tidyverse)
library(reshape2)
library(grid)
# library(cmapR)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

# senmayo
senmayo = readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Research//Aging/data/human_sen_genesets/senmayo.xlsx")
senmayo = str_to_upper(senmayo$`Gene(human)`)
senmayo.d = rep("up", length(senmayo))
senmayo = c(senmayo, c("CDKN1A", "CDKN2A", "TGFB1"))

# geneAge
genAge = read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/genage_human.csv")
genAge = str_to_upper(genAge$symbol)

# cellAge
cellAge = readr::read_tsv("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/cellage3.tsv")
cellAge.d = cellAge$`Senescence Effect`
cellAge.d[cellAge.d == "Unclear"] = NA
cellAge.d[cellAge.d == "Induces"] = "up"
cellAge.d[cellAge.d == "Inhibits"] = "down"
cellAge = str_to_upper(cellAge$`Gene symbol`)

# CSGene
csgene = read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/csgene_human.txt")
csgene = str_to_upper(csgene$GeneSymb)


# SASP
sasp =readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/SASP.xlsx")
sasp = unlist(as.vector(sasp[3:87, 1]))

# SenSig
sensig.up = readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/SenSig.xlsx", sheet = 1, skip = 12)
sensig.down = readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/SenSig.xlsx", sheet = 2, skip = 12)
sensig = c(sensig.up$`Gene symbol`, sensig.down$`Gene symbol`)
sensig.d = rep(c("up", "down"), c(nrow(sensig.up), nrow(sensig.down)))

# # KEGG cellular senescence
# pathway <- keggGet("hsa04218")[[1]]
# kegg <- pathway$GENE
# kegg = kegg[seq(2,312,2)]
# kegg = sub(";.*", "", kegg)
# # Reactome
# reac = event2Ids(event.id = "R-HSA-2559583")
# reac = reac$geneSymbol
# # sasp
# sasp = event2Ids(event.id = "R-HSA-2559582")
# sasp = sasp$geneSymbol

# Inflammatory network https://doi.org/10.1016/j.molmed.2010.03.003
network = c("CXCL8", "MMP1", "IL1A", "IL1B", "CXCL1", "CXCL3", "CXCL2", "CXCL6", "CXCL5", "AREG",
            "MMP10", "SERPINE1", "PLAU", "IL6", "ICAM1", "IGFBP1", "PLAT", "TNFRSF10C", "IGFBP4",
            "TNFRSF11B", "TNFRSF1B", "TIMP1", "CCL2", "IGFBP3", "SERPINE2", "IGFBP7", "CCL5",
            "IGFBP2", "FAS", "LIF", "FGF2", "VEGFA", "IGFBP6", "EGFR", "AXL", "IL6R", "IGF2",
            "WNT2", "HMGB2", "HMGB1", "HMGB3")

# sennet rec
rec = c("CDKN2A", "CDKN1A", "IL6", "TNF", "IL1A", "IL1B", "SERPINE1", "TGFB1", "TP53", "CCL2", "CCL5",
        "CXCL1", "CXCL8", "MMP12", "CCL3", "HMGB1", "MMP2", "CCL8", "CXCL10", "CXCL2", "GDF15", "IGFBP3",
        "MMP3", "MMP9", "CSF2", "IGFBP7", "LMNB1", "MMP1", "BCL2", "CCL4", "CCL7", "CDKN2B", "CSF1", "CXCL12",
        "FAS", "ICAM1", "IGF1", "IGFBP4", "IL17A", "MAPK14", "MMP13", "NFKB1", "TIMP2", "TNFRSF1B", "VEGFA",
        "CSF3", "CXCL16", "HGF", "ICAM3", "IGFBP2", "IGFBP5", "IL7", "MMP10", "PLAUR", "SPP1", "TNFRSF1A")

# SeneQuest
quest = as.data.frame(read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/data/human_sen_genesets/seneQuest.csv"))
rownames(quest) = quest$`Gene symbol`
quest1 = quest[quest$n_literature >= 4 & quest$percentage >= 0.7, ]
quest2 = quest[quest$n_literature >= 15,]
quest = rbind(quest1, quest2)
quest.d = quest$up_down
quest = rownames(quest)

# tranSig
trans.down = c( "MCUB", "FBL", "HMGB3", "HIST1H1D", "HIST1H1A", "FAM129A", "ANP32B", "PARP1", "LBR", "SSRP1", "TMSB15A", "CBS", "CDC7L", "HIST1H1E", "CBX2", "PTMA", "HIST2H2AB", "ITPRIPL1", "AC074135.1")
trans.up = c("TMEM159", "CHPF2", "SLC9A7", "PLD1", "FAM234B", "DHRS7", "SRPX", "SRPX2", "TNFSF13B", "PDLIM1", "ELMOD1", "CCND3", "TMEM30A", "STAT1", "RND3", "TMEM59", "SARAF", "SLCO2B1", "ARRDC4", "PAM", "WDR78", "CLSTN2", "WDR63", "NCSTN", "SLC16A14", "GPR155", "CLDN1", "JCAD", "BLCAP", "FILIP1L", "TAP1", "TNFRSF10C", "SAMD9L", "SMCO3", "POFUT2", "KIAA1671", "LRP10", "BMS1P9", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "DIO2", "MAP4K3-DT", "AC002480.1", "LINC02154", "TM4SF1-AS1", "PTCHD4", "H2AFJ", "PURPL")
trans = c(trans.up, trans.down)
trans.d = rep(c("up","down"), c(length(trans.up), length(trans.down)))

# construct gene sets
all.genes = unique(c(trans, network, sensig, senmayo, genAge, csgene, cellAge, sasp, quest))
results = data.frame(
  trans = all.genes %in% trans,
  # rec = all.genes %in% rec,
  network = all.genes %in% network,
  sensig = all.genes %in% sensig,
  Senmayo = all.genes %in% senmayo,
  geneAge = all.genes %in% genAge,
  cellAge = all.genes %in% cellAge,
  CSgene = all.genes %in% csgene,
  SASP = all.genes %in% sasp,
  Quest = all.genes %in% quest
)
rownames(results) = all.genes

# get n
results$n = rowSums(results)

# get direction
results$cellAge.d = NA; results$quest.d = NA; results$SenSig.d = NA; results$trans.d = NA
results[sensig, "SenSig.d"] = sensig.d
results[quest, "quest.d"] = quest.d
results[cellAge, "cellAge.d"] = cellAge.d
results[trans, "trans.d"] = trans.d

  
# database agreement
results$percent = rowMeans(results[,c("SenSig.d", "quest.d", "cellAge.d")] == "up", na.rm = T)
results$direction = NA
for (i in 1:nrow(results)) {
  if (is.na(results$percent[i])) {results$direction[i] = NA}
  else if (results$percent[i] > 0.5) {results$direction[i] = "up"}
  else if(results$percent[i] < 0.5) {results$direction[i] = "down"}
  else if(results$percent[i] == 0.5) {results$direction[i] = results$quest.d[i]}
}
results$percent = NULL

gs = results[,c("n", "direction")]
gs$gene_symbol = rownames(gs)
gs = gs %>% arrange(-gs$n)
colnames(gs)[1] = "occurrence"

# gs["JUN", "direction"] = "up"
# gs["EGFR", "direction"] = "up"
# gs["MDM2", "direction"] = "up" # https://www.nature.com/articles/s41598-018-20000-4
# gs["HMGB1", "direction"] = "up"
# gs["MIF", "direction"] = "up"
# gs["VEGFA", "direction"] = "up" # https://www.nature.com/articles/s41598-023-28000-9

write.csv(gs, "./data/coreGS_v2.csv", row.names = T)
write.csv(gs, "./code/DeepScence/DeepScence/data/coreGS_v2.csv", row.names = T)
saveRDS(results, "./data/coreGS_v2.rds")
