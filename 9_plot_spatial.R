library(tidyverse)
library(reshape2)
library(pROC)
library(Seurat)
library(RColorBrewer)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")

auc_mtx = NULL
cols = NULL
conversion = readRDS("./code/benchmark/method_names_conversion.rds")
datasets = c("mouse_notexin_d2", "mouse_notexin_d5", "human_micro")
for (d in datasets) {
  # read all
  print(d)
  part1 = read.csv(paste0("./data/VADLIATION_DATA/SPATIAL/Current/metas/", d, "_scores_part1.csv"), row.names = 1)
  part2 = read.csv(paste0("./data/VADLIATION_DATA/SPATIAL/Current/metas/", d, "_scores_part2.csv"), row.names = 1)
  part1 = part1[rownames(part2),]
  s = cbind(part1, part2)
  s = na.omit(s)
  
  # fix method names
  colnames(s) = plyr::mapvalues(colnames(s), names(conversion), conversion)
  cont_cols = colnames(s)[!colnames(s) %in% c("SID_binary", "y")]
  
  # calculate auc
  s_cont = s[,cont_cols]
  auc_values <- sapply(s_cont, function(col) {
    roc_obj <- roc(s$y, col, direction = "<") 
    auc(roc_obj) 
  })
  auc_mtx = cbind(auc_mtx, auc_values)
}
colnames(auc_mtx) = datasets 

rownames(auc_mtx)[c(21,22)] = c("SingleMarker:CDKN1A","SingleMarker:CDKN2A")
# auc_mtx = auc_mtx[!rownames(auc_mtx) %in% c("DeepScence:3","DeepScence:4","DeepScence:6"),]
colnames(auc_mtx) = c("Mouse muscle\nnotexin d2", "Mouse muscle\nnotexin d5", "Human brain\nmicroglia")
auc_mtx = auc_mtx[!rownames(auc_mtx) %in% c("DeepScence:Casella et al.", "DeepScence:Freund et al.", "DeepScence:De Cecco et al."),]

# get significance
pvals <- sapply(rownames(auc_mtx), function(method) {
  if(method == "DeepScence") {
    return(NA)  # ds_5. compared to itself gives NA
  } else {
    ttest_result <- t.test(auc_mtx[method, ], auc_mtx["DeepScence", ], paired = TRUE)
    return(ttest_result$p.value)
  }
})
# p.adj <- p.adjust(pvals, method = "bonferroni", n = length(pvals)-1)
p.adj = pvals 
p_sig <- ifelse(is.na(p.adj), NA,
                ifelse(p.adj < 0.001, "***",
                       ifelse(p.adj < 0.01, "**",
                              ifelse(p.adj < 0.05, "*", "NS"))))


df = as.data.frame(auc_mtx)
df$Mean <- apply(df, 1, mean, na.rm = TRUE)
df$method = rownames(df)
df <- df %>% 
  arrange(desc(-Mean))
df_long_auc <- melt(df, id.vars = 'method')
df_long_auc$method = factor(df_long_auc$method, levels = rownames(df))
df_long_auc$p_label <- NA

df_sig <- data.frame(
  method = rownames(df),
  variable = "Significance",  # new x-axis level for p-value annotation
  value = 0.5,         # constant fill value so tile is white
  p_label = p_sig[rownames(df)]    # significance stars (or NA for DeepScence)
)
df_long_auc_new <- rbind(df_long_auc, df_sig)
df_long_auc_new$variable <- factor(df_long_auc_new$variable,
                                   levels = c(setdiff(unique(df_long_auc_new$variable), "Significance"), "Significance"))


# create offset
levels_current <- levels(df_long_auc_new$variable)
n_levels <- length(levels_current)
positions <- 1:n_levels
positions[(n_levels-1):n_levels] <- positions[(n_levels-1):n_levels] + 0.5
position_lookup <- setNames(positions, levels_current)
df_long_auc_new$xpos <- position_lookup[as.character(df_long_auc_new$variable)]

df1 = df_long_auc_new[!df_long_auc_new$method %in% c("DeepScence:3" , "DeepScence:4", "DeepScence:6"),]
df2 = df_long_auc_new[df_long_auc_new$method %in% c("DeepScence:3" , "DeepScence:4", "DeepScence:6", "DeepScence"),]
df = df2

df = df[df$variable!="Significance",]
auc_plot <- ggplot() +
  # Layer for non-Significance rows (using the value for fill)
  geom_tile(data = subset(df, variable != "Significance"),
            aes(x = xpos, y = method, fill = value),
            color = 'black', alpha = 0.7, size = 0.4) +
  # Layer for Significance rows (fill set to white)
  geom_tile(data = subset(df, variable == "Significance"),
            aes(x = xpos, y = method),
            fill = "white",
            color = 'black', alpha = 0.7, size = 0.4) +
  geom_text(data = df,size=5,
            aes(x = xpos, y = method, 
                label = ifelse(variable == "Significance", p_label, sprintf("%.2f", value))),
            color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 11) +
  labs(title = NULL, fill = "AUROC", x = "Datasets", y = "Methods", caption = "Datasets") +
  theme_minimal() +
  scale_x_continuous(position="top", breaks = positions, labels = levels_current, expand = c(0,0)) +
  theme(axis.title.x = element_text(),
        axis.title.x.top = element_blank(),  
        axis.text.y = element_text(color = "black", size = 10, 
                                   face = ifelse(levels(factor(df$method)) == "DeepScence", "bold", "plain")),
        axis.text.x.top = element_text(color = "black", 
                                       face = c(rep('plain',3), "bold"), 
                                       size = 10, vjust = 0.3),
        axis.title.y = element_text(hjust = 0.8, size = 13),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.justification = c(1, 1.07),
        plot.caption = element_text(hjust = 0.5, size = 13)
  ) + 
  guides(fill="none")

auc_plot

# normal aging case
d = "mouse_aging"
adata = readRDS(paste0("./data/VADLIATION_DATA/SPATIAL/Current/rds/", d, ".rds"))
part1 = read.csv(paste0("./data/VADLIATION_DATA/SPATIAL/Current/metas/", d, "_scores_part1.csv"), row.names = 1)
part2 = read.csv(paste0("./data/VADLIATION_DATA/SPATIAL/Current/metas/", d, "_scores_part2.csv"), row.names = 1)
part1 = part1[rownames(part2),]
s = cbind(part1, part2)
cont_cols = colnames(s)[!colnames(s) %in% c("SID_binary", "y")]
s_cont = s[,cont_cols]
age_group = ifelse(adata$y==1, "Old", "Young")

get_or <- function(score_df, col, age_group, q) {
  score_df$age_group = age_group
  logORs <- c()
  f1s = c()
  rrs = c()
  for (quant in q) {
    threshold <- quantile(score_df[,col], probs = 1 - quant / 100, na.rm = TRUE)
    
    # Create contingency table
    a <- sum(score_df[,col] > threshold & score_df$age_group == "Old")
    b <- sum(score_df[,col] > threshold & score_df$age_group == "Young")
    c <- sum(score_df[,col] <= threshold & score_df$age_group == "Old")
    d <- sum(score_df[,col] <= threshold & score_df$age_group == "Young")
    
    # print table
    contingency_table <- data.frame(
      "Old" = c(a, c),
      "Young" = c(b, d),
      row.names = c("Top Quantile", "Bottom Quantile")
    )
    # get rr
    rr = (a/(a+b)) / (c/(c+d))
    rrs = c(rrs, rr)
  }
  return(rrs)
}
rr_table = data.frame()
q <- c(1, 1)
for (col in colnames(s_cont)) {
  rrs <- get_or(s_cont, col, age_group, q)
  rr_table = rbind(rr_table, rrs)
}
colnames(rr_table) = q
rownames(rr_table) = colnames(s_cont)
out = sort(rowMeans(rr_table), decreasing = T)
res = data.frame(
  method = names(out),
  avg_rrs = out
)
write.csv(res, "./data/NORMAL_AGING/ST_res.csv")


# stereoseq example
folder = "./data/additional_sc/stereo-seq/scored_meta/"
files <- list.files(folder, full.names = TRUE)
global_min <- Inf
global_max <- -Inf
for (f in files) {
  m <- read.csv(f, row.names = 1)
  # Compute min and max for ds_5. in this file (ignoring NA)
  local_min <- min(m$ds_5., na.rm = TRUE)
  local_max <- max(m$ds_5., na.rm = TRUE)
  if(local_min < global_min) global_min <- local_min
  if(local_max > global_max) global_max <- local_max
}
print(c(global_min, global_max))

l <- list()
samples <- c("Severe AD", "Normal Aging", "Normal Aging", "Severe AD")
names(samples) <- c("A02092E1_m", "B01806B6_m", "B01809A4_m", "B02008C6_m")
SpatialColors <- colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))
for (f in files) {
  sample <- gsub(".csv", "", basename(f))
  print(sample)
  m <- read.csv(f, row.names = 1)
  p <- ggplot() +
    geom_point(data = subset(m, is.na(ds_5.)),
               aes(x = x, y = y),
               color = "lightgrey", size = 0.1) +
    geom_point(data = subset(m, !is.na(ds_5.)),
               aes(x = x, y = y, color = ds_5.),
               size = 0.3) +
    scale_color_gradientn(colors = SpatialColors(100),
                          limits = c(global_min, global_max)) +
    theme_bw() +
    labs(x = NULL, y = NULL, color = "DeepScence",
         title = paste0(gsub("_m", "", sample), " (", samples[sample], ")"))+
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  l[[sample]] <- p
}
l = l[c("B01806B6_m","B02008C6_m", "B01809A4_m", "A02092E1_m")]


