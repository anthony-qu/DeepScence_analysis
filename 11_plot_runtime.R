library(tidyverse)
library(Seurat)
library(scDesign3)
library(SingleCellExperiment)
rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Aging/")



res = read.csv("./data/runtime_res.csv", row.names = 1)
res.gpu = read.csv("./data/runtime_res_gpu.csv", row.names = 1)

res_summary_cpu <- res %>%
  group_by(ncell) %>%
  summarize(
    mean_runtime_sec = mean(runtime_sec, na.rm = TRUE),
    mean_memory_mb   = mean(net_memory_mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    runtime_hr = mean_runtime_sec / 3600,
    memory_gb  = mean_memory_mb   / 1024,
    type = "CPU"
  )

# GPU summary with 'type'
res_summary_gpu <- res.gpu %>%
  group_by(ncell) %>%
  summarize(
    mean_runtime_sec = mean(runtime_sec, na.rm = TRUE),
    mean_memory_mb   = mean(net_memory_mb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    runtime_hr = mean_runtime_sec / 3600,
    memory_gb  = mean_memory_mb   / 1024,
    type = "GPU"
  )

# Combine both summaries
res_summary_all <- bind_rows(res_summary_cpu, res_summary_gpu)

# Plot runtime with color mapped to 'type'
p_runtime <- ggplot(res_summary_all, aes(x = ncell, y = runtime_hr, color = type)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  labs(
    title = "Runtime",
    x = "Number of cells",
    y = "Runtime (hours)",
    color ="Device"
  ) + 
  theme_bw()

# Plot memory with color mapped to 'type'
p_memory <- ggplot(res_summary_all, aes(x = ncell, y = memory_gb, color = type)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_line(aes(group = type)) +
  scale_x_log10() +
  labs(
    title = "Memory Usage",
    x = "Number of cells",
    y = "Memory (GB)",
    color ="Device"
  ) + 
  theme_bw()

p_runtime|p_memory


