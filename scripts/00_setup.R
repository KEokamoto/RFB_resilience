# Load required packages
library(stringr)
library(tidyverse)
library(readxl)
library(phyloseq)
library(vegan)
library(RColorBrewer)
library(patchwork)
library(DESeq2)
library(pheatmap)
library(scales)
library(grid)

# Global factor levels and styling
gen_levels <- c("G0", "G1", "G2", "G3")
gen_cols <- brewer.pal(4, "Set1")
names(gen_cols) <- gen_levels

trt_levels <- c("ControlA", "ExpA", "ControlB", "ExpB")
trt_shapes <- c("ControlA" = 16, "ExpA" = 17, "ControlB" = 1, "ExpB" = 2)

trt_cols <- c(
  ControlA = "#840032",
  ExpA     = "#f05006",
  ControlB = "#25998f",
  ExpB     = "#f36e98"
)

fill_vals <- c(
  "Pre-disturbance"  = "#d9d9d9",
  "Post-disturbance" = "#bdbdbd"
)

# reproducibility
set.seed(123)
