# Load required packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(readxl)
library(RColorBrewer)
library(patchwork)
library(DESeq2)
library(pheatmap)
library(scales)

# Define global variables
gen_levels <- c("G0","G1","G2","G3")
gen_cols <- brewer.pal(4, "Set1")
names(gen_cols) <- gen_levels

trt_levels <- c("ControlA","ExpA","ControlB","ExpB")
trt_shapes <- c("ControlA"=16, "ExpA"=17, "ControlB"=1, "ExpB"=2)
