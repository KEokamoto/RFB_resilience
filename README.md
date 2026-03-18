# RFB_resilience

Cross-generational gut microbiome dynamics in Tribolium castaneum under environmental and dietary perturbations.

This repository contains all analysis code used to generate figures and results for the manuscript:

"Resilience and reliability of gut microbiome recovery across generations in Tribolium castaneum"

---------------------------------------------------------------------

Project Overview

Microbial communities often experience environmental disturbances, yet recovery does not necessarily return to a single predictable state. This project investigates:

- Microbiome restructuring following transition from commercial to laboratory conditions
- Cross-generational stability (G0–G3)
- Effects of dietary perturbations (whole wheat vs oat)
- Distinction between resilience and reliability in microbiome recovery

---------------------------------------------------------------------

Repository Structure

data/
├── raw/                # Raw input data (not tracked if large)
├── processed/          # Processed phyloseq objects and tables

scripts/
├── 00_setup.R          # Load packages and define global variables
├── 01_processing.R     # Data import and preprocessing
├── 02_analysis.R       # Statistical analyses (alpha, beta, DESeq2)
├── 03_figures.R        # Figure generation (Fig 2–7)

figures/                # Final publication-ready figures
results/                # Statistical outputs and intermediate tables
---------------------------------------------------------------------

Workflow

1. Data processing
   - Filtering, taxonomic aggregation, phyloseq object construction  
   → scripts/01_processing.R

2. Statistical analysis
   - Alpha diversity (Observed, Shannon)
   - Beta diversity (Bray–Curtis, PERMANOVA)
   - Differential abundance (DESeq2)
   → scripts/02_analysis.R

3. Figure generation
   - Fig 2: Baseline alpha diversity (G0)
   - Fig 3: Order-level composition (flour vs gut)
   - Fig 4: Alpha diversity across generations
   - Fig 5: PCoA ordination
   - Fig 6: Differential abundance heatmaps
   - Fig 7: Resilience dynamics
   → scripts/03_figures.R

---------------------------------------------------------------------

Required R Packages

- phyloseq
- DESeq2
- vegan
- ggplot2
- dplyr
- tidyr
- pheatmap
- patchwork
- RColorBrewer
- scales

---------------------------------------------------------------------

Reproducibility

To reproduce all figures:

source("scripts/00_setup.R")
source("scripts/01_processing.R")
source("scripts/02_analysis.R")
source("scripts/03_figures.R")

---------------------------------------------------------------------

Author

Esther Okamoto  
Undergraduate Student | Agnes Scott College
PhD Student, Biology | California Institute of Technology  

---------------------------------------------------------------------

Notes

- Raw sequencing data are not included due to file size constraints
- Processed data objects may be provided upon request
