# RFB_resilience

Cross-generational restructuring of the *Tribolium castaneum* gut microbiome under environmental and dietary perturbations.

This repository contains all analysis code used to generate figures and statistical analyses for the manuscript:

**“Resilient but Unreliable: Cross-Generational Restructuring of the *Tribolium castaneum* Gut Microbiome under Contrasting Nutritional Regimes”**

---

## Project Overview

Microbial communities frequently experience environmental disturbance, yet recovery does not necessarily return communities to a single predictable state.

Using a four-generation experimental design (G0–G3), this project investigates:

- Rapid microbiome restructuring following transition from commercial to laboratory conditions
- Recovery dynamics following perturbation
- Effects of dietary shifts (whole wheat vs oat flour)
- Distinction between:
  - **Resilience** = stabilization/recovery following disturbance
  - **Reliability** = convergence toward a shared community endpoint
- Environmental filtering and history-dependent microbiome assembly

---

## Repository Structure

```text
data/
├── raw/                     # Raw input data (not tracked if large)
├── processed/               # Processed phyloseq objects and tables

scripts/
├── 00_setup.R
├── 01_processing.R
├── 02_analysis.R
├── 02_analysis_integrated.R
├── 03_figures.R
├── 03_figure_w_supp.R
├── 03_figures_integrated.R
├── run_all.R

figures/                     # Publication-ready figures
results/                     # Statistical outputs and intermediate tables
