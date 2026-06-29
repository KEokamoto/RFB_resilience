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
├── 00_setup.R                 # Package loading, global colors/factor levels
├── 01_processing.R            # Data import, metadata cleaning, phyloseq objects
├── 02_analysis_integrated.R    # Alpha/beta diversity, DA, SIMPER, overlap stats
├── 03_figures_integrated.R     # Publication figures
├── run_all.R                  # Runs full workflow

figures/                     # Publication-ready figures
results/                     # Statistical outputs and intermediate tables
