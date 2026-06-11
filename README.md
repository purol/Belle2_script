# Belle II $B \to X_s \nu \bar{\nu}$ Search: Advanced Statistical Modeling & ML Codebase

[![arXiv](https://img.shields.io/badge/arXiv-2511.10980-b31b1b.svg)](https://arxiv.org/abs/2511.10980)
[![Published in PRL](https://img.shields.io/badge/Published%20in-Physical%20Review%20Letters-blue)](https://journals.aps.org/prl/abstract/10.1103/kf73-hw61)
[![DOI](https://img.shields.io/badge/DOI-10.1103%2Fkf73--hw61-blue)](https://doi.org/10.1103/kf73-hw61)

## Overview
This repository contains the comprehensive research codebase for the rare decay search $B \to X_s \nu \bar{\nu}$ at the Belle II experiment. 

**This is the codebase that produced the results for the publication: [arXiv:2511.10980](https://arxiv.org/abs/2511.10980).**

Over the course of multiple years of continuous development, this project evolved to handle complex, high-dimensional physics data. It utilizes advanced Machine Learning techniques (FastBDT) and rigorous statistical modeling to extract incredibly faint signal signatures from massive background noise.

## Tech Stack & Core Competencies
- **Machine Learning:** `FastBDT` (Boosted Decision Trees) for binary classification of highly imbalanced data.
- **Statistical Inference:** HistFactory and ROOT-based Toy Monte Carlo (MC) simulations for profile likelihood fits and confidence interval calculations.
- **Languages:** C/C++ (core logic), Python (Data Handling & data visualization), Bash/Shell (Automation).

## Repository Structure
The project is structured into modular components for scalability and maintainability:

- **`src/`** : The core analysis directory, containing:
  - `Xsnunubar_analysis/`, `XsKsKs_analysis/`, `Xspp_analysis/`: Modularized C++ codebase for specific decay channel analyses.
  - `python/`, `matlab/`: Multi-language scripts for data analysis and processing.
- **`archive/`** : Deprecated scripts and legacy models securely version-controlled.

## Key Data Science & Engineering Features

### 1. Machine Learning Classification
- Designed and optimized Machine Learning pipelines to classify heavily imbalanced datasets (Signal vs. Massive Background).
- Conducted hyperparameter grid searches and over-training validation to maximize the AUC and Punzi Figure of Merit (FOM).

### 2. Rigorous Statistical Analysis
- Modeled complex Probability Density Functions (PDFs) and performed multi-dimensional maximum likelihood fits.
- Evaluated and integrated comprehensive systematic uncertainties into the final statistical model to ensure robust, publication-ready results.

### 3. Toy Monte Carlo Validations
- Designed and executed thousands of "Toy MC" pseudo-experiments to validate the mathematical stability of the fitting algorithms and calculate precise Upper Limits (CLs method).

## Impact & Achievements
- **Scientific Publication:** The codebase successfully processed Petabyte-scale experiment data and produced the robust statistical results published in **[arXiv:2511.10980](https://arxiv.org/abs/2511.10980)**.
- **Long-term Engineering & Persistence:** Built and iteratively refined a complex computational framework from scratch. Successfully drove a multi-year, large-scale data science project from initial data exploration to a final result.
