# MutSimABC: A Simulation-Based Approximate Bayesian Computational Framework for Mutation Rate Inference in Long-Lived Trees

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)


**MutSimABC** is a simulation-based inference framework for estimating somatic mutation rates in long-lived trees using Approximate Bayesian Computation (ABC). The framework explicitly accounts for mechanistic meristem dynamics (elongation and branching processes) that govern mutation accumulation and distribution across tree architectures.

This repository accompanies the paper:
> Grecu, A.M. et al. (2025). "MutSimABC: A Simulation-Based Approximate Bayesian Computational Framework for Mutation Rate Inference in Long-Lived Trees" *IEEE/ACM Transactions on Computational Biology and Bioinformatics* (in review).

---

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Repository Structure](#repository-structure)
- [Usage Guide](#usage-guide)
  - [ABC Inference on Empirical Data](#abc-inference-on-empirical-data)
  - [Validation Framework](#validation-framework)
- [Understanding the Code](#understanding-the-code)
- [Defining Custom Tree Topologies](#defining-custom-tree-topologies)
- [Understanding Model Parameters](#understanding-model-parameters)
- [Computational Requirements](#computational-requirements)
- [Reproducing Paper Results](#reproducing-paper-results)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

The rate and distribution of mutations accumulated in long-lived trees is governed by the structure and developmental dynamics of meristematic cells. Understanding and quantifying somatic mutation rates is essential for understanding plant evolution, optimizing breeding programs in forestry and agriculture, and managing genetic resources in conservation.
However, existing WGS methods struggled distinguishing true somatic mutations from sequencing errors and alignment artifacts. An alternative *Phylogenomic method* introduced by [Orr et al. 2020](https://doi.org/10.1098/rspb.2019.2364) aimed to estimate mutation rates by assuming mutations follow tree topology in order to reduce errors and presevere only variants with 'true biological signal'. However, this assumption is violated when meristem dynamics create patterns independent of branching structure. We introduce MutSimABC as an alternative prototype framework which addreses limitations of existing methods. 

**MutSimABC** addresses these limitations by:
1. **Simulating** mutation accumulation under mechanistic models of meristem elongation and branching (Tomimoto & Satake 2023)
2. **Inferring** mutation rates and developmental meristem parameters using likelihood-free ABC-Reject inference
3. **Validating** against synthetic data with known ground truth across diverse tree architectures

The framework extends the simulation model of [Tomimoto & Satake (2023)](https://doi.org/10.1016/j.jtbi.2023.111465) to handle arbitrary tree topologies and integrates it into a Approximate Bayesian Computation (ABC)- Reject/Accept inference pipeline.

## Installation

### Prerequisites
- Python 3.11.4 or higher
- pip package manager
- (Optional) Access to HPC cluster with SLURM for large-scale runs

### Clone Repository
```bash
git clone https://github.com/andreag186/MutSimABC.git
cd MutSimABC
```

### Install Dependencies 
```bash
pip install -r requirements.txt
```

### Verify Installation
```bash
python -c "import numpy, scipy, pandas, arviz; print('MutSimABC Dependencies Installed Succesfully')"
```

## Quick Start

### Example 1: Run ABC-Reject inference on E. melliodora pre-DNG data
```bash
python abc_non_dng.py
```
This will:

Sample 5,000 parameter sets from prior distributions
-Run simulations for each parameter set using the E. melliodora tree topology
-Accept parameter sets where simulated mutation distribution is within ε=20 of observed distribution of mutations prior to topological filtering 
-Save accepted samples to results/accepted_samples_*.csv

### Example 2: Run ABC inference on E. melliodora post-DNG (topologically filtered) data
```bash
python abc_phylo.py
```
Similar to above but uses the 90 mutations that remained after topological filtering.

### Example 3: Validate ABC framework with known ground truth
```bash
# To validate synthetic sample 1
python abc_validation.py 1
```
This will:

-Run ABC-Reject validation for sample 1 from abc_validation_samples.csv
-Sample parameters (StD, biasVar, input_mut) from prior distributions
-Simulate mutation distributions using the tree topology defined for that sample
-Compare simulated and observed mutation distributions, accepting parameter sets where the Euclidean distance < ε = 20
-Continue sampling (up to ~5,000 trials) until 100 accepted samples are obtained
-Save results to:
result_valid/accepted_samples_sample1_epsilon20_<timestamp>.csv

To validate a different sample, change the task ID argument (e.g., python abc_validation.py 42 for sample 42- there are 200 samples total(index starting at 0)).

