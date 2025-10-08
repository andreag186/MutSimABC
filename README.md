# MutSimABC: A Simulation-Based Approximate Bayesian Computational Framework for Mutation Rate Inference in Long-Lived Trees

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
- [License](#license)
- [Contact](#contact)

---

## Overview

The rate and distribution of mutations accumulated in long-lived trees is governed by the structure and developmental dynamics of meristematic cells. Understanding and quantifying somatic mutation rates is essential for understanding plant evolution, optimizing breeding programs in forestry and agriculture, and managing genetic resources in conservation.
However, existing WGS methods struggled distinguishing true somatic mutations from sequencing errors and alignment artifacts. An alternative *Phylogenomic method* introduced by [Orr et al. 2020](https://doi.org/10.1098/rspb.2019.2364)) aimed to estimate mutation rates by assuming mutations follow tree topology in order to reduce errors and presevere only variants with 'true biological signal'. However, this assumption is violated when meristem dynamics create patterns independent of branching structure. We introduce MutSimABC as an alternative prototype framework which addreses limitations of existing methods. 

**MutSimABC** addresses these limitations by:
1. **Simulating** mutation accumulation under mechanistic models of meristem elongation and branching (Tomimoto & Satake 2023)
2. **Inferring** mutation rates and developmental meristem parameters using likelihood-free ABC-Reject inference
3. **Validating** against synthetic data with known ground truth across diverse tree architectures

The framework extends the simulation model of [Tomimoto & Satake (2023)](https://doi.org/10.1016/j.jtbi.2023.111465) to handle arbitrary tree topologies and integrates it into a Bayesian inference pipeline.

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


