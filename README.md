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

### Example 1: Run ABC-Reject inference on *E. melliodora* pre-DNG data
```bash
python abc_non_dng.py
```
This will:

- Sample 5,000 parameter sets from prior distributions
- Run simulations for each parameter set using the E. melliodora tree topology
- Accept parameter sets where simulated mutation distribution is within ε=20 of observed distribution of mutations prior to topological filtering 
- Save accepted samples to results/accepted_samples_*.csv

### Example 2: Run ABC inference on *E. melliodora* post-DNG (topologically filtered) data
```bash
python abc_phylo.py
```
Similar to above but uses the 90 mutations that remained after topological filtering.

### Example 3: Validate ABC framework with known simulated 'observed data' (*synthetic sample*)
```bash
# To validate synthetic sample 1
python abc_validation.py 1
```
This will:

- Run ABC-Reject validation for sample 1 from abc_validation_samples.csv
- Sample parameters (StD, biasVar, input_mut) from prior distributions
- Simulate mutation distributions using the tree topology defined for that sample
- Compare simulated and observed mutation distributions, accepting parameter sets where the Euclidean distance < ε = 20
- Continue sampling (up to ~5,000 trials) until 100 accepted samples are obtained
- Save results to:
result_valid/accepted_samples_sample1_epsilon20_<timestamp>.csv

To validate a different sample, change the task ID argument (e.g., `python abc_validation.py 42` for sample 42- there are 200 samples total (index starting at 0).

## Repository Structure
```
MutSimABC/
├── README.md                      # This file
├── requirements.txt               # Python dependencies
│
├── abc_non_dng.py                # ABC inference for pre-DNG E. melliodora data
├── abc_phylo.py                  # ABC inference for post-DNG E. melliodora data
├── abc_validation.py             # Validation framework with known ground truth
│
├── abc_validation_samples_combined.csv    # Pre-generated validation samples (169 samples)
│
├── results/                      # Output directory for pre-DNG accepted samples
├── result_phylo/                 # Output directory for post-DNG accepted samples
└── result_valid/                 # Output directory for validation accepted samples
```
## Usage Guide 

### ABC Inference on Empirical *E. melliodora* Data
Each script ( `abc_non_dng.py` and `abc_phylo.py`) contains a complete ABC-Reject implementation that:
1. Samples parameters (μ = mutation rate, StD = meristem elongation stochasticity, σ = bias of branching) from prior distributions
2. Simulates mutation distributions using Tomimoto & Satake's mechanistic model
3. Compares simulated vs. observed distributions using Euclidean distance
4. Accepts parameter sets when distance ≤ ε (epsilon- set Euclidean distance at which samples are accepted)

**Key Parameters to Modify**
In `abc_non_dng.py` or `abc_phylo.py`, find the `abc_rejection()` function at the bottom:
```python
if __name__ == "__main__":
    num_samples = 5000      # Number of trials
    epsilon = 20            # Acceptance threshold
    job_id = os.getenv('SLURM_ARRAY_TASK_ID', 'default_job')
    batch_id = sys.argv[1] if len(sys.argv) > 1 else 'default_batch'
    
    abc_rejection(num_samples, epsilon, job_id, batch_id)
```
To customise:
- `num_samples`: Increase for more thorough posterior sampling (10,000 recommended for final analysis)
- `epsilon`: Adjust acceptance threshold (lower = stricter, higher = more permissive- can conduct own sensitivity analysis)

**Prior Distributions**
Modify the `sample_prior()` function to change parameter ranges:
```python
def sample_prior():
    StD = np.random.uniform(1, 5)              # Elongation parameter
    biasVar = np.random.uniform(0.5, 10)       # Branching bias
    input_mut = np.random.uniform(1e-11, 9e-9) # Mutation rate per site per year
    return StD, biasVar, input_mut
```
**Running on HPC with SLURM**
For parrallel execution on computing clusters 
```bash
#!/bin/bash
#SBATCH --job-name=abc_euc
#SBATCH --array=1-20              # Run 20 parallel jobs
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load Python/3.11.4
python abc_non_dng.py ${SLURM_ARRAY_TASK_ID}
```
