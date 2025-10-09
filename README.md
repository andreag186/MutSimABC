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
  - [ABC Inference on Empirical *E. melliodora* Data](#abc-inference-on-empirical-e-melliodora-data)
  - [Custom Tree Analysis](#custom-tree-analysis)
  - [Validation Framework](#validation-framework)
- [Understanding the Code](#understanding-the-code)
- [Defining Custom Tree Topologies](#defining-custom-tree-topologies)
- [Computational Requirements](#computational-requirements)
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
python abc_validation.py 0
```
This will:

- Run ABC-Reject validation for sample 1 from abc_validation_samples.csv
- Sample parameters (StD, biasVar, input_mut) from prior distributions
- Simulate mutation distributions using the tree topology defined for that sample
- Compare simulated and observed mutation distributions, accepting parameter sets where the Euclidean distance < ε = 20
- Continue sampling (up to ~5,000 trials) until 100 accepted samples are obtained
- Save results to:
result_valid/accepted_samples_sample1_epsilon20_<timestamp>.csv

To validate a different sample, change the task ID argument (e.g., `python abc_validation.py 42` for sample 43- there are 200 samples total (index starting at 0).

## Repository Structure
```
MutSimABC/
├── README.md                      # This file
├── requirements.txt               # Python dependencies
│
├── abc_non_dng.py                # ABC inference for pre-DNG E. melliodora data
├── abc_phylo.py                  # ABC inference for post-DNG E. melliodora data
├── abc_custom_tree.py            # Template for custom tree analysis
├── abc_validation.py             # Validation framework with known ground truth
│
├── abc_validation_samples_combined.csv    # Pre-generated validation samples (169 samples)
│
├── results/                      # Output directory for pre-DNG accepted samples
├── result_phylo/                 # Output directory for post-DNG accepted samples
├── results_my_tree/              # Output directory for custom tree analysis
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

For parallel execution on computing clusters 
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
### Custom Tree Analysis
To analyze your own tree data, use `abc_custom_tree.py` as a template:

#### Step 1: Define Your Tree Topology

Locate `MODIFY SECTION 1` in `abc_custom_tree.py` and add your tree structure:
```python
tree_topologies_dict = {
    "my_custom_tree": {
        "numBranch": 6,           # How many terminal branches you sampled
        "age": 150,               # Maximum root-to-tip distance (years)
        "s10": 10,                # Trunk before first split
        
        # Right side branches (see figure in "Defining Custom Tree Topologies")
        "b11": 40, "bb11": 25,    # Internal + terminal
        "b12": 40, "bb12": 25,
        "b13": 25,                # Last internal = terminal
        
        # Left side branches
        "s40": 40, "b41": 25,
        "s41": 40, "b42": 25,
        "s42": 25
    }
}
```
*See [Defining Custom Tree Topologies](#defining-custom-tree-topologies) for detailed nomenclature rules.*

#### Step 2: Set Your Tree in Simulation
In `MODIFY SECTION 2`, change the tree name:
```python
tree_list, tree_dict, numBranch, age = create_tree_list_and_dict(
    tree_topologies_dict['my_custom_tree']  # ← Use your tree name
)
```
#### Step 3: Add Your Observed Mutation Data
In MODIFY SECTION 3, replace with your empirical mutation counts per branch:
```python
def calculate_distance(unique_mutations):
    real_distribution = [
        45.2,    # Unique mutations in branch 1 (bb11)
        38.7,    # Branch 2 (bb12)
        52.1,    # Branch 3 (b13)
        41.5,    # Branch 4 (b41)
        35.8,    # Branch 5 (b42)
        48.3     # Branch 6 (s42)
    ]  # ← YOUR DATA HERE - order must match tree topology
```
**Important: Branch order must match the tree structure (right branches first, then left branches).**

### Step 4: Adjust Priors
Using prior knowledge of expected mutation rates of your long-lived trees, you can change the input_mut prior change. We reccomend a sensitivity analysis. Reffering to the methods/results you may wish to include or remove StD= 0 due to model-jumping effects. Adjust the priors by finding `MODIFY SECTION 4`:
```python
def sample_prior():
    StD = np.random.uniform(1, 5)              # Keep as default
    biasVar = np.random.uniform(0.5, 10)       # Keep as default
    input_mut = np.random.uniform(5e-11, 2e-9) # Narrow if you have expectations
    return StD, biasVar, input_mut
```

#### Step 5: Run Analysis 
```python
python abc_custom_tree.py
```
Results will be saved to `results_my_tree/accepted_samples_*.csv`.

#### Step 6: Analyze Posterior
Use ArviZ to compute summary statistics and visualize the posterior distribution of accepted parameter sets:
```python
import pandas as pd
import arviz as az

# Load accepted samples
df = pd.read_csv('results_my_tree/accepted_samples_default_job_default_batch_*.csv')

# Create ArviZ inference data object
idata = az.from_dict(posterior={
    'mu': df['input_mut'].values,
    'StD': df['StD'].values,
    'biasVar': df['biasVar'].values
})

# Summary statistics with 95% HPD intervals
print(az.summary(idata, hdi_prob=0.95))

# Plot posteriors
az.plot_posterior(idata)
```
What this does:

- `az.summary()` computes posterior mean, standard deviation, and 95% highest posterior density (HPD) intervals for each parameter. The HPD interval represents the range containing 95% of the posterior probability mass.
- `az.plot_posterior()` generates histograms showing the posterior distribution shape for each parameter, with HPD intervals marked.

This basic analysis provides point estimates (mean) and uncertainty quantification (HPD intervals) for the mutation rate and developmental parameters. For more advanced diagnostics (effective sample size, convergence checks), see the ArviZ documentation.

### Validation Framework

The validation framework tests the accuracy of MutSimABC on synthetic simulated samples with known 'true' observed mutations.

**Pre-Generated Validation Samples**

`abc_validation_samples_combined.csv` contains 200 pre-generated validation samples with columns:
- `Sample_ID`: Unique identifier (0-199)
- `input_mut`: True mutation rate
- `StD`: True elongation parameter
- `biasVar`: True branching bias
- `tree_topology`: Tree architecture code (e.g., "bS4", "ubL8")
- `Unique_Mutations`: true mutation distribution across branches

**Running Validation**

Single sample:
```bash
python abc_validation.py 0  # Runs validation for sample 1
```
Batch Validation with SLURM:
```bash
#!/bin/bash
#SBATCH --job-name=abc_valid
#SBATCH --array=0-199           # Run all 200 validation samples
#SBATCH --time=00-48:00:00 
#SBATCH --cpus-per-task=20  
#SBATCH --mem=100GB

module load Python/3.11.4
python abc_validation.py ${SLURM_ARRAY_TASK_ID}
```
Each validation run continues until 100 accepted samples are collected, ensuring sufficient posterior sampling for HPD interval calculation.

**Analysing Validation Results**

After running the validation, analyse the coverage and accuracy as below:
```python
# SET UP
import pandas as pd
import glob
import arviz as az

# LOAD ALL VALID RESULTS
result_files = glob.glob('result_valid/accepted_samples_sample*_epsilon20_*.csv')
validation_results = []

for file in result_files:
    df = pd.read_csv(file)
    sample_id = int(file.split('sample')[1].split('_')[0])
    
    # EXRTRACT TRUE VALUES
    true_values = pd.read_csv('abc_validation_samples_combined.csv').iloc[sample_id]
    
    # CALC 95% HPD INT
    idata = az.from_dict(posterior={
        'mu': df['input_mut'].values,
        'StD': df['StD'].values,
        'biasVar': df['biasVar'].values
    })
    summary = az.summary(idata, hdi_prob=0.95)
    
    # CHECK IF TRUE VALUE INSIDE HPD
    mu_in_hpd = (true_values['input_mut'] >= summary.loc['mu', 'hdi_2.5%']) & \
                (true_values['input_mut'] <= summary.loc['mu', 'hdi_97.5%'])
    StD_in_hpd = (true_values['StD'] >= summary.loc['StD', 'hdi_2.5%']) & \
                 (true_values['StD'] <= summary.loc['StD', 'hdi_97.5%'])
    biasVar_in_hpd = (true_values['biasVar'] >= summary.loc['biasVar', 'hdi_2.5%']) & \
                     (true_values['biasVar'] <= summary.loc['biasVar', 'hdi_97.5%'])
    
    validation_results.append({
        'sample_id': sample_id,
        'mu_in_hpd': mu_in_hpd,
        'StD_in_hpd': StD_in_hpd,
        'biasVar_in_hpd': biasVar_in_hpd
    })

results_df = pd.DataFrame(validation_results) 
print(f"Coverage: μ={results_df['mu_in_hpd'].mean():.1%}, "
      f"StD={results_df['StD_in_hpd'].mean():.1%}, "
      f"σ={results_df['biasVar_in_hpd'].mean():.1%}")
```
The expected coverage across samples should be approximately 95% for μ, StD and σ . A more in-depth statistical analysis can be conducted via `results_df` & `validation_results`.

## Understanding the Code

### Core Simulation Functions (from Tomimoto & Satake 2023)

These functions simulate meristem dynamics and mutation accumulation:

1. `mutInStemCells(num_stem, t, mu_0, Genom, st_d)`
  - Simulates elongation in the main stem before first branch split
  - Parameters: 5 stem cells, time t, mutation rate, genome size, StD parameter
  - Returns: mutation history of stem cells

2. `mutInBrStemCells(num_stem, stemCells, t, mu_0, Genom, st_d)`
  - Simulates elongation along a branch
  - Takes initial stem cells from parent meristem
  - Returns: mutation history along branch

3. `sample_mutations(stem_cells, num_stem, num_div, mu_0, Genom, bias_d)`
  - Simulates branching event with spatial bias (`bias_d`)
  - Uses wrapped normal distribution to model cell sampling
  - Returns: stem cells for new axillary meristem

4. `simulate_somatic_mutations(tree_dict, NumStem, NumTime, mu_0, GenSize, StD, NumDiv, nDiv, biasVar)`
  - Main simulation orchestrator
  - Builds mutation matrix across entire tree
  - Returns: binary matrix of mutations × branches

### Tree Topology Handling 

`create_tree_list_and_dict(tree_topology)`
- Parses tree topology dictionary
- Extracts branch ages and connectivity
- Returns: processed tree structure for simulation

`tree_topologies_dict`
- Dictionary containing predefined tree architectures
- Keys: tree codes (e.g., "bS4", "ubL8", "test")
- Values: nested dictionaries with branch specifications (age/num left or right)

### ABC Framework Functions
1. `sample_prior()`
  - Samples parameters from uniform prior distributions
  - Returns: (StD, biasVar, input_mut)

2. `run_simulation(StD, biasVar, input_mut)`
  - Runs full simulation pipeline with sampled parameters
  - Calls `simulate_somatic_mutations()` → `mut_dist_func()` → `gen_matrices()` → `calc_variants()`
  - Returns: unique mutations per branch, estimated mutation rate

3. `calculate_distance(unique_mutations)`
  - Computes Euclidean distance between simulated and observed distributions
  - Returns: scalar distance value

4. `abc_rejection(num_samples, epsilon, job_id, batch_id)`
  - Main ABC loop
  - Samples → simulates → compares → accepts/rejects
  - Saves accepted samples every 100 iterations

### Mutation Analysis Functions
1. `mut_dist_func(resultList_2)`
  - Analyzes mutation patterns across branches
  - Calculates shared, unique, and total mutations
  - Returns: mutation distribution statistics

2. `gen_matrices(tree_dict, unique_mutations, shared_mutations, mutations_counts)`
  - Constructs genetic and physical distance matrices
  - Returns: pairwise distance matrices

3. `calc_variants(genetic_matrix, physical_matrix, tree_dict, GenSize)`
  - Performs linear regression of genetic vs. physical distance
  - Estimates mutation rate per site per year
  - Returns: variant count, regression equation, mutation rate

## Defining Custom Tree Topologies

### Understanding the Nomeclature
The tree topology system follows the nomenclature from Tomimoto & Satake (2023), originally developed for *Populus trichocarpa*.

<img width="450" height="600" alt="image" src="https://github.com/user-attachments/assets/58009964-f507-47e8-b0bb-89566d2b96dc" />

*Tree topology nomenclature system. (a) Wild P. trichocarpa trees used in Hoffmeister et al. (2020). (b) Branch ages at terminals (black) and junctions (grey). (c) Coding system: right branches (b11-bxx, bb11-bbxx) and left branches (s40-sxx, b41-bxx). Figure adapted from Hoffmeister et al. (2020) and Tomimoto & Satake (2023).*

**Key Points:**

- Ages represent years (10 cm/year growth assumed)
- Right branches: b11-bxx (internal), bb11-bbxx (terminal)
- Left branches: s40-sxx (internal), b41-bxx (terminal)
- Last internal branch without corresponding terminal = treated as terminal


### Example: Balanced 4-Branch Tree
```python
"bS4": {
    "numBranch": 4,
    "age": 123,
    "s10": 10,
    "b11": 75, "bb11": 38, "b12": 38,  # Right: 1 internal, 2 terminals
    "s40": 75, "b41": 38, "s41": 38    # Left: 1 internal, 2 terminals
}
#### **(bS4)** Balanced Tree 4 Terminal Branches (terminal branches < internal branches)
```
Visual Structure
```
          ┌─bb11 (38y)
    ┌─b11 (75y)
    │    └─b12 (38y)
s10 ┤
(10y)│    ┌─b41 (38y)
    └─s40 (75y)
         └─s41 (38y)
```

To test your tree has been added correctly, run the following code:
```python
from abc_custom_tree import create_tree_list_and_dict, tree_topologies_dict

tree_list, tree_dict, numBranch, age = create_tree_list_and_dict(
    tree_topologies_dict['my_custom_tree']
)

print(f"Branches: {numBranch}, Max age: {age}")
print(f"Right: {len(tree_dict['Rt'])} terminals")
print(f"Left: {len(tree_dict['Lt'])} terminals")
# Total terminals should equal numBranch
```

## Computational Requirements
**Per ABC Run (5,000 trials)**

Runtime: ~6-12 hours 
Memory: 8-16 GB RAM
CPU: Single core (parallelizable via SLURM array jobs)
Storage: ~50-100 MB per accepted sample CSV

**Validation Suite (200 samples)**

Runtime: ~1-2 weeks on HPC cluster with 50 parallel jobs
Memory: 16 GB per job
Storage: ~20 GB total for all validation results

**Tested Environments**

- NeSI (New Zealand eScience Infrastructure) - Mahuika cluster
- Local workstation - Ubuntu 22.04, 32GB RAM, 16-core CPU

## Citation

If you use MutSimABC in your research, please cite:
```bibtex
@article{grecu2025mutsimabc,
  title={MutSimABC: Likelihood-Free Inference of Somatic Mutation Rates Accounting for Meristem Dynamics in Long-Lived Trees},
  author={Grecu, Andrea Maria and [Co-Authors]},
  journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  year={2025},
  note={In review}
}
```

Please also cite the original Tomimoto & Satake (2023) meristem model:
```bibtex
@article{Tomimoto_2023,
title={Modelling somatic mutation accumulation and expansion in a long-lived tree with hierarchical modular architecture},
volume={565}, ISSN={0022-5193},
url={http://dx.doi.org/10.1016/j.jtbi.2023.111465},
DOI={10.1016/j.jtbi.2023.111465},
journal={Journal of Theoretical Biology},
publisher={Elsevier BV},
author={Tomimoto, Sou and Satake, Akiko},
year={2023},
month=may,
pages={111465} }
```

And the empirical data source:
```bibtex
@article{orr2020population,
  title={A phylogenomic approach reveals a low somatic mutation rate in a long-lived plant},
  author={Orr, Adam J and Padovan, Amanda and Kainer, David and Külheim, Carsten and Bromham, Lindell and Bustos-Segura, Carlos and Foley, William and Haff, Tonya and Hsieh, Ji-Fan and Morales-Suarez, Alejandro and others},
  journal={Proceedings of the Royal Society B},
  volume={287},
  number={1922},
  pages={20193364},
  year={2020},
  publisher={The Royal Society}
}
```

## Contact
### **Andrea Maria Grecu**
MSc Graduate, University of Auckland

Email: amgstar86@gmail.com
