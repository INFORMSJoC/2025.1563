# Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules

This repository contains the source code and data associated with the paper

**Milad Hasanzadeh and Amin Kargarian**,  
"Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules,"  
*INFORMS Journal on Computing*.

## Cite

The citation information for the final INFORMS Journal on Computing code and data repository will be added after the repository DOI is assigned.

## Description

This repository contains the computational implementation used in the paper.

The code implements distributed DC and AC optimal power flow using the
Alternating Direction Method of Multipliers (ADMM), together with the proposed
asymptotic, finite-time, and fixed-time update rules for the primal and dual
variables.

The numerical experiments compare the proposed update rules with classical
ADMM on benchmark power-system test cases.

## Repository Structure

The repository is organized as follows:

- `src/`: source code used to run the distributed DC and AC OPF experiments.
- `data/`: input data and test-system files required by the source code.
- `AUTHORS`: author information.
- `LICENSE`: software license.

## Software Requirements

The source code requires the software and packages used in the computational
experiments.

Please install the required dependencies before running the code.

The exact software versions used in the paper are:

- MATLAB
- YALMIP
- IPOPT

## Data

The `data/` directory contains the input files used in the computational
experiments.

The experiments consider the benchmark test systems reported in the paper,
including:

- 48-bus system
- 72-bus system
- 118-bus system
- 2383wp-bus system
- 2869pegase-bus system
- 9241pegase-bus system
- 24464goc-bus system


## Source Code

The `src/` directory contains the MATLAB source code used to reproduce the
computational experiments reported in the paper.

## Running the Code

1. Install the required software and dependencies.

2. Download or clone this repository.

3. Open MATLAB.

4. Set the repository directory as the MATLAB working directory.

5. Add the source directory to the MATLAB path:

   ```matlab
   addpath(genpath('src'))