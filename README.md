[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper
**Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules**
by Milad Hasanzadeh and Amin Kargarian.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

Paper DOI: To be added when assigned by INFORMS Journal on Computing.

https://doi.org/10.1287/ijoc.2025.1563.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{HasanzadehKargarian2026,
  author =        {Milad Hasanzadeh and Amin Kargarian},
  publisher =     {INFORMS Journal on Computing},
  title =         {Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules},
  year =          {2026},
  doi =           {10.1287/ijoc.2025.1563.cd},
  url =           {https://github.com/INFORMSJoC/2025.1563},
  note =          {Available for download at https://github.com/INFORMSJoC/2025.1563},
}
```

## Description

This repository contains the MATLAB source code and data used in the
computational experiments reported in the paper.

The code implements distributed DC and AC optimal power flow using the
Alternating Direction Method of Multipliers (ADMM), together with the proposed
asymptotic, finite-time, and fixed-time update rules for the primal and dual
variables.

The numerical experiments compare the proposed update rules with classical
ADMM on benchmark power-system test cases.

The proposed computational methods include classical ADMM, asymptotic
primal-variable update rules, finite-time primal-variable update rules,
fixed-time primal-variable update rules, finite-time dual-variable update
rules, and fixed-time dual-variable update rules.

The methods are evaluated for both distributed DC OPF and distributed AC OPF.

The repository is organized as follows:

- `src/`: MATLAB source code used to run the distributed DC and AC OPF experiments.
- `data/`: input data and test-system files required by the source code.
- `AUTHORS`: author and contact information.
- `LICENSE`: software license.

The computational experiments use MATLAB, YALMIP, and IPOPT.

The `data/` directory contains the input files used in the computational
experiments. The benchmark test systems reported in the paper include the
48-bus, 72-bus, 118-bus, 2383wp-bus, 2869pegase-bus, 9241pegase-bus, and
24464goc-bus systems.

The `src/` directory contains the MATLAB source code used to reproduce the
computational experiments reported in the paper.

## Building

The source code is implemented in MATLAB and does not require compilation.

Before running the code, install MATLAB, YALMIP, IPOPT, and all other required
dependencies.

After downloading or cloning the repository, open MATLAB and set the repository
root directory as the current working directory.

Add the source directory to the MATLAB path using

```
addpath(genpath('src'));
```

Make sure the required test-system files are available in the `data` directory.

Run the appropriate MATLAB main file in the `src` directory for the desired
test system and OPF formulation. Use the classical, asymptotic, finite-time, or
fixed-time ADMM configuration corresponding to the computational experiment
being reproduced.

## Results

The source code reproduces the computational experiments reported in the paper.

The numerical experiments compare the classical ADMM method with the proposed
asymptotic, finite-time, and fixed-time primal- and dual-variable update rules
for both distributed DC and AC OPF.

The reported computational results include ADMM iteration counts,
inter-regional primal residuals, optimality gaps, computation times, and
convergence trajectories.

The experiments use the following benchmark test systems:

- 48-bus system
- 72-bus system
- 118-bus system
- 2383wp-bus system
- 2869pegase-bus system
- 9241pegase-bus system
- 24464goc-bus system

## Replicating

To replicate the computational experiments reported in the paper:

1. Install MATLAB, YALMIP, IPOPT, and all required dependencies.
2. Clone or download this repository.
3. Open MATLAB and set the repository root directory as the working directory.
4. Add the `src` directory to the MATLAB path:

```
addpath(genpath('src'));
```

5. Verify that the required test-system files are available in the `data`
   directory.
6. Run the corresponding MATLAB main file in the `src` directory for the desired
   DC or AC OPF test case.
7. Use the classical, asymptotic, finite-time, or fixed-time ADMM configuration
   corresponding to the experiment being reproduced.

The resulting outputs can be used to reproduce the iteration counts, convergence
behavior, optimality gaps, computation times, and convergence trajectories
reported in the paper.

## Support

For support in using this software, please contact the authors listed in the
`AUTHORS` file.
