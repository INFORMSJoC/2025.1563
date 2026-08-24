# Data

This directory contains the benchmark power-system data used in the
computational experiments reported in the paper

**Alternating Direction Method of Multipliers-Based Distributed Optimal Power Flow with Fixed- and Finite-Time Update Rules**

by Milad Hasanzadeh and Amin Kargarian.

## Data Source

The benchmark data used in this repository are based on DPLib:

https://github.com/LSU-RAISE-LAB/DPLib

DPLib is a MATLAB-based benchmark library developed by the LSU RAISE Lab for
distributed power-system analysis and optimization.

Some of the benchmark datasets used in this study were obtained directly from
the benchmark cases provided in DPLib. The remaining datasets were generated
using the power-system partitioning tool provided by DPLib.

The partitioning tool converts centralized benchmark power-system cases into
multi-region datasets suitable for distributed optimization and distributed
optimal power flow studies.

The benchmark systems used in the computational experiments reported in the
paper include:

- 48-bus system
- 72-bus system
- 118-bus system
- 2383wp-bus system
- 2869pegase-bus system
- 9241pegase-bus system
- 24464goc-bus system

## File Format

All benchmark datasets in this directory are stored as MATLAB `.mat` files.

Each `.mat` file contains the power-system and regional-decomposition
information required to formulate the corresponding distributed power-system
optimization problem.

The internal data structures and variables follow the benchmark-data format
used in DPLib.

## Data Contents

Depending on the benchmark system, each `.mat` file contains information
describing the original power network and its multi-region decomposition.

This information includes, as applicable:

- bus data;
- generator data;
- branch and transmission-line data;
- active and reactive power demands;
- generator operating limits;
- voltage-magnitude and voltage-angle information;
- line-flow limits;
- generator cost information;
- assignments of buses to regions;
- assignments of generators to regions;
- internal transmission lines for each region;
- inter-regional tie-lines;
- boundary buses;
- neighboring-region relationships; and
- other information required to construct the distributed DC and AC OPF
  formulations.

The exact organization and meaning of the variables stored in the `.mat` files
follow the definitions and data structures described in DPLib.

For additional details on the contents and organization of the benchmark files,
please refer to:

https://github.com/LSU-RAISE-LAB/DPLib

## Regional Decomposition

The distributed OPF formulation requires each benchmark power system to be
partitioned into multiple interconnected regions.

For the datasets generated using the DPLib partitioning tool, the tool was used
to create the regional decomposition of the corresponding centralized benchmark
system.

The resulting `.mat` files contain the regional information needed by the source
code in this repository, including local network components, inter-regional
tie-lines, boundary buses, and neighboring-region connections.

Some of the partitioned benchmark datasets were already available in DPLib,
while others were generated using the same DPLib partitioning procedure for the
computational experiments reported in the paper.

## Use in This Repository

The MATLAB programs in the `src/` directory load the corresponding `.mat`
benchmark file from this directory and use its network and regional information
to construct the distributed DC or AC OPF problem.

The same benchmark data are used to evaluate:

- classical ADMM;
- asymptotic primal-variable update rules;
- finite-time primal-variable update rules;
- fixed-time primal-variable update rules;
- finite-time dual-variable update rules; and
- fixed-time dual-variable update rules.

## Data Generation and Reproducibility

The benchmark files that were not taken directly from the existing DPLib
datasets were generated using the partitioning tool provided by DPLib.

The partitioning methodology and corresponding implementation are available in
the DPLib repository:

https://github.com/LSU-RAISE-LAB/DPLib

Therefore, the regional benchmark datasets can be reproduced from the
corresponding centralized benchmark cases using the DPLib partitioning tool.
