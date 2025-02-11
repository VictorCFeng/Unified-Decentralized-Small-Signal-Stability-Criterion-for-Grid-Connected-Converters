# Unified Decentralized Small-Signal Stability Criterion for Grid-Connected Converters - Code

## Overview
This repository contains the simulation code for the paper *Unified Decentralized Small-Signal Stability Criterion for Grid-Connected Converters*. The simulations are built using Simulink in MATLAB 2024b. The repository is organized into three main folders:

### 1. New England 39 Bus System Electromagnetic Transient Simulation (`NewEngland39`)
This folder contains electromagnetic transient (EMT) simulation files for the New England 39-bus system. Each `.slx` Simulink file has a corresponding `.m` data file.

- Ensure that the system starts with small-signal stable parameters before switching to another parameter set to verify small-signal stability.
- Recommended startup method: High GFM active power damping and slow GFL operation.
- The simulation uses a discrete time step of `1e-5` (same as PSCAD's time step). The step size can be increased, but this may introduce errors.

### 2. Small-Signal Model Validation (`SmallSignalModel`)
This folder contains files for theoretical derivation, frequency scan simulations, and result validation:

- Run the frequency scan script `dqscan.m` first to obtain the actual impedance model of the inverter.
- Then, execute `xxx_compare.m` to compare the results and validate the small-signal model.

### 3. DW Calculation Functions (`DWShell`)
This folder includes functions for computing the DW shell and numerical range using sampling and rotation methods.

- The computation utilizes parallel processing (`parfor`). If parallel execution is not supported, replace `parfor` with `for`.
- Example usage cases are included in the folder.

## Requirements
- MATLAB 2024b
- Simulink
- Parallel Computing Toolbox (optional, for `parfor` functionality)

## Usage
1. Navigate to the respective folder based on the type of analysis.
2. Follow the instructions provided for each simulation task.
3. Adjust parameters as needed to explore system stability under different conditions.


