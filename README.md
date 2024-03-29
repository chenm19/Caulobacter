# Modeling the temporal dynamics of master regulators and CtrA proteolysis in *Caulobacter crescentus* cell cycle
*Caulobacter crescentus* cell undergoes polar morphogenesis and asymmetric cell division driven by precise interactions and regulation of proteins, making it an ideal model organism for studying bacterial cell development and differentiation. Here, we propose a comprehensive model to accurately characterize the mechanisms of cell cycle regulation that takes into account: a) chromosome replication; b) interactions between five master regulatory proteins including DnaA, GcrA, CcrM, CtrA, and SciP, as well as their corresponding mRNAs; c) cell cycle-dependent proteolysis of CtrA through hierarchical protease complexes. Our model is able to replicate the result of mounting observational studies and reproduce the main phenotype of seven mutant *Caulobacter crescentus* strains, and therefore be used as a tool to study the metabolism of proteins with similar behaviors.

## Structure of the repository
```
┌──MOPGA-Caulo
├──Model-Caulo
├──Model-Caulo-QSSA
└──VTMOP-Caulo
```

- `MOPGA-Caulo`: contains the scripts for running and plotting Multiple Objective Genetic Algorithm for optimizing parameters.
- `Model-Caulo`: contains the scripts for running and plotting the numerical simulation of Caulobacter cell cycle model for wild types and mutant cases.
- `Model-Caulo-QSSA`: contains the quasi-steady-state assumption of Complex 3 with the same set of parameters with Model-Caulo for running and plotting the numerical simulation of Caulobacter cell cycle model for wild types and mutant cases.
- `VTMOP-Caulo`: contains the scripts for running VTMOP Multiple Objective Optimization tool.

## Software Requirement
MATLAB, Fortran

## Contact
If you have any further questions or suggestions, please contact chenm@wfu.edu or rachelxu@vt.edu
