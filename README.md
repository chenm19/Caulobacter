# Modeling the temporal dynamics of master regulators and CtrA proteolysis in *Caulobacter crescentus* cell cycle

Matlab code for simulating Caulobacter cell cycle model 
Run main5RegCPLX.m

## Code Structure for Model-Caulo
There are total of 9 script files, with the following dependecy structure

    main5RegCPLX.m
    ├──load_para.m
    ├──load_data.m
    ├──events5RegCPLX.m
    ├──odes5RegCPLX.m
    ├──mutantPlotterFinal.m
    └──graphCellCycle.m / graphCellCycleFinal.m
    
- `main5RegCPLX.m`: numericaly simulation of Caulobacter cell cycle model.
- `events5RegCPLX.m`: script for detecting specific events of the cell cycle.
- `odes5RegCPLX.m`: script for ordinary differential equations of Caulobacter cell cycle model
- `load_para.m` and `load_data.m`: load the data and parameters used in the model
- `mutantPlotterFinal.m` and `graphCellCycle.m / graphCellCycleFinal.m`: plotting simulations for mutant cases as well as regular proteins

### Usage:
For running the Caulobacter cell cycle model `main5RegCPLX.m`, please load `opt_result.mat` to load the optimized parameters.
`events5RegCPLX.m` and `odes5RegCPLX.m` are the ODEs and events for the model, please make sure they are under the same folder as `main5RegCPLX.m`.

## Code Structure for MOPGA-Caulo

    ├──load_para.m
    ├──load_data.m
    ├──main5RegCPLX.m
    |  ├──events5RegCPLX.m
    |  ├──odes5RegCPLX.m
    ├──caulo_MOPGA.m
    |  ├──fitness.m
    ├──mutantPlotterFinal.m
    └──graphCellCycle.m / graphCellCycleFinal.m
    
- `caulo_MOPGA.m`: script for running Genetic Algorithm for parameter optimization.
- `fitness.m`: fitness function used by Genetic Algorithm for parameter optimization.
- `pareto_plot.m / pareto_plotFinal.m`: script for plotting pareto front and calculate the sensitivity and goodness.

### Usage:
For running parameter optimization in `caulo_MOPGA.m`, please make sure that `events5RegCPLX.m`, `odes5RegCPLX.m`, and `fitness.m` are under the same folder as `caulo_MOPGA.m`.

Load different starting points or change options in below for better Multi-obj GA performance:
```
options = optimoptions('gamultiobj','MaxGenerations', 30,'PopulationSize', 50, 'InitialPopulation', para, 'PlotFcn', @gaplotpareto);
```

For plotting experimental data vs. simulated data for proteins, please load `model_result.mat` for outputs.

For plotting mutant cases in `mutantPlotterFinal.m`, please load the outputed time and value for each mutant cases first (for example, loading `CcrM_Mutant.mat` if you want to plotting mutant case for CcrM), similarly for different predictions.
After loading the outputs for certain protein, please specify features below to see the outputed cell cycle:
```
% set to 1 to save figures
save_figs = 0;

% here to specify the mutant case to plot
% 1 gcra, 2 ccrm, 3 cdG, 4 PleD, 5 PdeA, 6 DnaA, 7 SM921, 8 ctrA3, 9 Prediction1, 10 Prediction2, 11 Prediction3, 12 Prediction4
mut = 1
```

### VTMOP-Caulo


## Software Requirement
MATLAB

## Citations
Please refer to our paper-------

## Contact
If you have any further questions or suggestions, please contact chenm@wfu.edu
