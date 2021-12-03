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
    ├──graphCellCycle.m / graphCellCycleFinal.m

    
- `caulo_MOPGA.m`: script for running Genetic Algorithm for parameter optimization
- `fitness.m`: fitness function used by Genetic Algorithm for parameter optimization




- `pareto_plot.m / pareto_plotFinal.m`: script for plotting 




## Contact
If you have any further questions or suggestions, please contact chenm@wfu.edu
