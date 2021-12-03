## Code Structure for Model-Caulo
There are total of 7 script files, with the following dependecy structure

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

Load the optimized parameter output (`opt_result.mat`), run `main5RegCPLX.m` to load the model, and the resulted output will be stored under the name `model_result.mat`. `events5RegCPLX.m` and `odes5RegCPLX.m` are the ODEs and events for the model, please make sure they are under the same folder as `main5RegCPLX.m`.

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
