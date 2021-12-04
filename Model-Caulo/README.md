## Code Structure for Model-Caulo
There are total of 8 script files, with the following dependecy structure:

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

### Usage
(For specific instructions on the file loaded, please refer to the comments in the code. All the output matrices are stored under folder `LoadMatrices`)

For running the Caulobacter cell cycle model `main5RegCPLX.m`, please ensure the correct optimized parameters is loaded `load('./LoadMatrices/opt_result.mat','para');`.
`events5RegCPLX.m` and `odes5RegCPLX.m` are the ODEs and events for the model, please make sure they are under the same folder as `main5RegCPLX.m`.

After running the Caulobacter cell cycle model `main5RegCPLX.m`, the resulted output will be stored under the name `model_result.mat` in `LoadMatrices` folder.

If you want to run the model for mutant cases, please change the conditions in `odes5RegCPLX.m` and save the outputs under different name as specified in `main5RegCPLX.m`.

For plotting experimental data vs. simulated data for regular proteins using `graphCellCycle.m` or `graphCellCycleFinal.m`, please ensure the correct outputs are loaded `load("./LoadMatrices/model_result.mat","tout","yout")`. `graphCellCycle.m` combines graphs together, and `graphCellCycleFinal.m` plots separately.

For plotting mutant cases in `mutantPlotterFinal.m`, please ensure the correct mutant case is loaded `load("./LoadMatrices/GcrA_Mutant.mat","tout","yout")` for gcrA mutant, similarly for different predictions. After loading the outputs for certain protein, please specify features below to see the outputed cell cycle:
```
% here to specify the mutant case to plot
% 1 gcra, 2 ccrm, 3 cdG, 4 PleD, 5 PdeA, 6 DnaA, 7 SM921, 8 ctrA3, 9 Prediction1, 10 Prediction2, 11 Prediction3, 12 Prediction4
mut = 1
```

For all plottings, please change code below to save the output plots.
```
% set to 1 to save figures
save_figs = 0;
```
