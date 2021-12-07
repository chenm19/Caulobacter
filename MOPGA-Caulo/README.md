# MOPGA-Caulo
This folder contains all the code used to run Multiple Objective Genetic Algorithm for parameter optimization.

## Code Structure for MOPGA-Caulo
There are total of 9 script files, with the following dependecy structure

    caulo_MOPGA.m
    ├──load_para.m
    ├──load_data.m
    ├──fitness.m
    |  ├──events5RegCPLX.m
    |  ├──odes5RegCPLX.m
    └──pareto_plot.m
    
- `caulo_MOPGA.m`: script for running Genetic Algorithm for parameter optimization.
- `events5RegCPLX.m`: script for detecting specific events of the cell cycle.
- `odes5RegCPLX.m`: script for ordinary differential equations of Caulobacter cell cycle model
- `load_para.m` and `load_data.m`: contains the data and parameters used in the model
- `fitness.m`: fitness function used by Genetic Algorithm for parameter optimization.

### Usage:
For running parameter optimization in `caulo_MOPGA.m`, please make sure that `events5RegCPLX.m`, `odes5RegCPLX.m`, and `fitness.m` are under the same folder as `caulo_MOPGA.m`.

Load different starting points or change options in below for better Multi-obj GA performance:
```
% load good starting point
% load('best4S.mat')
% para = val(2,:);

options = optimoptions('gamultiobj','MaxGenerations', 30,'PopulationSize', 50, 'InitialPopulation', para, 'PlotFcn', @gaplotpareto);
```

The outputed parameters will be stored in `opt_result.mat` in `LoadMatrices` folder.

(For details of `events5RegCPLX.m`, `odes5RegCPLX.m`, please refer to `Model-Caulo`)
