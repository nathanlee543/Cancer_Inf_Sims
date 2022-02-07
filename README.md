# Code and simulations for "Inferring parameters of cancer evolution from sequencing and clinical data," Lee & Bozic.
 
## Description

This repository contains the code for the Monte Carlo simulations used to simulate the growth of tumors and accumulation of mutations (`simulations` directory). The code and input data used to generate estimates for CLL patients is in the `CLL` directory. PhylogicNDT outputs of the CCF posterior distributions as well as clone structures (the `.html` files) are included in the `data` directory.

## Getting Started

### Dependencies

* Monte Carlo simulations:
    * g++ compiler
    * recommend running on a server to take full advantage of parallelization. Was run on server with Ubuntu 18.04.1 LTS.
* analysis for figures 2 and 3
    * Jupyter Notebook
    * Python3
    * Python pacakges pandas, numpy, os, pathlib, sys, seaborn, matplotlib, math, and importlib.metadata
* CLL analysis
    * R
    * R packages dplyr, pracma, MASS, Rmpfr, ggplot2, tidyr


### Installing

```
git clone https://github.com/nathanlee543/Cancer_Inf_Sims
```
```
cd Cancer_Inf_Sims
```

## How to run

### Monte Carlo simulation
* Go to `simulation` directory (`cd simulation`)
* Change numeric values in the input file `run_inputs.txt`, or use the existing values, for:
    * b, birth rate of type-0 population
    * b1, birth rate of type-1 population
    * d, death rate of type-0 population
    * d1, death rate of type-1 population
    * t1, time to drive mutation
    * t, time between driver mutation and diagnosis
    * delta, time between two sequencing timepoints
    * u, mutation rate
    * num_threads, number of threads to parallelize over
    * runs, number of simulation runs
    * max_clones, maximum number of clones allowed
    * 0 for multithreading, 1 for just a single thread
* Compile and run: 
```
make parallel_paper
./MCSim run_inputs.txt
```

### Code accompanying figures 2 and 3
The jupyter notebook `simulations/figures_2_and_3.ipynb` computes estimates from the results of the Monte Carlo simulation, simulates sequencing reads, and comparing corrections for t<sub>1</sub> and γ.

### Generating estimates for CLL patients

* Go to the CLL directory (`cd CLL`)
* Run the `generate_CIs.R` script, e.g. from the command line,
```
Rscript generate_CIs.R
```

### Generating fishplots

* Run from `CLL` directory.
* To generate fishplots for each patient, run the `making_fish_diagrams.R` script, e.g. from the command line,
```
Rscript making_fish_diagrams.R
```



