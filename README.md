# lhcMR-simulations

:grey\_exclamation: Once these scripts get integrated into the main `lhcMR` package (found [here](https://github.com/LizaDarrous/lhcMR)), this repo will be removed. :grey\_exclamation:

The following scripts require two data files to be present in the working directory:
- the correlation matrix file generated from chromosome 10 (referenced in [Supplementary Methods 1.3](https://www.nature.com/articles/s41467-021-26970-w#Sec18) of *Darrous et al.*) can be downloaded from [here](https://drive.google.com/file/d/17Q51_kHFWL6tpEWPryEn9vr-2Cki0P8r/view?usp=sharing).
- the `pi1` and `sig1` columns corresponding to the proportion and variance of the effective SNPs in the local LD of each genotyped SNP also referenced in [Supplementary Methods 1.3](https://www.nature.com/articles/s41467-021-26970-w#Sec18), and which can be downloaded [here](https://drive.google.com/file/d/1ZToHtACcbOfww5d5QAbG_VITzjhqZFR1/view?usp=sharing). 

**Note** that these files are quiet large, and the simulation scripts will require a lot of space on a server to run (master node). 

There are two scripts involved in the LHC-MR simulations:

#### 1 - LHC-MR_SimulationGeneration
After setting the working directory, the scenario number and the number of data generations needed, 3 booleans must be set to further aid the generation:
- `detailedSave` will save a larger `.Rdata` file containing direct effects on X,U and Y
- `kurtosi`s will enable a non-normal distribution / skewed distribution with kurtosis set up at a later if-clause for the direct effects
- `twoUs` will enable the generation of two underlying confounders, the parameters of which will be set up at a later if-clause

Note that both `kurtosis` and `twoUs` can't be TRUE at the same time (haven't been tested).  
Afterwards all the parameters will have to be updated to match your desired scenario, e.g.: no reverse causal effect (`bet=0`), no causal effect (`alp=0`), negative confounder effect (`qX/qY = -ve`), etc...

This script will generate `nsimi` data generations for this scenario, which will be read into the second script for optimisation and parameter estimation.

#### 2 - LHC-MR_SimulationOptimisation
After setting the working directory, the scenario number, the number of data generations created from the previous script, and the number of starting points that *each* data generation will optimise `nSP`, 2 parameters should be set:
- `saveSP`: boolean that will save the starting points used for each data generation
- `paral_method`: either "rlurm" or "lapply" to choose how the optimisation can be paralleled over the number of starting points.

The script will then load each generated data simulation, run a 2-step analysis where the pi, and i are first optimised, followed by the optimisation of the rest of the parameters.

The outputs include a detailed optimisation for each data generation, as well as a summary of the results for all the data generations where the best parameter estimates are stored.


