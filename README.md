# LPDM_Hseih
Code to model the LSPDM footprint model based on Hsieh et al. 2000

## Setting up the environment
Activate environment using the yml file provided. 

## Usage
Pick the case for running the model for Unstable, Neutral or Stable conditions by modifying the variable caseA.
The parameters used for measurement height, roughness length, and Monin-Obukhov length can also be changed if desired.
The function that calculates fetch for the corresponding input parameters is located in the module module_fetch.py
The final parameters D and P are obtained which can be compared with the same parameters computed by the authors of Hsieh et al. 2000.
