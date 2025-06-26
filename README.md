This repository contains the simulation and analysis scripts for the paper "Re-evaluating local genetic correlation analysis in the HDL framework", in response to "An enhanced framework for local genetic correlation analysis" by Li et al. ([link](https://www.nature.com/articles/s41588-025-02123-3)). 

The content is as follows, see comments inside each file for further details:  
- generate.r: main simulation script, generates summary statistics per genomic block
- generate_functions.r: helper functions for the simulation script
- analyse_lava.r: perform LAVA local genetic correlation analysis on the simulated summary statistics
- analyse_hdl.r: perform HDL-L local genetic correlation analysis on the simulated summary statistics
