# DAISY
Example for WinterWheat. 

File structure: 
1, Morris_parameter_generation_output_analysis
  This file contains the script to create morris parameters. First, this requires an excel sheet or .txt file containing lower bound and uppper bound value for each parameter. 
  In this case, +-10%, and WOFOST reference are considered. 

  Also contains functions to compute morris indices.

  Also contains eFAST summarized and eFAST temporal analysis functions. This part not used in final output.

  Also contains differentiation of behavioural sets for different weather classes from eFAST parameters (as the number of parameter vectors generated was high enough from eFAST). 
    This was also not used in final output. 

2, Morris_wwheat_exec_params
  This contains the script to run all parameter vectors from 1. 
  Also contains script that compiles all outputs. This is used as input for further processing in 1. 

4, NSGA_II_wwheat_eval_functions_pc
  This contains all functions taken as helper functions for 5. 
  Here, year classification of crop is input manually. In the end of each years management file, the croptype can be changed. 
  This writes exec files for each parameter vector, reads their output, compiles all outputs, and sends the rmse matrix to 5 for analysis.


  
5, NSGA_II_wwheat_exec_functions_pc
  This executes the genetic algorithm. 
  For each generation, this determines the parameter vectors. These individuals are executed and measured using 4.
  For creating next generation, the output from 4 is used as input.
  
  After running for required number of generations or closing criteria, this reruns the pareto front and saves the outputs. 
  Then, this creates 3d plot, pareto parameter distribution plot from the pareto front. For each pareto individual also plots time series of yield (simulated, observed) and irrigaiton.

6, Final_Outputs
  The pareto front of 5 simulations are combined. The MATLAB code for this was provided. Then this plots the final pareto into an ensemble.


As exec files for crop are hardcoded, different files of the same structure were created for Potato and Sugar Beet. 
