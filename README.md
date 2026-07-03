# DAISY
The text is structured as follows: 
First, the output examples are briefly discussed. Then, the file structure is discussed. 
<br><br>
**A. Method of Morris**<br> 
In this case, +-20% of parameter default value, and WOFOST reference were considered to create bounds.
PBIAS, and RMSE of irrigation amount, yield amount, and soil moisture to 30 cm were taken as error functions to calculate elemntary effects. Different objective functions were taken as sensitivity indices themselves are sensitive to error functions. PBIAS and RMSE were chosen to represent error direction and magnitude respectively.

<figure>
  <img width="779" height="803" alt="image" src="https://github.com/user-attachments/assets/b4adbec5-4680-44b4-a3fe-b2920c3dc0d7" />
  <br>
  <figcaption><em>Fig.1. Covariance plot of elementary effects for Winter Wheat for RMSE of Yield. All parameters with both mean and standard deviation less than 0.2 were considered to be not sensitive.</em></figcaption>
</figure>

<figure>
  <img width="779" height="803" alt="image" src="https://github.com/user-attachments/assets/0207bc8e-bc2f-47bc-a1b3-97e6f5adcf29" />
  <br>
  <figcaption><em>Fig.2. Covariance plot of elementary effects for Winter Wheat for PBias of Yield.</em></figcaption>
</figure>
<br><br>

**B. eFAST** <br>
eFAST was used to check if the parameters sensitivities behaved differently in different weather classification (wet, dry, normal) for crop years. First, S1 (first order index, measuring direct contribution of parameter to output sensitivities) and ST (total order index, including interactions) indices for all parameters were calculated; parameters with S1 < 0.01 and ST < 0.05 for all time series were not considered to be sensitive. The parameters were ranked based on S1 index for each year, and the correlation between parameter sensitivities beteween all years was calculated. As seen in fig.3. parameter sensitivities behaved similarly across different weather types. 

<figure>
  <img width="722" height="590" alt="image" src="https://github.com/user-attachments/assets/9e882fcb-c481-4349-9da2-f07a48e0c2be" />
  <br>
  <figcaption><em>Fig.3. Kendall rank matrix for parameters screened from Morris sensitivity test for Winter Wheat based on Soil Moisture RMSE.</em></figcaption>
</figure>

The parameter similarity was counter intuitive to initial hypothesis; thus, this was checked by comparing parameter distribution across all weather types. Exclusive sets of behavioural parameter sets (rmse < 1Mg/ha) from the eFAST samples were obtained for all weather types. Then, for each parameter, normality test was done. Afterwards, if any parameter distribution was non-normal Kruskal Wallis test was done to check similarity between groups. An example of parameter distribution is shown in Fig.4.

The results showed similar parameter distribution for all weather types. Then, the model was calibrated for Winter Wheat yield using limits of acceptance (abs(PBIAS) < 10%) using DREAMzs algorithm in MATLAB. After multiple calibration runs, the model did not converge to a single parameter vector i.e. the calibration algorithm could not find a set of parameters equally suitable for all weather types. 
Due to this reason, despite the above conclusions, time-variant parameter was adopted for all crops. Likewise, the calibration was shifted to multiobjective calibration using NSGA-II.

<figure>
  <img width="736" height="775" alt="image" src="https://github.com/user-attachments/assets/5d02750a-72da-460b-90b2-2a4b17d79154" />
  <br>
  <figcaption><em>Fig.4. Comparison of distribution of parameter PenPar2 between behavioural parameter sets for different weather types. PenPar2 is the parameter determining the minimum temperature below which root growth will be 0. </em></figcaption>
</figure>
<br><br>

**C. NSGA-II calibration** <br>

<figure>
  <img width="1600" height="700" alt="image" src="https://github.com/user-attachments/assets/664e2529-c0a0-41c3-bb2a-0295e12fa520" />
  <br>
  <figcaption><em>Fig.5. Final pareto front of Winter Wheat after integrating 5 NSGA-II runs. The comparison is made with yield here. </em></figcaption>
</figure>

<figure>
  <img width="1600" height="700" alt="image" src="https://github.com/user-attachments/assets/ba054bc3-e33e-4534-a48b-a90ec51c7d33" />
  <br>
  <figcaption><em>Fig.6. Final pareto front of Winter Wheat after integrating 5 NSGA-II runs. The comparison is made with irrigation here. The bias in irrigation during validation period is possibly because the groundwater abstraction rule for irrigation was not implemented in the model. Reduced irrigation after 2018 drought was hence not factored in. </em></figcaption>
</figure>
<br><br>

**D. Auto-irrigation rule** <br>
Auto-irrigation rule obtained from the final pareto front were comparable to field regulation. 
<figure>
  <img width="640" height="480" alt="image" src="https://github.com/user-attachments/assets/99f867be-0a8e-436d-9dc5-27517de51f5b" />
  <br>
  <figcaption><em>Fig.7. Distribution of irrigation thresholds for different weather types. Irrigation threshold was calculated as % drop in net field capacity(nFK). For optimal irrigaiton, irrigation is triggered when nFK drops to 50%. </em></figcaption>
</figure>

<br><br>

**File structure** <br>
Example for WinterWheat. 

**1, Morris_parameter_generation_output_analysis**
  This file contains the script to create morris parameters. First, this requires an excel sheet or .txt file containing lower bound and uppper bound value for each parameter.  
  Also contains functions to compute and plot Elementary Effects.
  Also contains eFAST summarized and eFAST temporal analysis functions.
  Also contains differentiation of behavioural sets for different weather classes from eFAST parameters (as the number of parameter vectors generated was high enough from eFAST). 

**2, Morris_wwheat_exec_params**
  This contains the script to run all parameter vectors from 1. 
  Also contains script that compiles all outputs. This is used as input for further processing in 1. 

**4, NSGA_II_wwheat_eval_functions_pc**
  This contains all functions taken as helper functions for 5. 
  Here, year classification of crop is input manually. In the end of each years management file, the croptype can be changed. 
  This writes exec files for each parameter vector, reads their output, compiles all outputs, and sends the rmse matrix to 5 for analysis.

  
**5, NSGA_II_wwheat_exec_functions_pc**
  This executes the genetic algorithm. 
  For each generation, this determines the parameter vectors. These individuals are executed and measured using 4.
  For creating next generation, the output from 4 is used as input.
  
  After running for required number of generations or closing criteria, this reruns the pareto front and saves the outputs. 
  Then, this creates 3d plot, pareto parameter distribution plot from the pareto front. For each pareto individual also plots time series of yield (simulated, observed) and irrigaiton.

**6, Final_Outputs**
  The pareto front of 5 simulations are combined. The MATLAB code for this was provided. Then this plots the final pareto into an ensemble.


As exec files for crop are hardcoded, different files of the same structure were created for Potato and Sugar Beet. 
