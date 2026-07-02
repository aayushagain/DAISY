# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 18:16:08 2025

@author: Adhikari
"""

import sys
sys.path.append(r'C:\Users\Adhikari\.spyder-py3')
# sys.path.append(r'D:\Aayush\WinterWheat NSGA\nsgaII_wwheat')
sys.path.append(r'C:\ProgramData\spyder-6')
sys.path.append(r'C:\ProgramData\spyder-6\Lib\site-packages\pymoo')
from nsga_II_wwheat_eval_functions_pc import irrigation_threshold_calculator, param_blocks, write_management_nsga_II, daisy_hame, read_sm, read_harvest, read_irrigation
from nsga_II_wwheat_eval_functions_pc import daisy_read_n_eval, daisy_kame
from nsga_II_wwheat_eval_functions_pc import param_blocks, daisy_kame
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.sampling.lhs import LatinHypercubeSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.termination import get_termination
from pymoo.optimize import minimize
import numpy as np
from pymoo.util.running_metric import RunningMetricAnimation
import pickle
import os
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymoo.indicators.hv import Hypervolume
import pandas as pd
import shutil
#from multiprocessing import starmap





class Wwheat(Problem):
    def __init__(self, main_dir, obs_yld_irr_path, obs_sm_30_60_path, cal_years, 
                 parameter_names, parameter_bounds, objfs, n_ieq, **kwargs):
        self.main_dir = main_dir
        self.obs_yld_irr_path = obs_yld_irr_path
        self.obs_sm_30_60_path = obs_sm_30_60_path
        self.cal_years = cal_years
        self.parameter_names = parameter_names
        self.lower_bounds = parameter_bounds[0,:]
        self.upper_bounds = parameter_bounds[1,:]
        self.len_var = parameter_bounds.shape[1]
        self.n_obj = len(objfs)
        self.n_ieq = n_ieq
        #self.gen_count = 0
        super().__init__(n_var = self.len_var, n_obj = self.n_obj, 
                         n_ieq_constr = self.n_ieq, 
                          xl = self.lower_bounds, 
                          xu = self.upper_bounds,
                          **kwargs)     #parameter limits hardcoded for each crop                          
        
    
    def _evaluate(self, X, out, *args, **kwargs):
        #print('Number of parameters sent for _evaluate',X.shape)
        #print(self.main_dir)
        #print(X)
        param_sm = param_blocks(X, 20)
        #print('Number of parameter blocks',len(param_sm))
        #print('Parameter block 0 sent from NSGAII population: ',param_sm[0])
        #print(param_sm)
        #print(f'Generation {gen_count} is being evaluated')
        #gen_count += 1
        F = daisy_kame(self.main_dir, self.parameter_names, param_sm, self.obs_yld_irr_path,
                       self.obs_sm_30_60_path , self.cal_years, True, True)
        out['F'] = F
        # out['G'] = G
    
#### analyzing convergence of GA
def analyse_convergence(all_evolution_save_path):
    with open(all_evolution_save_path, 'rb') as file:
        res = pickle.load(file)
    running = RunningMetricAnimation(delta_gen=5,
                        n_plots=8,
                        key_press=False,
                        do_show=True)

    for algorithm in res.history:
        running.update(algorithm)
        
def analyze_pareto(res_path,all_evolution_save_path, parameter_names):
    
    with open(all_evolution_save_path, 'rb') as file:
        result = pickle.load(file)
    obj = result.F
    var = result.X
    # three objectives
    x = obj[:,0]            #yield
    y = obj[:,1]            #irrigation
    z = obj[:,2]*100            #soil moisture, in mm/mm (converted to 100% for later)
    # t = obj[:,3]            #irrigation water use efficiency
    print('Variable shape', var.shape)
    df = pd.DataFrame(var, columns = parameter_names)
    df['RMSE_Yield'] = x
    df['RMSE_Irrigation'] = y
    df['RMSE_SM%'] = z
    # df['IWUE[tons/dm]'] = t 
    print(df)
    save_name = 'pareto_front_params_objfs.xlsx'
    save_path = os.path.join(res_path, save_name)
    df.to_excel(save_path, index = False)
    # print('qunatiles of nfk threshold used are: ',
    #       df['nfk_threshold'].quantile([0.1,0.25, 0.5, 0.75, 0.9]))
    
    
    
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot each axis as separate scatter 
    scatter = ax.scatter(x, y, z, marker = "o", color = "blue", label = "Objective functions at Pareto Front")
    #scatter = ax.scatter(x, y, z*100, marker = "o", color = "blue", label = "Objective functions at Pareto Front")
    #scatter1 = ax.scatter(x_gen1, y_gen1, z_gen1, marker = "o", color = "Orange", label = "First Gen")
    #scatter25 = ax.scatter(x_gen15, y_gen15, z_gen15, marker = "o", color = "red", label =  "15th Gen")
    
    # ax.set_xticks(x)
    # ax.set_yticks(y)
    ax.set_xlim(0, x.max())
    ax.set_ylim(0, y.max())
    # ax.set_zlim(0, t.max())
    ax.set_zlim(0, z.max())
    
    ax.set_title('Pareto Front')
    ax.set_xlabel('Yield RMSE (tDM/ha)')
    ax.set_ylabel('Irrigation RMSE (mm)')
    # ax.set_zlabel('IWUE [tons/dm]')
    ax.set_zlabel('Soil Moisture RMSE [%]')
    #ax.view_init(elev=25, azim=20)
    
    # ax.view_init(elev=25, azim=135)
    plt.legend()
    
    plt.tight_layout()
    file_name = 'pareto_front.png'
    plt.savefig(os.path.join(res_path,file_name))
    plt.show()

def analyze_hypervolume(res_path, all_evolution_save_path):
    with open(all_evolution_save_path, 'rb') as file:
        result = pickle.load(file)
    obj = result.F
    # no of evaluations arraz
    n_evals = np.array([e.evaluator.n_eval for e in result.history])
    
    # objs result in each gen
    hist_F = []
    for gen in result.history:
        f = gen.opt.get('F')
        hist_F.append(f)
    
    # worst and best result
    ideal_F = obj.min(axis = 0)
    nadir_F = obj.max(axis = 0)
    
    # hyper volume calculation
    metrics = Hypervolume(
        ref_point = np.array([1.0, 1.0, 1.0]),
        norm_ref_point = False,
        zero_to_one = True,
        ideal = ideal_F,
        nadir = nadir_F,
    
    )
    
    hvs = [metrics.do(_F) for _F in hist_F]
    
    # plot
    plt.figure(figsize=(12, 5))
    plt.plot(n_evals, hvs,  color='black', lw=0.7, label="Algorithm Performance")
    plt.scatter(n_evals, hvs,  facecolor="none", edgecolor='black', marker="p")
    plt.title("Convergence")
    plt.xlabel("Function Evaluations")
    plt.ylabel("Hypervolume")
    plt.tight_layout()
    file_name = 'Hypervolume.png'
    plt.savefig(os.path.join(res_path,file_name))
    plt.show()    


def plot_column_stats(df, group_size):
    n_groups = len(df) // group_size
    n_cols = df.shape[1]
    
    # Create a figure with subplots
    fig, axes = plt.subplots(n_cols, 1, figsize=(10, 3*n_cols))
    
    
    for i, col in enumerate(df.columns):
        # Reshape the column into groups of 195 rows
        grouped = df[col].values.reshape(n_groups, group_size)
        
        # Calculate statistics
        means = np.quantile(grouped, 0.5, axis=1)#grouped.mean(axis=1)
        mins = np.quantile(grouped, 0.1, axis=1)#grouped.min(axis=1)
        maxs = np.quantile(grouped, 0.9, axis=1)#grouped.max(axis=1)
        
        # Create x-axis values (group indices)
        x = np.arange(n_groups)
        
        # Plot mean line
        axes[i].plot(x, means, label='Median', color='blue')
        
        # Plot shaded area for min-max range
        axes[i].fill_between(x, mins, maxs, alpha=0.3, color='blue', label='Q10-Q90 range')
        
        # Add labels and title
        axes[i].set_title(f'{col}')
        axes[i].set_xlabel('Generation')
        axes[i].set_ylabel('Parametetr Value')
        axes[i].legend()
        axes[i].grid(True)
    
    plt.tight_layout()
    file_name = 'parameter_convergence.png'
    plt.savefig(os.path.join(res_path,file_name))
    plt.show()
def execute_pareto(main_dir, res_path, pareto_param_save_path, parameter_names, yld_obs, sm_obs, all_years, NSGA = False, calibration = False):
    pareto_parameters = np.load(pareto_param_save_path)
    print(pareto_parameters.shape)
    parameter_blocks = param_blocks(pareto_parameters, 20)
    # execute pareto front
    daisy_kame(main_dir, parameter_names, parameter_blocks, yld_obs, sm_obs, all_years, NSGA, calibration)
    pareto_simualtions = os.path.join(res_path, 'save_dir')
    if os.path.exists(pareto_simualtions):
        shutil.rmtree(pareto_simualtions)
    shutil.copytree(os.path.join(main_dir, 'save_dir'), pareto_simualtions)
    
def plot_pareto_simulation_outputs(main_dir, yld_path, sm_path, years, NSGA = False):
    #read the pareto outputs
    yield_par, irr_par, sm_par = daisy_read_n_eval(main_dir, yld_path, sm_path, years, NSGA = False)
    years = np.array(years)
    #read the observed values for all objectives
    df_obsyld = pd.read_excel(yld_path, header = [0])
    df_obsyld = df_obsyld [df_obsyld ['years'].isin(years)]
    obsyld = df_obsyld.iloc[:,1].to_numpy() 
    obsirr = df_obsyld.iloc[:,2].to_numpy()
    df_obssm = pd.read_excel(sm_path, header = [0])
    obssm = df_obssm['obs moist_30 cm'].to_numpy()                        #index hardcoded: index hardcoded
    sm_dates = df_obssm['Date'].to_numpy()                                 #index hardcoded: index hardcoded
    print(obsyld)
    
    for pareto_ind in range(yield_par.shape[0]):
            fig, ax = plt.subplots(figsize=(9, 9))
            ax_irr = ax.twinx()
            yield_ind = yield_par[pareto_ind,:]
            irr_ind = irr_par[pareto_ind,:]
            yield_min = np.max(yield_par, axis = 0)
            scatter1 = ax.plot(years,yield_ind, label='Simulated Yield', marker='o', linestyle='-', zorder=3)
            # scatter4 = ax.plot(years, yield_ind/irr_ind *100/2,label='Simulated IWUE', marker='x', linestyle='-.', zorder=3)
            # scatter5 = ax.plot(years, obsyld/obsirr *100/2,label='Observed IWUE', marker='x', linestyle='-.', zorder=3)
            scatter6 = ax.axvline(2013.5, linestyle = '-.', label = 'Validation period')
            bar1 = ax_irr.bar(years-0.125, irr_ind/10, width=0.25, alpha=0.5, label='Simulated Irrigation', zorder=2, color = 'orange')
            scatter2 = ax.plot(years,obsyld, label='Observed Yield', marker='o', linestyle='-',zorder=3)
            bar2 = ax_irr.bar(years+0.125, obsirr/10, width=0.25, alpha=0.5, label='Observed Irrigation', zorder=2, color = 'blue')
            ax.set_xticks(years)
            ax.set_ylabel('Yield [ton/ha]')
            ax.set_yticks(np.arange(0,21,1))
            ax_irr.set_ylabel('Irrigation [cm]')
            ax_irr.set_yticks(np.arange(0,29,2.5))
            plt.title(f'Pareto Parameter Set {pareto_ind}')
            yld_rmse1 = np.sqrt(np.mean(yield_ind-obsyld)**2)
            irr_rmse1 = np.sqrt(np.mean(irr_ind-obsirr)**2)
            txt = f'yield rmse = {yld_rmse1}\nirrigation rmse = {irr_rmse1}'
            #ax.text((2013, 10), txt, zorder=3)
            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax_irr.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
            #print(yield_ind/irr_ind *100)
            file = f'Pareto Parameter Set {pareto_ind}.png'
            file_path = os.path.join(res_path, file)
            
            # plt.text(1.1,0.5, txt)
            plt.savefig(file_path)
            plt.show()





def analyze_parameter_history_convergence(res_path, all_evolution_save_path, parameter_names):
    
    #1. collect all parameter values sequentially
    with open(all_evolution_save_path, 'rb') as file:
        result = pickle.load(file)
    all_parameter_dict = dict()
    for parameter in parameter_names:
        #print(parameter)
        all_parameter_dict[parameter] = dict()
        par_id = parameter_names.index(parameter)
        for generation in range(len(result.history)):
            #print(generation)
            for individual in range(len(result.history[generation].pop)):
                individual_count = generation * len(result.history[generation].pop) + individual
                #print(generation,individual, individual_count, end = ':')
                #print(parameter, individual_count)
                all_parameter_dict[parameter][individual_count] = result.history[generation].pop[individual].X[par_id]
    
    #print(all_parameter_dict)
    df_param_history = pd.DataFrame(all_parameter_dict)
    #print(df_param_history)
   
    #.plot parameter values individually
    plot_column_stats(df_param_history, len(result.history[0].pop))
    
    for parameter in df_param_history.columns:
        print('parameter range in pareto front: ', parameter, df_param_history[parameter].min(), df_param_history[parameter].max())
    

def analyze_pareto_parameters_spread(res_path, pareto_param_save_path, parameter_names):
    res_par = np.load(pareto_param_save_path)
    fig, axes = plt.subplots(res_par.shape[1], 1, figsize=(10, 3*res_par.shape[1]))
    for i in range(res_par.shape[1]):
        par_values = res_par[:,i]
        # median = np.quantile(par_values, 0.5)
        # q10 = np.quantile(par_values, 0.1)#grouped.min(axis=1)
        # q90 = np.quantile(par_values, 0.9)#grouped.max(axis=1)
        
        
        # x = np.arange(par_values.shape[0])
        axes[i].boxplot(par_values)
        # axes[i].fill_between(x, q10, q90, alpha=0.3, color='blue', label='Q10-Q90 range')
        
        
        axes[i].set_title(parameter_names[i])
    file_name = 'pareto_spread.png'
    plt.savefig(os.path.join(res_path,file_name))

if __name__ == "__main__":
    main_dir = r'C:\Users\Adhikari\Desktop\wwheat_stepwise\2_a_wwheat_nsga\random_try' #where nsga will be executed
    
    #observed yield and irrigation path
    observed_yield_irr_path = r'C:\Users\Adhikari\Desktop\Thesis\Thesis Shared Cloud\all_crop_field_results.xlsx'
    observed_sm_30_60_path = r'C:\Users\Adhikari\Desktop\Thesis\Thesis Shared Cloud\SoilMoisture_Observed_OP_30and60cm.xlsx'
    
    #years to calibrate for
    years_calibration = [_ for _ in range(2014,2023) if _ != 2012]
    #calibration + validation years
    years_all = [_ for _ in range(2008,2023) if _ != 2012]
    year_classes = ['normal', 'wet', 'dry']
    #parameters to change
    parameter_names_single = ['EmrTSUM','DSRate1', 'DSRate2', #devel
                        'Fm','Qeff',  #photosynthesis
                        'SPLAI', #canopy 
                        "PenPar1", "PenPar2", #root 
                        'fRoot_emr','fRoot_elong', 'fLeaf_emr', 'fLeaf_elong', 'fLeaf_flag', #partition
                        'E_leaf','r_Sorg', 'LfDR_anthesis', 'LfDR_fill', #Production
                        "nfkthres"] #irrigation
    parameters_lower_bounds = [ 80, 0.0208, 0.02, 4, 0.04, 0.0176, 0.1,   0, 0.45, 0.405,0.648,0.522,0.477,0.67,0.01, 0.0216,0.052, 40]
    parameters_upper_bounds = [120, 0.0312, 0.03, 6, 0.06, 0.0264, 0.225, 5, 0.55, 0.495,0.792,0.638,0.583,0.72,0.02, 0.0324,0.078, 60]
    
    #creates parameter names and parameter ranges for normal, wet, and dry year. 
    parameters_names = [i + '_' + j for j in year_classes for i in parameter_names_single] #this assumes that all parameters are time variant
    parameter_bounds = np.array([parameters_lower_bounds * len(year_classes), parameters_upper_bounds * len(year_classes)])
    
    
    #objective functions used in nsga, name only for ease access, these objfs are hardcoded in the read_n_eval function
    objfs = ['rmse_yield[tdm/ha]', 'rmse_irr [mm]', 'rmse_sm [%]']
    #number of constraints, none
    n_ieq = 0
    
    
    #problem definition
    problem = Wwheat(main_dir, observed_yield_irr_path, observed_sm_30_60_path, years_calibration,
                     parameters_names, parameter_bounds, objfs, n_ieq)
    
    #nsga_II hyperparameters
    algorithm = NSGA2(pop_size= 200, n_offsprings= 200,
                          sampling=LatinHypercubeSampling(),
                          crossover=SBX(prob=0.5, eta=15),
                          mutation=PM(prob=0.2, eta=5),
                          eliminate_duplicates=True)
    #termination criteria
    termination = get_termination("n_gen", 40)
    
    #execute the nsga for 5 times at different random initializations
    for i in range(5):
        print(i)
        #initialize nsga at random location in parameter space
        random_seed = int(time.time())
        parameter_names = problem.parameter_names
        #run nsga
        res = minimize(problem, algorithm, termination,
                        seed=random_seed, save_history=True,
                        verbose=True)
        #save the results                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        res_path = os.path.join(main_dir,str(i))
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        #save pareto parameter set
        pareto_param_save_path = os.path.join(res_path, 'pareto_parameters.npy')
        np.save(pareto_param_save_path, res.X)
        #save pareto objective functions
        pareto_obfs_save_path = os.path.join(res_path, 'pareto_objfs.npy')
        np.save(pareto_obfs_save_path, res.F)
        #save all population from all generations
        all_evolution_save_path = os.path.join(res_path, 'all_gen_all_population.pkl')
        with open(all_evolution_save_path, 'wb') as f:
            
            pickle.dump(res, f)
        
        
        
        #analyze the nsga run 
        analyse_convergence(all_evolution_save_path)
        
        analyze_pareto(res_path, all_evolution_save_path, parameter_names)
        
        analyze_hypervolume(res_path,all_evolution_save_path)
        
        analyze_parameter_history_convergence(res_path, all_evolution_save_path, parameter_names)
        
        analyze_pareto_parameters_spread(res_path, pareto_param_save_path, parameter_names)
        execute_pareto(main_dir, res_path, pareto_param_save_path, parameter_names,
                       observed_yield_irr_path, observed_sm_30_60_path, years_all, NSGA = False)
        plot_pareto_simulation_outputs(res_path, observed_yield_irr_path, observed_sm_30_60_path, years_all, NSGA = False)
        # break