#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
import math
from scipy.optimize import curve_fit
import scipy
from scipy.stats import t
import sys
np.random.seed(30)
random.seed(30)
from scipy.stats import beta
from scipy.stats import norm
from scipy.stats import uniform
import seaborn as sns
import pandas as pd


# In[ ]:


## Deterministic Take 3
def Fkl(sigma_m,k,l):
    if (k+l)==0:
        cell_fitness = 0
    else: 
        cell_fitness = (k)/(k+l)*np.exp(sigma_m) + l/(k+l)
    return cell_fitness

def Pkl(mutant_MOI,wild_MOI,k,l):
    return poisson.pmf(k,mutant_MOI) * poisson.pmf(l,wild_MOI)

#Main loop here
def Main_Deterministic(MOI,sigma_m,initial_freq,generations):
    freq = initial_freq
    record = np.array([freq])
    
    #LOOP STARTS FROM HERE
    for i in np.arange(generations)+1:
        #Updating each generation b/c mutant freq is changing
        mutant_MOI = MOI * freq
        wild_MOI = MOI * (1-freq)

        k_max = 0
        while (poisson.cdf(k_max,mutant_MOI)<0.9999999):
            k_max = k_max + 1
        l_max = 0
        while (poisson.cdf(l_max,wild_MOI)<0.9999999):
            l_max = l_max + 1
        
        #Calculate mutant fitness of each generation enumerating all possible coinfection scenarios
        mutant_fitness_sum = 0 
        sum_check_m = 0
        k = 0
        #while(k<100):
        #while ((k+1)<5 or (k+1)<MOI or poisson.pmf(k+1,MOI)>1e-3):        #some check on poisson.pmf(k+1,MOI) to avoid infinite sums
        while (k<=k_max):
            l = 0
            #while(l<100):
            while (l<=l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
                #calculate mutant fitness under this (k+1,l) scenario
                sum_check_m = sum_check_m + Pkl(mutant_MOI,wild_MOI,k,l)
                mutant_fitness_sum = mutant_fitness_sum + k * Pkl(mutant_MOI,wild_MOI,k,l) * Fkl(sigma_m,k,l)
                l = l + 1
            k = k + 1
        mutant_fitness_sum = mutant_fitness_sum + k * Fkl(sigma_m,(k_max+1),(l_max+1))* 0.000001 * 0.000001
        
        #Calculate wild type fitness of each generation enumerating all possible coinfection scenarios
        sum_check_w = 0
        wild_fitness_sum = 0 
        k = 0
        while (k<=k_max):       #some check on poisson.pmf(k+1,MOI) to avoid infinite sums
            l = 0
            while (l<=l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
                #calculate mutant fitness under this (k+1,l) scenario
                sum_check_w = sum_check_w + Pkl(mutant_MOI,wild_MOI,k,l)
                wild_fitness_sum = wild_fitness_sum + l * Pkl(mutant_MOI,wild_MOI,k,l) * Fkl(sigma_m,k,l)
                l = l + 1
            k = k + 1
        wild_fitness_sum = wild_fitness_sum + l * Fkl(sigma_m,(k_max+1),(l_max+1))* 0.000001 * 0.000001
            
            
        mutant_fitness = mutant_fitness_sum / mutant_MOI
        wild_fitness = wild_fitness_sum / wild_MOI
    
        #Determine mutant freq into next generation
        freq = freq * mutant_fitness / (freq * mutant_fitness + (1-freq)* wild_fitness)  
        record = np.concatenate((record,np.array([freq]))) 
    #LOOP ENDS HERE

    return record


# In[ ]:


wilker_site788_data = np.zeros((4,2))

#REPLACE WITH EXACT VALUES WHEN WE HEAT BACK FROM TOM
wilker_site788_data[0,:] = np.array([0.057,0.32])
wilker_site788_data[1,:] = np.array([0.188,0.449])
wilker_site788_data[2,:] = np.array([0.353,0.632])
wilker_site788_data[3,:] = np.array([0.202,0.721])   #Initial Frequency: 0.044


# In[ ]:


#Without using Day 0 data
def LogFitFour_Wilker(Init_freq1,Init_freq2,Init_freq3,Init_freq4,log_MOI,log_Fitness):
    v = 100 #Change noise level here
    #v = 25
    #v = 400
    
    generation_time = 8  #Change generation time here
    #generation_time = 6
    #generation_time = 12
    
    generation_multiplier = math.floor(24 / generation_time)
    simulated1 = Main_Deterministic(np.exp(log_MOI),log_Fitness,Init_freq1,3*generation_multiplier+1)  
    pdf0 = beta.pdf(0.044,simulated1[0]*v,(1-simulated1[0])*v)    #Change mock data accordingly
    pdf1 = beta.pdf(wilker_site788_data[0,0],simulated1[generation_multiplier]*v,(1-simulated1[generation_multiplier])*v)    #Change mock data accordingly
    pdf3 = beta.pdf(wilker_site788_data[0,1],simulated1[3*generation_multiplier]*v,(1-simulated1[3*generation_multiplier])*v)
    
    simulated2 = Main_Deterministic(np.exp(log_MOI),log_Fitness,Init_freq2,3*generation_multiplier+1)  
    pdf20 = beta.pdf(0.044,simulated2[0]*v,(1-simulated2[0])*v)
    pdf21 = beta.pdf(wilker_site788_data[1,0],simulated2[generation_multiplier]*v,(1-simulated2[generation_multiplier])*v)    #Change mock data accordingly
    pdf23 = beta.pdf(wilker_site788_data[1,1],simulated2[3*generation_multiplier]*v,(1-simulated2[3*generation_multiplier])*v)
   
    simulated3 = Main_Deterministic(np.exp(log_MOI),log_Fitness,Init_freq3,3*generation_multiplier+1)  
    pdf30 = beta.pdf(0.044,simulated3[0]*v,(1-simulated3[0])*v)
    pdf31 = beta.pdf(wilker_site788_data[2,0],simulated3[generation_multiplier]*v,(1-simulated3[generation_multiplier])*v)    #Change mock data accordingly
    pdf33 = beta.pdf(wilker_site788_data[2,1],simulated3[3*generation_multiplier]*v,(1-simulated3[3*generation_multiplier])*v)
    
    simulated4 = Main_Deterministic(np.exp(log_MOI),log_Fitness,Init_freq4,3*generation_multiplier+1)  
    pdf40 = beta.pdf(0.044,simulated4[0]*v,(1-simulated4[0])*v)
    pdf41 = beta.pdf(wilker_site788_data[3,0],simulated4[generation_multiplier]*v,(1-simulated4[generation_multiplier])*v)    #Change mock data accordingly
    pdf43 = beta.pdf(wilker_site788_data[3,1],simulated4[3*generation_multiplier]*v,(1-simulated4[3*generation_multiplier])*v)

    return pdf0*pdf1*pdf3*pdf20*pdf21*pdf23*pdf30*pdf31*pdf33*pdf40*pdf41*pdf43


# In[ ]:


def logMOIprior(logMOI):
    return norm.pdf(logMOI,np.log(4),0.4)


# In[ ]:


def MCMC(freq1, freq2, freq3, freq4,log_moi0,sigmam0,iterations):     #Initialize theta
    noise = 100 ##Change noise level here
    #noise =25
    #noise =400
    cur_likelihood = LogFitFour_Wilker(freq1,freq2,freq3,freq4,log_moi0,sigmam0)
    #theta_table = np.zeros((iterations+1,7))
    theta_table[0,:] = np.array((freq1,freq2,freq3,freq4,log_moi0,sigmam0,cur_likelihood))
    
    for i in np.arange(iterations):
        print(i)
        cur_likelihood = theta_table[i,6]
        #propose new theta set
        freq1_new = np.random.normal(theta_table[i,0], 0.01)
        freq2_new = np.random.normal(theta_table[i,1], 0.01)
        freq3_new = np.random.normal(theta_table[i,2], 0.01)
        freq4_new = np.random.normal(theta_table[i,3], 0.01)
        logmoi_new = np.random.normal(theta_table[i,4], np.sqrt(0.02))
        sigmam_new = np.random.normal(theta_table[i,5], np.sqrt(0.002))
        new_likelihood = LogFitFour_Wilker(freq1_new,freq2_new,freq3_new,freq4_new,logmoi_new,sigmam_new)
        A = logMOIprior(logmoi_new) * new_likelihood / (logMOIprior(theta_table[i,4])*cur_likelihood)
        #Generate a random value between 0 and 1
        a = random.random()
        if (A >a):
            theta_table[i+1,:] = np.array((freq1_new,freq2_new,freq3_new,freq4_new,logmoi_new,sigmam_new,new_likelihood))
        else:
            theta_table[i+1,:] = theta_table[i,:]
    return theta_table


# In[ ]:


init1 = 0.08
init2 = 0.08
init3 = 0.08
init4 = 0.08
logmoi = np.log(1.5)
logfitness = np.log(2)


# In[ ]:


iterations= 20000
theta_table = np.zeros((iterations+1,7))   #theta_table: program can be interrupted
mcmc_table = MCMC(init1,init2,init3,init4,logmoi,logfitness,iterations)


# In[ ]:


df = pd.DataFrame(theta_table, columns = ['initfreq1','initfreq2','initfreq3','initfreq4','logmoi', 'logfitness','likelihood']) 
df.to_csv('wilkerMCMC.csv')

