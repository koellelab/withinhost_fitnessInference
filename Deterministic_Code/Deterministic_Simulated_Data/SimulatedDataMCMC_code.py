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
            while (l<l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
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
            while (l<l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
                #calculate mutant fitness under this (k+1,l) scenario
                sum_check_w = sum_check_w + Pkl(mutant_MOI,wild_MOI,k,l)
                wild_fitness_sum = wild_fitness_sum + l * Pkl(mutant_MOI,wild_MOI,k,l) * Fkl(sigma_m,k,l)
                l = l + 1
            k = k + 1
        wild_fitness_sum = wild_fitness_sum + l * Fkl(sigma_m,(k_max+1),(l_max+1))** 0.001 * 0.001
            
            
        mutant_fitness = mutant_fitness_sum / mutant_MOI
        wild_fitness = wild_fitness_sum / wild_MOI
    
        #Determine mutant freq into next generation
        freq = freq * mutant_fitness / (freq * mutant_fitness + (1-freq)* wild_fitness)  
        record = np.concatenate((record,np.array([freq]))) 
    #LOOP ENDS HERE

    return record


# In[ ]:


trend = Main_Deterministic(2,np.log(1.5),0.1,25)
V = 100
sample=np.random.beta(trend*V,(1-trend)*V)
time_point = np.array([5,10,15,20,25])
data = sample[time_point]
data


# In[ ]:


def LogFitFour50(Init_freq,log_MOI,log_fitness):
    v = 100 #True value here
    simulated = Main_Deterministic(np.exp(log_MOI),log_fitness,Init_freq,25)
    prob = 1
    for i in np.arange(len(time_point)):
        prob = prob * (beta.pdf(sample[time_point[i]],simulated[time_point[i]]*v,(1-simulated[time_point[i]])*v))
    return prob   


# In[ ]:


def logMOIprior(logMOI):
    return norm.pdf(logMOI,np.log(2),0.5)


# In[ ]:


def MCMC(freq0,log_moi0,sigmam0,iterations):     #Initialize theta
    noise = 100 #true value here
    cur_likelihood = LogFitFour50(freq0,log_moi0,sigmam0)
    theta_table = np.zeros((iterations+1,4))
    theta_table[0,:] = np.array((freq0,log_moi0,sigmam0,cur_likelihood))
    
    for i in np.arange(iterations):
        #print(i)
        cur_likelihood = theta_table[i,3]
        #propose new theta set
        freq_new = np.random.normal(theta_table[i,0], np.sqrt(0.0001))
        logmoi_new = np.random.normal(theta_table[i,1], np.sqrt(0.02))
        sigmam_new = np.random.normal(theta_table[i,2], np.sqrt(0.002))
        new_likelihood = LogFitFour50(freq_new,logmoi_new,sigmam_new)
        A = logMOIprior(logmoi_new) * new_likelihood / (logMOIprior(theta_table[i,1])*cur_likelihood)
        #Generate a random value between 0 and 1
        a = random.random()
        if (A >a):
            theta_table[i+1,:] = np.array((freq_new,logmoi_new,sigmam_new,new_likelihood))
        else:
            theta_table[i+1,:] = theta_table[i,:]
    return theta_table


# In[ ]:


mcmc_table = MCMC(0.08,np.log(1.5),np.log(2),20000)


# In[ ]:


df = pd.DataFrame(mcmc_table, columns = ['initfreqe', 'logmoi', 'logfitness','likelihood']) 
df.to_csv('MCMC20000.csv')

