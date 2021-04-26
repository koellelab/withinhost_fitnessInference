#!/usr/bin/env python
# coding: utf-8

# In[1]:


import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
import math
from scipy.stats import beta
np.random.seed(30)
random.seed(30)


# In[2]:


df = pd.read_csv('MCMC20000.csv')
df


# In[3]:


init_freq = df[['initfreqe']].values
logmoi = df[['logmoi']].values
logfitness = df[['logfitness']].values

convertmoi = np.exp(logmoi)
convertfitness = np.exp(logfitness)


# In[4]:


fig, ax = plt.subplots(figsize=(12,3), dpi= 360, facecolor='w', edgecolor='k')
plt.subplot(1,3,1)
plt.plot(np.arange(20003),init_freq)
plt.text(0.5,0.18,"A",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.axhline(0.1, color="red")
plt.xlabel("MCMC iteration")
plt.ylabel("Initial frequency")

plt.subplot(1,3,2)
plt.plot(np.arange(20003),convertmoi)
plt.text(0.5,7.7,"B",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.axhline(2, color="red")
plt.xlabel("MCMC iteration")
plt.ylabel("MOI")

plt.subplot(1,3,3)
plt.plot(np.arange(20003),convertfitness)
plt.text(0.5,4.25,"C",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.axhline(1.5, color="red")
plt.xlabel("MCMC iteration")
plt.ylabel("Variant fitness")

plt.show()
fig.savefig("MCMC_trace.jpg",bbox_inches = 'tight')


# In[5]:


n = 20003


# In[6]:


interval = np.arange(2000,20050,50)
sample = np.zeros(len(interval)-1)
for i in np.arange(len(interval)-1):
    sample[i] = random.randint(interval[i],interval[i+1])
sample = sample.astype(int)


# In[7]:


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
        while (k<k_max):
            l = 0
            #while(l<100):
            while (l<l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
                #calculate mutant fitness under this (k+1,l) scenario
                sum_check_m = sum_check_m + Pkl(mutant_MOI,wild_MOI,k,l)
                mutant_fitness_sum = mutant_fitness_sum + k * Pkl(mutant_MOI,wild_MOI,k,l) * Fkl(sigma_m,k,l)
                l = l + 1
            k = k + 1
        mutant_fitness_sum = mutant_fitness_sum + k * Fkl(sigma_m,k_max,l_max)* 0.000001 * 0.000001
        
        #Calculate wild type fitness of each generation enumerating all possible coinfection scenarios
        sum_check_w = 0
        wild_fitness_sum = 0 
        k = 0
        while (k<k_max):       #some check on poisson.pmf(k+1,MOI) to avoid infinite sums
            l = 0
            while (l<l_max):   #some check on poisson.pmf(k+l+1,MOI) to avoid infinite sums
                #calculate mutant fitness under this (k+1,l) scenario
                sum_check_w = sum_check_w + Pkl(mutant_MOI,wild_MOI,k,l)
                wild_fitness_sum = wild_fitness_sum + l * Pkl(mutant_MOI,wild_MOI,k,l) * Fkl(sigma_m,k,l)
                l = l + 1
            k = k + 1
        wild_fitness_sum = wild_fitness_sum + l * Fkl(sigma_m,k_max,l_max)* 0.001 * 0.001
            
            
        mutant_fitness = mutant_fitness_sum / mutant_MOI
        wild_fitness = wild_fitness_sum / wild_MOI
    
        #Determine mutant freq into next generation
        freq = freq * mutant_fitness / (freq * mutant_fitness + (1-freq)* wild_fitness)  
        record = np.concatenate((record,np.array([freq]))) 
    #LOOP ENDS HERE

    return record


# In[8]:


main_trend = Main_Deterministic(2,np.log(1.5),0.1,25)
V = 100
sample2=np.random.beta(main_trend*V,(1-main_trend)*V)
time_point = np.array([5,10,15,20,25])
data = sample2[time_point]
data


# In[9]:


index = np.zeros(10)
for i in np.arange(10):
    index[i] = int(random.randint(2000, 20000))
index


# In[10]:


trend = np.zeros((10,26))
for i in np.arange(10):
    current_init_freq= init_freq[int(index[i])]
    current_logmoi = logmoi[int(index[i])]
    current_logfitness = logfitness[int(index[i])]
    current = Main_Deterministic(np.exp(current_logmoi),current_logfitness,current_init_freq,25)
    for j in np.arange(26):
        trend[i,j] = current[j]


# In[11]:


df = pd.read_csv('LL100.csv')
df.drop_duplicates(['moi','fitness'], inplace=True)
log = df[['log']].values
fitness = df[['fitness']].values
moi = df[['moi']].values


# In[12]:


partition_MOI = 0.1
partition_fitness = 0.01

MOI_range = np.hstack(np.arange(0.1,5.1,partition_MOI)).round(decimals=2)
Fitness_range = np.arange(1.01,3.01,partition_fitness).round(decimals=2)

maxlog = max(log)
max_index=np.argmax(log)
minlog = min(log)
low_bound = maxlog-1.92


# In[13]:


fig, ax = plt.subplots(figsize=(8, 3), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(1,2,1)
for i in np.arange(10):
    plt.plot(np.arange(26),trend[i,:],alpha=0.3,label="Simulated",c="black")
plt.plot(time_point,data,'o',label="Simulated + Noise, V=50",c="firebrick")
plt.plot(np.arange(26),main_trend,label="Simulated + Noise, V=50",c="firebrick")
plt.ylim(0,1)
plt.xlabel("Generation")
plt.ylabel("Variant Frequency")
plt.text(0.05,0.9,"A",fontsize=15)

ax = plt.subplot(1,2,2)
extent = np.min(Fitness_range),np.max(Fitness_range),np.min(MOI_range),np.max(MOI_range)
pivot3 = df.pivot(index='moi', columns='fitness', values='log')
log = ax.imshow(pivot3.iloc[::-1], cmap="Blues", extent=extent)  #interpolation='nearest',
ct = plt.contour(Fitness_range,MOI_range,pivot3,levels=[low_bound],colors = 'white',extent = extent)
#plt.clabel(ct, colors = 'white', fmt = '%2.1f', fontsize=12)
plt.plot(1.5,2,"o",label="True",c="Red",markersize=10)
plt.plot(fitness[max_index],moi[max_index],"o",c="Orange",label="Max",markersize=10)
plt.gca().invert_yaxis()
ax.set_aspect('auto')
fig.colorbar(log)
plt.legend(loc="best")
plt.xlabel("Variant fitness")
plt.ylabel("MOI")
plt.text(1.1,0.55,"B",fontsize=15,color="white")
plt.show()
fig.savefig("fig2a.jpg",bbox_inches = 'tight')


# In[14]:


fig, ax = plt.subplots(figsize=(12, 3), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(1,3,1)
plt.hist(init_freq[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(init_freq[sample],0.5), color="black",alpha=0.8 )
plt.axvline(np.quantile(init_freq[sample],0.025), color="black", alpha=0.8, linestyle='--')
plt.axvline(np.quantile(init_freq[sample],0.975), color="black", alpha=0.8,linestyle='--')
plt.axvline(0.1, color="red")
plt.xlabel("Initial frequency")
plt.ylabel("Proportion")
plt.text(0.033,0.17,"C",fontsize=15)

plt.subplot(1,3,2)
plt.hist(convertmoi[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertmoi[sample],0.5), color="black",alpha=0.8)
plt.axvline(np.quantile(convertmoi[sample],0.025), color="black", alpha=0.8,linestyle='--')
plt.axvline(np.quantile(convertmoi[sample],0.975), color="black",alpha=0.8, linestyle='--')
plt.xlabel("MOI")
plt.axvline(2, color="red")
plt.text(0.7,0.215,"D",fontsize=15)

plt.subplot(1,3,3)
plt.hist(convertfitness[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertfitness[sample],0.5), color="black",alpha=0.8)
plt.axvline(np.quantile(convertfitness[sample],0.025), color="black", alpha=0.8,linestyle='--')
plt.axvline(np.quantile(convertfitness[sample],0.975), color="black", alpha=0.8,linestyle='--')
plt.xlabel("Variant fitness")
plt.axvline(1.5, color="red")
plt.text(3.9,0.32,"E",fontsize=15)

plt.show()
fig.savefig("fig2b.jpg",bbox_inches = 'tight')

