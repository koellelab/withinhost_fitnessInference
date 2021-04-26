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
from matplotlib.ticker import PercentFormatter

np.random.seed(30)
random.seed(30)


# In[2]:


wilker_site788_data = np.zeros((4,4))

#REPLACE WITH EXACT VALUES WHEN WE HEAT BACK FROM TOM
wilker_site788_data[0,:] = np.array([0.044,0.057,0.32,0.72])
wilker_site788_data[1,:] = np.array([0.044,0.188,0.449,0.683])
wilker_site788_data[2,:] = np.array([0.044,0.353,0.632,0.809])
wilker_site788_data[3,:] = np.array([0.044,0.202,0.721,0.59])   #Initial Frequency: 0.0447 , 0.94212871])
time_point = np.array([0,1,3,5])


# In[3]:


df = pd.read_csv('wilkerMCMC.csv')


# In[4]:


init_freq1 = df[['initfreq1']].values
init_freq2 = df[['initfreq2']].values
init_freq3 = df[['initfreq3']].values
init_freq4 = df[['initfreq4']].values
logmoi = df[['logmoi']].values
logfitness = df[['logfitness']].values

convertmoi = np.exp(logmoi)
convertfitness = np.exp(logfitness)


# In[5]:


fig, ax = plt.subplots(figsize=(16, 9), dpi= 360, facecolor='w', edgecolor='k')
plt.subplot(2,3,1)
plt.plot(np.arange(20002),init_freq1)
plt.text(0.5,0.072,"A",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("Initial frequency for ferret 13")

plt.subplot(2,3,2)
plt.plot(np.arange(20002),init_freq2)
plt.text(0.5,0.095,"B",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("Initial frequency for ferret 15")


plt.subplot(2,3,3)
plt.plot(np.arange(20002),init_freq3)
plt.text(0.5,0.162,"C",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("Initial frequency for ferret 17")


plt.subplot(2,3,4)
plt.plot(np.arange(20002),init_freq4)
plt.text(0.5,0.13,"D",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("Initial frequency for ferret 21")

plt.subplot(2,3,5)
plt.plot(np.arange(20002),convertmoi)
plt.text(0.5,8.3,"E",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("MOI")

plt.subplot(2,3,6)
plt.plot(np.arange(20002),convertfitness)
plt.text(0.5,9.3,"F",fontsize=15)
plt.axvline(2000, color="black", linestyle='--')
plt.xlabel("MCMC iteration")
plt.ylabel("Variant fitness")

plt.show()
fig.savefig("wilkerMCMC_trace.jpg",bbox_inches = 'tight')


# In[6]:


interval = np.arange(2000,20050,50)
sample = np.zeros(len(interval)-1)
for i in np.arange(len(interval)-1):
    sample[i] = random.randint(interval[i],interval[i+1])
sample = sample.astype(int)


# In[13]:


fig, ax = plt.subplots(figsize=(10, 7), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(2,2,1)
plt.hist(init_freq1[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(init_freq1[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(init_freq1[sample],0.025),color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(init_freq1[sample],0.975),color="black",alpha = 0.8, linestyle='--')
plt.text(0.003,0.215,"A",fontsize=15)
plt.xlabel("Initial frequency for ferret 13")
plt.ylabel("Proportion")

plt.subplot(2,2,2)
plt.hist(init_freq2[sample],weights=np.ones(len(sample)) / len(sample))
plt.text(0.017,0.23,"B",fontsize=15)
plt.axvline(np.quantile(init_freq2[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(init_freq2[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(init_freq2[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("Initial frequency for ferret 15")


plt.subplot(2,2,3)
plt.hist(init_freq3[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(init_freq3[sample],0.5), color="black",alpha = 0.8)
plt.text(0.06,0.24,"C",fontsize=15)
plt.axvline(np.quantile(init_freq3[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(init_freq3[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("Initial frequency for ferret 17")
plt.ylabel("Proportion")

plt.subplot(2,2,4)
plt.hist(init_freq4[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(init_freq4[sample],0.5), color="black",alpha = 0.8)
plt.text(0.05,0.215,"D",fontsize=15)
plt.axvline(np.quantile(init_freq4[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(init_freq4[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("Initial frequency for ferret 21")
plt.show()
fig.savefig("wilker_freq_post.jpg",bbox_inches = 'tight')


# In[14]:


fig, ax = plt.subplots(figsize=(20, 3), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(1,4,1)
plt.plot(np.array([0,1,3]),wilker_site788_data[0,:3],'o-',label="ferret13",color="blue")
plt.plot(np.array([3,5]),wilker_site788_data[0,2:],'o--',color="blue")
plt.plot(np.array([0,1,3]),wilker_site788_data[1,:3],'o-',label="ferret15",color="orange")
plt.plot(np.array([3,5]),wilker_site788_data[1,2:],'o--',color="orange")
plt.plot(np.array([0,1,3]),wilker_site788_data[2,:3],'o-',label="ferret17",color="red")
plt.plot(np.array([3,5]),wilker_site788_data[2,2:],'o--',color="red")
plt.plot(np.array([0,1,3]),wilker_site788_data[3,:3],'o-',label="ferret21",color="green")
plt.plot(np.array([3,5]),wilker_site788_data[3,2:],'o--',color="green")
plt.legend(loc=4)
plt.xlabel("Day")
plt.ylabel("G788A variant frequency")
plt.text(0.07,0.73,"A",fontsize=15)

plt.subplot(1,4,2)
plt.hist(convertmoi[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertmoi[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(convertmoi[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(convertmoi[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("MOI")
plt.text(1.8,0.19,"B",fontsize=15)
plt.ylabel("Proportion")


plt.subplot(1,4,3)
plt.hist(convertfitness[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertfitness[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(convertfitness[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(convertfitness[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("Variant fitness")
plt.text(9,0.29,"C",fontsize=15)
plt.ylabel("Proportion")

plt.subplot(1,4,4)
plt.plot(convertfitness,convertmoi,".")   #MOI; Fitness
plt.xlabel("Variant fitness")
plt.ylabel("MOI")
plt.text(9.2,2.2,"D",fontsize=15)
plt.gca().invert_yaxis()
plt.show()

fig.savefig("fig3.jpg",bbox_inches = 'tight')


# In[15]:


print(np.quantile(convertfitness[sample],0.5))
print(np.quantile(convertfitness[sample],0.025))
print(np.quantile(convertfitness[sample],0.975))


# In[16]:


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


# In[17]:


def NoCoinfection(a_fitness,init_freq):
    initial_a_frequency = init_freq
    frequency_count = np.array([1-initial_a_frequency,initial_a_frequency])
    other_fitness =1
    t = 25 #next twenty generations
    record = np.zeros(26)
    for x in np.arange(t)+1:
        record[x] = frequency_count[1]
        count_a = frequency_count[1] * a_fitness
        count_others = frequency_count[0] * other_fitness
        norm_factor = count_a + count_others
        frequency_count[1] = count_a / norm_factor
        frequency_count[0] = 1 - frequency_count[1]
    return record


# In[18]:


index = np.zeros(10)
for i in np.arange(10):
    index[i] = int(random.randint(200, 20000))
index


# In[19]:


trend1 = np.zeros((10,26))
noncoinfect1 = np.zeros((10,26))
trend2 = np.zeros((10,26))
noncoinfect2 = np.zeros((10,26))
trend3 = np.zeros((10,26))
noncoinfect3 = np.zeros((10,26))
trend4 = np.zeros((10,26))
noncoinfect4 = np.zeros((10,26))

for i in np.arange(10):
    current_init_freq1= init_freq1[int(index[i])]
    current_init_freq2= init_freq2[int(index[i])]
    current_init_freq3= init_freq3[int(index[i])]
    current_init_freq4= init_freq4[int(index[i])]
    current_logmoi = logmoi[int(index[i])]
    current_logfitness = logfitness[int(index[i])]
    current = Main_Deterministic(np.exp(current_logmoi),current_logfitness,current_init_freq1,25)
    current2 = NoCoinfection(np.exp(current_logfitness),current_init_freq1)
    for j in np.arange(26):
        trend1[i,j] = current[j]
        noncoinfect1[i,j] = current2[j]
    current = Main_Deterministic(np.exp(current_logmoi),current_logfitness,current_init_freq2,25)
    current2 = NoCoinfection(np.exp(current_logfitness),current_init_freq2)
    for j in np.arange(26):
        trend2[i,j] = current[j]
        noncoinfect2[i,j] = current2[j]
    current = Main_Deterministic(np.exp(current_logmoi),current_logfitness,current_init_freq3,25)
    current2 = NoCoinfection(np.exp(current_logfitness),current_init_freq3)
    for j in np.arange(26):
        trend3[i,j] = current[j]
        noncoinfect3[i,j] = current2[j]
    current = Main_Deterministic(np.exp(current_logmoi),current_logfitness,current_init_freq4,25)
    current2 = NoCoinfection(np.exp(current_logfitness),current_init_freq4)
    for j in np.arange(26):
        trend4[i,j] = current[j]
        noncoinfect4[i,j] = current2[j]


# In[21]:


fig, ax = plt.subplots(figsize=(18, 3), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(1,4,1)
for i in np.arange(10):
    plt.plot(np.arange(16)/3,trend1[i,:16],alpha=0.5,label="Simulated",c="black")
    plt.plot(np.arange(16)/3,noncoinfect1[i,1:17],alpha=0.5,label="Simulated",c="purple")
plt.plot(time_point,wilker_site788_data[0,:],'o',label="Simulated + Noise, V=50",c="firebrick")
plt.ylim(0,1)
plt.xlabel("Day")
plt.ylabel("Variant frequency for ferret 13")
plt.text(0.17,0.85,"A",fontsize=15)


plt.subplot(1,4,2)
for i in np.arange(10):
    plt.plot(np.arange(16)/3,trend2[i,:16],alpha=0.5,label="Simulated",c="black")
    plt.plot(np.arange(16)/3,noncoinfect2[i,1:17],alpha=0.5,label="Simulated",c="purple")
plt.plot(time_point,wilker_site788_data[1,:],'o',label="Simulated + Noise, V=50",c="firebrick")
plt.ylim(0,1)
plt.xlabel("Day")
plt.ylabel("Variant frequency for ferret 15")
plt.text(0.17,0.85,"B",fontsize=15)

plt.subplot(1,4,3)
for i in np.arange(10):
    plt.plot(np.arange(16)/3,trend3[i,:16],alpha=0.5,label="Simulated",c="black")
    plt.plot(np.arange(16)/3,noncoinfect3[i,1:17],alpha=0.5,label="Simulated",c="purple")
plt.plot(time_point,wilker_site788_data[2,:],'o',label="Simulated + Noise, V=50",c="firebrick")
plt.ylim(0,1)
plt.xlabel("Day")
plt.ylabel("Variant frequency for ferret 17")
plt.text(0.17,0.85,"C",fontsize=15)

plt.subplot(1,4,4)
for i in np.arange(10):
    plt.plot(np.arange(16)/3,trend4[i,:16],alpha=0.5,label="Simulated",c="black")
    plt.plot(np.arange(16)/3,noncoinfect4[i,1:17],alpha=0.5,label="Simulated",c="purple")

plt.plot(time_point,wilker_site788_data[3,:],'o',label="Simulated + Noise, V=50",c="firebrick")
plt.ylim(0,1)
plt.xlabel("Day")
plt.ylabel("Variant frequency for ferret 21")
plt.text(0.17,0.85,"D",fontsize=15)
plt.show()
fig.savefig("fig4.jpg",bbox_inches = 'tight')


# In[ ]:




