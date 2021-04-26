#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
import math
random.seed(30)


# In[2]:


def stochastic(N,C,initial_freq,fitness,generation):
    record = np.zeros(generation + 1)
    record[0] = initial_freq
    i = 1
    Nm = N * initial_freq
    Nw = N - Nm
    while (i<= generation and Nm>0 and Nw>0):
        mutant_infect = np.random.multinomial(Nm, [1/C]*C, size=1)
        #print("mutant infection: ",mutant_infect)
        wild_infect = np.random.multinomial(Nw, [1/C]*C, size=1)
        #print("wild type infection: ",wild_infect)
        total_number = mutant_infect + wild_infect
        #print("total infection numbers: ", total_number)
        total_fitness = mutant_infect * fitness + wild_infect
        average_fitness = total_fitness / total_number
        average_fitness = np.nan_to_num(average_fitness, nan=0.0)
        #print(average_fitness)
        mutant_average_fitness = np.sum((average_fitness * mutant_infect))/Nm
        wild_average_fitness = np.sum((average_fitness * wild_infect))/Nw
        #print(mutant_average_fitness)
        #print(wild_average_fitness)
        p_binom = record[i-1] * mutant_average_fitness / (record[i-1] * mutant_average_fitness + (1 - record[i-1]) * wild_average_fitness)
        record[i] = np.random.binomial(N,p_binom)/N
        Nm = round( N * record[i] )
        Nw = N - Nm
        i = i + 1
        #print(i<= generation and Nm>0 and Nm>0)
    if(Nw ==0):
        record[i:generation+1] = 1

    return record
        


# In[3]:


def FcA(c,current_fitness,mean_fitness):
    """
    returns the expected ﬁtness of a cell carrying mutant a and (c − 1) other virions.
    """
    return (current_fitness + (c-1)*mean_fitness)/c


# In[4]:


def FA(current_fitness,mean_fitness, NCRatio):
    """
    returns the effective fitness of a mutant a accounting for all possible coinfection scenarios
    """
    fitness = 0
    c = 1
    p_infection=1-poisson.pmf(0,NCRatio)
    while (c<(NCRatio) or c<5 or poisson.pmf(c,NCRatio)>1e-2):
        fitness = fitness + poisson.pmf(c,NCRatio) *FcA(c,current_fitness,mean_fitness)/p_infection
        c = c+1
    return fitness 


# In[5]:


def TestNCRatio(NCRatio,a_fitness):
    initial_a_frequency = 0.4
    frequency_count = np.array([1-initial_a_frequency,initial_a_frequency])
    other_fitness =1
    t = 101 #next twenty generations
    record = []
    for x in np.arange(t)+1:
        record = np.concatenate((record,np.array([frequency_count[1]])))
        mean_fitness = other_fitness*frequency_count[0]+ a_fitness*frequency_count[1]  #mean fitness of the population
        count_a = frequency_count[1] * FA(a_fitness,mean_fitness, NCRatio)
        count_others = frequency_count[0] * FA(other_fitness,mean_fitness,NCRatio)
        norm_factor = count_a + count_others
        frequency_count[1] = count_a / norm_factor
        frequency_count[0] = 1 - frequency_count[1]
    return record


# In[6]:


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
        while (poisson.cdf(k_max,mutant_MOI)<0.99999999):
            k_max = k_max + 1
        l_max = 0
        while (poisson.cdf(l_max,wild_MOI)<0.99999999):
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
        mutant_fitness_sum = mutant_fitness_sum + k * Fkl(sigma_m,k_max,l_max)* 0.00000001 * 0.00000001
        
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


# In[7]:


def NoCoinfection(a_fitness):
    initial_a_frequency = 0.4
    frequency_count = np.array([1-initial_a_frequency,initial_a_frequency])
    other_fitness =1
    t = 101 #next twenty generations
    record = []
    for x in np.arange(t)+1:
        record = np.concatenate((record,np.array([frequency_count[1]])))
        count_a = frequency_count[1] * a_fitness
        count_others = frequency_count[0] * other_fitness
        norm_factor = count_a + count_others
        frequency_count[1] = count_a / norm_factor
        frequency_count[0] = 1 - frequency_count[1]
    return record


# In[8]:


benIllingworth = NoCoinfection(1.1)
ben01 = Main_Deterministic(0.1,np.log(1.1),0.4,100)
ben1 = Main_Deterministic(1,np.log(1.1),0.4,100)
ben5 = Main_Deterministic(5,np.log(1.1),0.4,100)
ben20 = Main_Deterministic(20,np.log(1.1),0.4,100)


# In[9]:


delIllingworth = NoCoinfection(0.9)
del01 = Main_Deterministic(0.1,np.log(0.9),0.4,100)
del1 = Main_Deterministic(1,np.log(0.9),0.4,100)
del5 = Main_Deterministic(5,np.log(0.9),0.4,100)
del20 = Main_Deterministic(20,np.log(0.9),0.4,100)


# In[10]:


ben_det_3 = Main_Deterministic(5,np.log(1.1),0.4,100)
del_det_3 = Main_Deterministic(5,np.log(0.9),0.4,100)


# In[11]:


fig, ax = plt.subplots(figsize=(8*1.2, 6*1.2), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(2,2,1)
plt.plot(np.arange(101)+1,benIllingworth,label="Illingworth")
plt.plot(np.arange(101)+1,ben01,label="0.1")
plt.plot(np.arange(101)+1,ben1,label="1")
plt.plot(np.arange(101)+1,ben5,label="5")
plt.plot(np.arange(101)+1,ben20,label="20")
plt.text(0.5,0.9,"A",fontsize=15)
plt.xlabel("Generation")
plt.ylabel("Variant frequency")
plt.ylim(0,1)


plt.subplot(2,2,2)
plt.plot(np.arange(101)+1,delIllingworth,label="No Coinfection")
plt.plot(np.arange(101)+1,del01,label="0.1")
plt.plot(np.arange(101)+1,del1,label="1")
plt.plot(np.arange(101)+1,del5,label="5")
plt.plot(np.arange(101)+1,del20,label="20")
plt.text(0.5,0.9,"B",fontsize=15)
plt.xlabel("Generation")
plt.ylim(0,1)
#ax.set(ylim=(0, 1))
plt.legend(title = "MOI",loc=1,prop={'size': 8})

plt.subplot(2,2,3)
plt.plot(np.arange(101),ben_det_3,color = "red")
plt.plot(np.arange(101),del_det_3,color = "blue")
for i in np.arange(10):
    ben_sto = stochastic(1000,200,0.4,1.1,100)
    del_sto = stochastic(1000,200,0.4,0.9,100)
    plt.plot(np.arange(101),ben_sto,"--",color = "red",alpha = 0.5)
    plt.plot(np.arange(101),del_sto,"--",color = "blue",alpha = 0.5)
plt.ylim(0,1)
plt.text(0.5,0.9,"C",fontsize=15)
plt.xlabel("Generation")
plt.ylabel("Variant frequency")

plt.subplot(2,2,4)
plt.plot(np.arange(101),ben_det_3,color = "red")
plt.plot(np.arange(101),del_det_3,color = "blue")
for i in np.arange(10):
    ben_sto = stochastic(100,20,0.4,1.1,100)
    del_sto = stochastic(100,20,0.4,0.9,100)
    plt.plot(np.arange(101),ben_sto,"--",color = "red",alpha = 0.5)
    plt.plot(np.arange(101),del_sto,"--",color = "blue",alpha = 0.5)
plt.ylim(0,1)
plt.text(0.5,0.9,"D",fontsize=15)
plt.xlabel("Generation")
plt.show()
fig.savefig("fig1.jpg",bbox_inches = 'tight')

