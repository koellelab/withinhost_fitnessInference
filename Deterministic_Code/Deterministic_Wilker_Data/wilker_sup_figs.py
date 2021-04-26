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


#df = pd.read_csv('wilkerMCMCnoise25.csv')
#df = pd.read_csv('wilkerMCMCnoise400.csv')
#df = pd.read_csv('wilkerMCMC6h.csv')
df = pd.read_csv('wilkerMCMC12h.csv')


# In[3]:


df2 = pd.read_csv('wilkerMCMC.csv')


# In[4]:


init_freq1 = df[['initfreq1']].values
init_freq2 = df[['initfreq2']].values
init_freq3 = df[['initfreq3']].values
init_freq4 = df[['initfreq4']].values
logmoi = df[['logmoi']].values
logfitness = df[['logfitness']].values

logmoi2 = df2[['logmoi']].values
logfitness2 = df2[['logfitness']].values

convertmoi = np.exp(logmoi)
convertfitness = np.exp(logfitness)

convertmoi2 = np.exp(logmoi2)
convertfitness2 = np.exp(logfitness2)


# In[5]:


interval = np.arange(2000,10001,50)
sample = np.zeros(len(interval)-1)
for i in np.arange(len(interval)-1):
    sample[i] = random.randint(interval[i],interval[i+1])
sample = sample.astype(int)


# In[20]:


fig, ax = plt.subplots(figsize=(15, 3), dpi= 180, facecolor='w', edgecolor='k')
plt.subplot(1,3,1)
plt.hist(convertmoi[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertmoi[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(convertmoi[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(convertmoi[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("MOI")
#plt.text(10.7,0.22,"A",fontsize=15)
plt.text(5,0.165,"A",fontsize=15)
plt.ylabel("Proportion")


plt.subplot(1,3,2)
plt.hist(convertfitness[sample],weights=np.ones(len(sample)) / len(sample))
plt.axvline(np.quantile(convertfitness[sample],0.5), color="black",alpha = 0.8)
plt.axvline(np.quantile(convertfitness[sample],0.025), color="black",alpha = 0.8, linestyle='--')
plt.axvline(np.quantile(convertfitness[sample],0.975), color="black",alpha = 0.8, linestyle='--')
plt.xlabel("Variant fitness")
#plt.text(11.8,0.395,"B",fontsize=15)
#plt.text(3.67,0.203,"B",fontsize=15)
#plt.text(4.3,0.2,"B",fontsize=15)
plt.text(9.5,0.178,"B",fontsize=15)
plt.ylabel("Proportion")

plt.subplot(1,3,3)
plt.plot(convertfitness2,convertmoi2,".",alpha=0.02,color="red")
plt.plot(convertfitness,convertmoi,".")  #MOI; Fitness
plt.xlabel("Variant fitness")
plt.ylabel("MOI")
#plt.text(13,2.4,"C",fontsize=15)
#plt.text(9.4,2.1,"C",fontsize=15)
#plt.text(9.4,2.1,"C",fontsize=15)
plt.text(9.4,2.1,"C",fontsize=15)
plt.gca().invert_yaxis()
plt.show()

fig.savefig("figS7.jpg",bbox_inches = 'tight')


# In[ ]:




