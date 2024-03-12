#!/usr/bin/env python
# coding: utf-8

# In[25]:


#import ROOT
import uproot
import matplotlib.pyplot as plt
import seaborn
import numpy as np
import math
from scipy.fft import fft, fftfreq
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import numpy as np
#from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
pyo.init_notebook_mode()
import time
from tqdm import tqdm
from multiprocessing import Pool

file = uproot.open("/Users/danielcarber/Documents/SBND/Noise Analysis/hitdumper/gen_g4_detsim_reco1-61193afb-d94b-460d-a307-3dd6417a9b80_HitDumper-20230523T200522.root")


# In[13]:


file['Events;1'].keys()


# In[14]:


ADC = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fADC']
Channel = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fChannel']
Samples = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fSamples']
Pedestal = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fPedestal']
Sigma = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fSigma']
Compression = file['Events;1']['raw::RawDigits_daq__DetSim./raw::RawDigits_daq__DetSim.obj/raw::RawDigits_daq__DetSim.obj.fCompression']


# In[15]:


#ADC.all_members)
ADC_array = ADC.array()


# In[16]:


Pedestal_array = Pedestal.array()


# In[17]:


print(np.shape(ADC_array))
channel = np.array(range(11224))
print(channel)
#np.shape(Pedestal_array)


# In[18]:


def signal_removal(channel):
    good_channels = []
    for r in range(34):
        max_sig = max(channel[r*100:(r+1)*100])
        min_sig = min(channel[r*100:(r+1)*100])
        #print(max_sig,min_sig)
        signal_diff= (max_sig - min_sig)
        if signal_diff >30:
            good_channels.extend([0]*100)
            continue
        #print(channel[r*100:(r+1)*100])
        good_channels.extend(channel[r*100:(r+1)*100])
    #print(good_channels)
    return good_channels


# In[75]:


def correlation_func(event,i = 0):
    correlations = np.zeros((11224,11224))
    channels_y = []
    for c_y,channel_y in enumerate((event)):
            #print(i, c_y)
            channel_y = np.subtract(channel_y,Pedestal_array[i][c_y])
            noise_channels_y = signal_removal(channel_y)
            channels_y.extend(noise_channels_y)
            #r =  np.corrcoef(noise_channels_x, noise_channels_y)
            #print(r)
    for c_x, channel_x in enumerate(tqdm(event)):
        #print(np.shape(Pedestal_array[i][c_x]))
        #print(np.shape(channel_x))
        channel_x = np.subtract(channel_x,Pedestal_array[i][c_x])
        noise_channels_x = signal_removal(channel_x)
        #print(channel_x)
        
        #print(r)
        for j,channel_y in enumerate((event)):
        #    #print(i, c_y)
        #    channel_y = np.subtract(channel_y,Pedestal_array[i][c_y])
        #    noise_channels_y = signal_removal(channel_y)
        #    #r =  np.corrcoef(noise_channels_x, noise_channels_y)
        #    #print(r)
            
            a, b, c, d, e = 0,0,0,0,0
            
            #for time_tick in (range(len(noise_channels_x))):
            #    a += noise_channels_x[time_tick]*noise_channels_y[time_tick]
            #    c += noise_channels_x[time_tick]**2
            #    e += noise_channels_y[time_tick]**2
            a = np.sum(np.diag(np.outer(noise_channels_x,channels_y[j])))
            c = np.sum(np.square(noise_channels_x))
            e = np.sum(np.square(channels_y[j]))
            #b = np.sum(noise_channels_x)
            #d = np.sum(channel_y[i])
            
            #correlations[c_x][c_y] = r[0,1]
            #correlations[c_y][c_x] = r[1,0]
    #print(r)
            #correlations[c_x][c_y] += (3400*a - b*d)/(math.sqrt((3400*c-(b**2))*(3400*e - (d**2))))   
            correlations[c_x][c_y] += (a)/(math.sqrt(c)*math.sqrt(e))  
    return correlations


# In[77]:


def main():
    correlation = np.zeros((11224,11224))
    #for i, event in enumerate(ADC_array):
    #print(i)
    start = time.time()
    #cor = correlation_func(event,i)
    with Pool(5) as p:
        test_1 = set(p.map(correlation_func,ADC_array[0:4]))

    end = time.time()
    print("elapsed",end-start)
    print(cor)
    correlation = np.divide(correlation, len(ADC_array))
    print(len(ADC_array))


# In[78]:


if __name__ == '__main__':
    main()


# In[ ]:




