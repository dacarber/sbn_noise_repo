###################################################
#### Made by Dan Carber on 08/19/2024          ####               
#### This script will make plots of FFTs where ####
#### the wire planes are on the same plot      ####
###################################################



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
import pandas as pd
from datetime import datetime, timedelta
from tqdm import tqdm
import os
from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
import sys
from operator import add

files = []
Run_num = input("Enter the Run Number: ")

fig = make_subplots(rows=1,cols =1,subplot_titles = (f'FFT Spectrum {Run_num}',))

#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_runsim.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_coh_run14784_fft.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_coh_runsim_fft.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_run14784_20.root"))





#print((raw_rms[0]))
#print(len((raw_rms)))
#df = {'total_0':[0]*1708,'total_1':[0]*1708,'total_2':[0]*1708}
df= {}
for f in range(len(files)):
    raw_rms = files[f][f'{files[f].keys()[0]}']['coh_FFT'].array()
    df[f'total_U_{f}'] = [0]*1709
    df[f'total_V_{f}'] = [0]*1709
    df[f'total_Y_{f}'] = [0]*1709
    

#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i in tqdm(range(11264)):
        if raw_rms[1709*(i+1)-1]==0:
            df[f'{i}'] = 'skip'
            continue
        df[f'{i}'] =list(raw_rms[i*1709:(i+1)*1709]/raw_rms[1709*(i+1)-1])
    skipped_U = 0
    skipped_Y = 0
    skipped_V = 0
    
    for channel in tqdm(range(0,11264,1)):
        if channel <1984:
            if df[f'{channel}'] == 'skip':
                skipped_U +=1
                continue
            df[f'total_U_{f}'] =[df[f'total_U_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_U_{f}']))]
        elif channel<3968:
            if df[f'{channel}'] == 'skip':
                skipped_V +=1
                continue
            df[f'total_V_{f}'] =[df[f'total_V_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_V_{f}']))]
        elif channel<5632:
            if df[f'{channel}'] == 'skip':
                skipped_Y +=1
                continue
            df[f'total_Y_{f}'] =[df[f'total_Y_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_Y_{f}']))]
        elif channel<7616:
            if df[f'{channel}'] == 'skip':
                skipped_U +=1
                continue
            df[f'total_U_{f}'] =[df[f'total_U_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_U_{f}']))]
        elif channel<9600:
            if df[f'{channel}'] == 'skip':
                skipped_V +=1
                continue
            df[f'total_V_{f}'] =[df[f'total_V_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_V_{f}']))]
        elif channel<11264:
            if df[f'{channel}'] == 'skip':
                skipped_Y +=1
                continue
            df[f'total_Y_{f}'] =[df[f'total_Y_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_Y_{f}']))]
        

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['blue','green','orange']
    
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_U_{f}'],1984*2-skipped_U),marker_color = color[0],opacity = 1/(f+1),name = f'First Induction'),row = 1, col = 1)
    
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_V_{f}'],1984*2-skipped_V),marker_color = color[1],opacity = 1/(f+1),name = f'Second Induction'),row = 1, col = 1)
   
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_Y_{f}'],1664*2-skipped_Y),marker_color = color[2],opacity = 1/(f+1),name = f'Collection'),row = 1, col = 1)

fig.update_xaxes(title_text = "Frequency [MHz]",row = 1, col = 1)
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = .01))
fig.update_layout(height = 600, width = 1800,showlegend = True)
fig.show()
