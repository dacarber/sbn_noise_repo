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



#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_runsim.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_coh_run14784_fft.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_full_14784.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_coh_runsim_fft.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_run18253NF.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_sim_test.root"))


event_length = 3432#3427 #3418

#print((raw_rms[0]))
#print(len((raw_rms)))
#df = {'total_0':[0]*1708,'total_1':[0]*1708,'total_2':[0]*1708}
df= {}
i=0
for f in range(len(files)):
    
    raw_rms = files[f][f'{files[f].keys()[0]}']['raw_FFT'].array()
    df[f'total_U_{f}'] = [0]*(int(event_length/2))
    df[f'total_V_{f}'] = [0]*(int(event_length/2))
    df[f'total_Y_{f}'] = [0]*(int(event_length/2))
    

#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i in tqdm(range(11264)):
        if raw_rms[(int(1716)*(i+1))-2]==0 or np.isnan(raw_rms[(int(1716)*(i+1))-2]):
            df[f'{i}'] = 'skip'
            continue
        df[f'{i}'] =list((raw_rms[i*(int(event_length/2)):(i+1)*(int(event_length/2))]/raw_rms[(int(1716)*(i+1))-2])*(1800/4095))
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
        elif channel<5638:
            if df[f'{channel}'] == 'skip':
                skipped_Y +=1
                continue
            df[f'total_Y_{f}'] =[df[f'total_Y_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_Y_{f}']))]
        elif channel<7622:
            if df[f'{channel}'] == 'skip':
                skipped_U +=1
                continue
            df[f'total_U_{f}'] =[df[f'total_U_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_U_{f}']))]
            print(df[f'total_U_{f}'][1], channel)
        elif channel<9606:
            if df[f'{channel}'] == 'skip':
                skipped_V +=1
                continue
            df[f'total_V_{f}'] =[df[f'total_V_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_V_{f}']))]
        elif channel<11276:
            if df[f'{channel}'] == 'skip':
                skipped_Y +=1
                continue
            df[f'total_Y_{f}'] =[df[f'total_Y_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_Y_{f}']))]

df_orig= {}
i_orig=0
for f in range(len(files)):
    
    orig_rms = files[f][f'{files[f].keys()[0]}']['orig_FFT'].array()
    df_orig[f'total_U_{f}'] = [0]*(int(event_length/2))
    df_orig[f'total_V_{f}'] = [0]*(int(event_length/2))
    df_orig[f'total_Y_{f}'] = [0]*(int(event_length/2))
    

#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i_orig in tqdm(range(11264)):
        if orig_rms[(int(1716)*(i_orig+1))-2]==0:
            df_orig[f'{i_orig}'] = 'skip'
            continue
        df_orig[f'{i_orig}'] =list((orig_rms[i_orig*(int(event_length/2)):(i_orig+1)*(int(event_length/2))]/orig_rms[(int(1716)*(i_orig+1))-2])*(1800/4095))
    skipped_U_orig = 0
    skipped_Y_orig = 0
    skipped_V_orig = 0
    for channel in tqdm(range(0,11264,1)):
        if channel <1984:
            if df_orig[f'{channel}'] == 'skip':
                skipped_U_orig +=1
                continue
            df_orig[f'total_U_{f}'] =[df_orig[f'total_U_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_U_{f}']))]
        elif channel<3968:
            if df_orig[f'{channel}'] == 'skip':
                skipped_V_orig +=1
                continue
            df_orig[f'total_V_{f}'] =[df_orig[f'total_V_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_V_{f}']))]
        elif channel<5638:
            if df_orig[f'{channel}'] == 'skip':
                skipped_Y_orig +=1
                continue
            df_orig[f'total_Y_{f}'] =[df_orig[f'total_Y_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_Y_{f}']))]
        elif channel<7622:
            if df_orig[f'{channel}'] == 'skip':
                skipped_U_orig +=1
                continue
            df_orig[f'total_U_{f}'] =[df_orig[f'total_U_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_U_{f}']))]
        elif channel<9606:
            if df_orig[f'{channel}'] == 'skip':
                skipped_V_orig +=1
                continue
            df_orig[f'total_V_{f}'] =[df_orig[f'total_V_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_V_{f}']))]
        elif channel<11276:
            if df_orig[f'{channel}'] == 'skip':
                skipped_Y_orig +=1
                continue
            df_orig[f'total_Y_{f}'] =[df_orig[f'total_Y_{f}'][j] + df_orig[f'{channel}'][j] for j in range(len(df_orig[f'total_Y_{f}']))]
        
freq = list(range(len(df['100'])))
freq = (np.add(freq,.5))*2/(int(event_length))
color = ['blue','green','orange','purple', 'black','red']
print(np.divide(df[f'total_U_{f}'],1670*2-skipped_U))
fig_u = make_subplots(rows=1,cols =1,subplot_titles = (f"<span style='font-size: 26px;'>FFT Spectrum U Plane {Run_num}</span>",))
fig_u.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_U_{f}'],1984*2-skipped_U),marker_color = color[0],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>First Induction NF</span>"),row = 1, col = 1)
fig_u.add_trace(go.Scatter(x=freq,y = np.divide(df_orig[f'total_U_{f}'],1984*2-skipped_U_orig),marker_color = color[3],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>First Induction</span>"),row = 1, col = 1)

fig_u.update_xaxes(title_text = "<span style='font-size: 20px;'>Frequency [MHz]</span>",row = 1, col = 1)
fig_u.update_yaxes(title_text = "<span style='font-size: 20px;'>Magnitude [mV/0.59 kHz]</span>",range = (0,100),row = 1, col = 1)
fig_u.update_layout(xaxis = dict(tickmode = 'linear',dtick = .1))
fig_u.update_layout(height = 600, width = 1800,showlegend = True)
fig_u.show()


fig_v = make_subplots(rows=1,cols =1,subplot_titles = (f"<span style='font-size: 26px;'>FFT Spectrum V Plane {Run_num}</span>",))
fig_v.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_V_{f}'],1984*2-skipped_V),marker_color = color[1],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>Second Induction NF</span>"),row = 1, col = 1)
fig_v.add_trace(go.Scatter(x=freq,y = np.divide(df_orig[f'total_V_{f}'],1984*2-skipped_V_orig),marker_color = color[4],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>Second Induction</span>"),row = 1, col = 1)

fig_v.update_xaxes(title_text = "<span style='font-size: 20px;'>Frequency [MHz]</span>",row = 1, col = 1)
fig_v.update_yaxes(title_text = "<span style='font-size: 20px;'>Magnitude [mV/0.59 kHz]</span>",range = (0,100),row = 1, col = 1)
fig_v.update_layout(xaxis = dict(tickmode = 'linear',dtick = .1))
fig_v.update_layout(height = 600, width = 1800,showlegend = True)
fig_v.show()
print(np.divide(df[f'total_Y_{f}'],1670*2-skipped_Y))

fig_w = make_subplots(rows=1,cols =1,subplot_titles = (f"<span style='font-size: 26px;'>FFT Spectrum W Plane {Run_num}</span>",))

fig_w.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_Y_{f}'],1670*2-skipped_Y),marker_color = color[2],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>Collection NF</span>"),row = 1, col = 1)
fig_w.add_trace(go.Scatter(x=freq,y = np.divide(df_orig[f'total_Y_{f}'],1670*2-skipped_Y_orig),marker_color = color[5],opacity = 1/(f+1),name = f"<span style='font-size: 16px;'>Collection</span>"),row = 1, col = 1)

fig_w.update_xaxes(title_text = "<span style='font-size: 20px;'>Frequency [MHz]</span>",row = 1, col = 1)
fig_w.update_yaxes(title_text = "<span style='font-size: 20px;'>Magnitude [mV/0.59 kHz]</span>",range = (0,100),row = 1, col = 1)
fig_w.update_layout(xaxis = dict(tickmode = 'linear',dtick = .1))
fig_w.update_layout(height = 600, width = 1800,showlegend = True)
fig_w.show()