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
fig = make_subplots(rows=3,cols =1,subplot_titles = ('East First Induction FFT Spectrum','East Second Induction FFT Spectrum','East Collection FFT Spectrum'))

files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/fft_output_runsim.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_run11995.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_run14784_signal.root"))
#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_run14784_20.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_sim_test.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_full_14784.root"))
run_number = ["Old Sim","New Simulated","Run 14784"]



#print((raw_rms[0]))
#print(len((raw_rms)))
#df = {'total_0':[0]*1708,'total_1':[0]*1708,'total_2':[0]*1708}
df= {}
for f in range(len(files)):
    raw_rms = files[f][f'{files[f].keys()[0]}']['avg_FFT'].array()
    df[f'total_UB_{f}'] = [0]*1709
    df[f'total_VB_{f}'] = [0]*1709
    df[f'total_YB_{f}'] = [0]*1709
    

#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i in tqdm(range(11264)):
        if raw_rms[1709*(i+1)-1]==0:
            df[f'{i}'] = 'skip'
            continue
        df[f'{i}'] =list(raw_rms[i*1709:(i+1)*1709]/raw_rms[1709*(i+1)-1])
    skipped = 0
    for channel in tqdm(range(0,1984,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        #print(len(df[f'total_UB_{f}']))
        #print(len(df[f'{channel}']))
        df[f'total_UB_{f}'] =[df[f'total_UB_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'total_UB_{f}']))]
        #print(df[f'total_UB_{f}'])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_UB_{f}'],1984-skipped),marker_color = color[f],opacity = 1/(f+1),name = f'{run_number[f]}'),row = 1, col = 1)
    
    skipped = 0
    for channel in tqdm(range(1984,3968,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        df[f'total_VB_{f}'] =[df[f'total_VB_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'{channel}']))]

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_VB_{f}'],1984-skipped),marker_color = color[f],opacity = 1/(f+1),showlegend=False),row = 2, col = 1)
    
    skipped = 0
    for channel in tqdm(range(3968,5632,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        df[f'total_YB_{f}'] =[df[f'total_YB_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'{channel}']))]

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_YB_{f}'],1664-skipped),marker_color = color[f],opacity = 1/(f+1),showlegend=False),row = 3, col = 1)
#fig.update_layout(xaxis = dict(range = [1*0,.996*1]))
fig.update_xaxes(title_text = "Frequency [MHz]",row = 3, col = 1)
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = .01),xaxis2 = dict(tickmode = 'linear',dtick = .01),xaxis3 = dict(tickmode = 'linear',dtick = .01))
fig.update_layout(height = 900, width = 1500,showlegend = True)
fig.show()






fig = make_subplots(rows=3,cols =1,subplot_titles = ('West First Induction FFT Spectrum','West Second Induction FFT Spectrum','West Collection FFT Spectrum'))
df= {}
for f in range(len(files)):
    raw_rms = files[f][f'{files[f].keys()[0]}']['avg_FFT'].array()
    df[f'total_UA_{f}'] = [0]*1709
    df[f'total_VA_{f}'] = [0]*1709
    df[f'total_YA_{f}'] = [0]*1709
    
#print(df['total'])
    channel = -1
    for i in tqdm(range(11264)):
        if raw_rms[1709*(i+1)-1]==0:
            df[f'{i}'] = 'skip'
            continue
        df[f'{i}'] =list(raw_rms[i*1709:(i+1)*1709]/raw_rms[1709*(i+1)-1])
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    skipped = 0
    for channel in tqdm(range(5632,7616,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        df[f'total_UA_{f}'] =[df[f'total_UA_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'{channel}']))]

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_UA_{f}'],1984-skipped),marker_color = color[f],opacity = 1/(f+1),name = f'Run {run_number[f]}'),row = 1, col = 1)
    
    skipped = 0
    for channel in tqdm(range(7616,9600,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        df[f'total_VA_{f}'] =[df[f'total_VA_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'{channel}']))]

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_VA_{f}'],1984-skipped),marker_color = color[f],opacity = 1/(f+1),showlegend=False),row = 2, col = 1)
    
    skipped = 0
    for channel in tqdm(range(9600,11264,1)):
        if df[f'{channel}'] == 'skip':
            skipped +=1
            continue
        df[f'total_YA_{f}'] =[df[f'total_YA_{f}'][j] + df[f'{channel}'][j] for j in range(len(df[f'{channel}']))]
        
    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_YA_{f}'],1664-skipped),marker_color = color[f],opacity = 1/(f+1),showlegend=False),row = 3, col = 1)

#fig.update_layout(xaxis = dict(range = [1*0,.996*1]))
fig.update_xaxes(title_text = "Frequency [MHz]",row = 3, col = 1)
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = .01),xaxis2 = dict(tickmode = 'linear',dtick = .01),xaxis3 = dict(tickmode = 'linear',dtick = .01))
fig.update_layout(height = 900, width = 1500,showlegend = True)
fig.show()
#mask = Noise_df['wire_plane'] == 'UB'
#median = np.median(Noise_df['Raw_rms'][mask])
#mean = np.mean(Noise_df['Raw_rms'][mask])
#fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'red',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 1, col =2)
#fig.add_trace(go.Scatter(x=freq,y = np.divide(df['total_0'],1000),marker_color = 'red'),row = 1, col = 1)
#fig.add_trace(go.Scatter(x=freq,y = np.divide(df['total_1'],1000)-30,marker_color = 'green'),row = 1, col = 1)
#fig.add_trace(go.Scatter(x=freq,y = np.divide(df['total_2'],1000)-60,marker_color = 'blue'),row = 1, col = 1)


#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=120,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=110,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

