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

files = []
fig = make_subplots(rows=3,cols =1,subplot_titles = ('UB FFT Spectrum','VB FFT Spectrum','YB FFT Spectrum'))

#files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_run11888.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_11995.root"))
files.append(uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/fft_output_11541.root"))




#print((raw_rms[0]))
#print(len((raw_rms)))
#df = {'total_0':[0]*1708,'total_1':[0]*1708,'total_2':[0]*1708}
df= {}
for f in range(len(files)):
    raw_rms = files[f]['tpc_noise;3']['avg_FFT'].array()
    df[f'total_UB_{f}'] = [0]*1708
    df[f'total_VB_{f}'] = [0]*1708
    df[f'total_YB_{f}'] = [0]*1708
#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i in tqdm(range(1708*0,1708*1984,1)):
        if i%1708 == 0:
            #print(channel)
            channel +=1
            df[f'{channel}'] = [raw_rms[i]]  
            df[f'total_UB_{f}'][0] +=raw_rms[i]

        else:
            df[f'{channel}'].append(raw_rms[i])
            df[f'total_UB_{f}'][i%1708] +=raw_rms[i]
        #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_UB_{f}'],1984),marker_color = color[f],opacity = 1/(f+1)),row = 1, col = 1)

    for i in tqdm(range(1708*1984,1708*3968,1)):
            if i%1708 == 0:
                #print(channel)
                channel +=1
                df[f'{channel}'] = [raw_rms[i]]  
                df[f'total_VB_{f}'][0] +=raw_rms[i]

            else:
                df[f'{channel}'].append(raw_rms[i])
                df[f'total_VB_{f}'][i%1708] +=raw_rms[i]
            #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_VB_{f}'],1984),marker_color = color[f],opacity = 1/(f+1)),row = 2, col = 1)

    for i in tqdm(range(1708*3968,1708*5632,1)):
            if i%1708 == 0:
                #print(channel)
                channel +=1
                df[f'{channel}'] = [raw_rms[i]]  
                df[f'total_YB_{f}'][0] +=raw_rms[i]

            else:
                df[f'{channel}'].append(raw_rms[i])
                df[f'total_YB_{f}'][i%1708] +=raw_rms[i]
            #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_YB_{f}'],1664),marker_color = color[f],opacity = 1/(f+1)),row = 3, col = 1)
#fig.update_layout(xaxis = dict(range = [1*0,.996*1]))
fig.update_xaxes(title_text = "Frequency [MHz]",row = 3, col = 1)
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = .01),xaxis2 = dict(tickmode = 'linear',dtick = .01),xaxis3 = dict(tickmode = 'linear',dtick = .01))
fig.update_layout(height = 600, width = 1800,showlegend = False)
fig.show()

fig = make_subplots(rows=3,cols =1,subplot_titles = ('UA FFT Spectrum','VA FFT Spectrum','YA FFT Spectrum'))
df= {}
for f in range(len(files)):
    raw_rms = files[f]['tpc_noise;3']['avg_FFT'].array()
    df[f'total_UA_{f}'] = [0]*1708
    df[f'total_VA_{f}'] = [0]*1708
    df[f'total_YA_{f}'] = [0]*1708
#print(df['total'])
    channel = -1
    #for i in tqdm(range(len(raw_rms))):
    #0-1984 UB, 1984-3968 VB, 3968-5632 YB, 5632-7616 UA, 7616-9600 VA, 9600-11264 YA
    for i in tqdm(range(1708*5632,1708*7616,1)):
        if i%1708 == 0:
            #print(channel)
            channel +=1
            df[f'{channel}'] = [raw_rms[i]]  
            df[f'total_UA_{f}'][0] +=raw_rms[i]

        else:
            df[f'{channel}'].append(raw_rms[i])
            df[f'total_UA_{f}'][i%1708] +=raw_rms[i]
        #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_UA_{f}'],1984),marker_color = color[f],opacity = 1/(f+1)),row = 1, col = 1)

    for i in tqdm(range(1708*7616,1708*9600,1)):
            if i%1708 == 0:
                #print(channel)
                channel +=1
                df[f'{channel}'] = [raw_rms[i]]  
                df[f'total_VA_{f}'][0] +=raw_rms[i]

            else:
                df[f'{channel}'].append(raw_rms[i])
                df[f'total_VA_{f}'][i%1708] +=raw_rms[i]
            #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_VA_{f}'],1984),marker_color = color[f],opacity = 1/(f+1)),row = 2, col = 1)

    for i in tqdm(range(1708*9600,1708*11264,1)):
            if i%1708 == 0:
                #print(channel)
                channel +=1
                df[f'{channel}'] = [raw_rms[i]]  
                df[f'total_YA_{f}'][0] +=raw_rms[i]

            else:
                df[f'{channel}'].append(raw_rms[i])
                df[f'total_YA_{f}'][i%1708] +=raw_rms[i]
            #print(raw_rms[i])

    freq = list(range(len(df['0'])))
    freq = (np.add(freq,.5))*2/3415
    color = ['red','green','blue']
    fig.add_trace(go.Scatter(x=freq,y = np.divide(df[f'total_YA_{f}'],1664),marker_color = color[f],opacity = 1/(f+1)),row = 3, col = 1)

#fig.update_layout(xaxis = dict(range = [1*0,.996*1]))
fig.update_xaxes(title_text = "Frequency [MHz]",row = 3, col = 1)
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = .01),xaxis2 = dict(tickmode = 'linear',dtick = .01),xaxis3 = dict(tickmode = 'linear',dtick = .01))
fig.update_layout(height = 600, width = 1800,showlegend = False)
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

