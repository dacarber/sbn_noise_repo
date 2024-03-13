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

import os
from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
import sys


Run_num = int(input("Enter the Run Number: "))
directory = f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/run{Run_num}/"
if not os.path.exists(directory):
    os.mkdir(directory)
files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_run{Run_num}.root")
files['tpc_noise;1'].keys()

raw_rms = files['tpc_noise;1']['raw_rms'].array().to_list()
Noise_df = {'Channel_id':[],'Raw_rms':[],'wire_plane':[]}
for r,rms in enumerate(raw_rms):
    Noise_df['Channel_id'].append(r)
    Noise_df['Raw_rms'].append(rms)
    if r <1984:
        Noise_df['wire_plane'].append('UB')
    elif r<3968:
        Noise_df['wire_plane'].append('VB')
    elif r<5632:
        Noise_df['wire_plane'].append('YB')
    elif r<7616:
        Noise_df['wire_plane'].append('UA')
    elif r<9600:
        Noise_df['wire_plane'].append('VA')
    else:
        Noise_df['wire_plane'].append('YA')
Noise_df = pd.DataFrame(Noise_df)

filename = f"RMS_plots_planeB_{Run_num}.png"

fig = make_subplots(rows=3,cols=2,column_widths = [0.7,0.3],subplot_titles = (f'UB Channel RMS Run {Run_num}','UB RMS',f'VB Channel RMS Run {Run_num}','VB RMS',f'YB Channel RMS Run {Run_num}','YB RMS',))
mask = Noise_df['wire_plane'] == 'UB'
median = np.median(Noise_df['Raw_rms'][mask])
mean = np.mean(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'red',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 1, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'red'),row = 1, col = 1)
fig.update_layout(xaxis2 = dict(range = [0,median+5]))
fig.update_layout(margin = dict(r=200))
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=120,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=110,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)


mask = Noise_df['wire_plane'] == 'VB'
median = np.median(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'purple',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 2, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 2, col = 1)
fig.update_layout(xaxis4 = dict(range = [0,median+5]))

fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-120,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-130,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

mask = Noise_df['wire_plane'] == 'YB'
median = np.median(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'blue',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 3, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'blue'),row = 3, col = 1)
fig.update_layout(xaxis6 = dict(range = [0,median+5]))

fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-370,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-380,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 1, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 2, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 3, col = 1)
fig.update_xaxes(title_text = "RMS [ADC]",row = 1, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 2, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 3, col = 2)
#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
fig.update_layout(yaxis = dict(range = [0,5]),yaxis3 = dict(range = [0,5]),yaxis5 = dict(range = [0,5]))
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = 64),xaxis3 = dict(tickmode = 'linear',dtick = 64),xaxis5 = dict(tickmode = 'linear',dtick = 64))
fig.update_layout(height = 800, width = 1200,showlegend = False)

fig.write_image(directory+filename)
fig.show()

filename = f"RMS_plots_planeA_{Run_num}.png"

fig = make_subplots(rows=3,cols=2,column_widths = [0.7,0.3],subplot_titles = (f'UA Channel RMS Run {Run_num}','UA RMS',f'VA Channel RMS Run {Run_num}','VA RMS',f'YA Channel RMS Run {Run_num}','YA RMS',))
mask = Noise_df['wire_plane'] == 'UA'
median = np.median(Noise_df['Raw_rms'][mask])
mean = np.mean(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'red',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 1, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'red'),row = 1, col = 1)
fig.update_layout(xaxis2 = dict(range = [0,median+5]))
fig.update_layout(margin = dict(r=200))
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=60,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=50,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)


mask = Noise_df['wire_plane'] == 'VA'
median = np.median(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'purple',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 2, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 2, col = 1)
fig.update_layout(xaxis4 = dict(range = [0,median+5]))

fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-180,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-190,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

mask = Noise_df['wire_plane'] == 'YA'
median = np.median(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask],marker_color = 'blue',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 3, col =2)
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'blue'),row = 3, col = 1)
fig.update_layout(xaxis6 = dict(range = [0,median+5]))

fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-420,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=-430,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 1, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 2, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 3, col = 1)
fig.update_xaxes(title_text = "RMS [ADC]",row = 1, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 2, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 3, col = 2)
#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
fig.update_layout(yaxis = dict(range = [0,5]),yaxis3 = dict(range = [0,5]),yaxis5 = dict(range = [0,5]))
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = 64),xaxis3 = dict(tickmode = 'linear',dtick = 64),xaxis5 = dict(tickmode = 'linear',dtick = 64))
fig.update_layout(height = 800, width = 1200,showlegend = False)

fig.write_image(directory+filename)
fig.show()
