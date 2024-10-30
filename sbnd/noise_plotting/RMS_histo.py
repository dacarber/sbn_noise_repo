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
import time
import os
from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
import sys


Run_num = input("Enter the Run Number: ")
noise = input("What Noise do you want (raw, coh, or int): ")
directory = f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/run{Run_num}/"
if not os.path.exists(directory):
    os.mkdir(directory)
files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_run{Run_num}.root")
files['tpc_noise;1'].keys()

raw_rms = files['tpc_noise;1'][f'{noise}_rms'].array().to_list()
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
    elif r<11264:
        Noise_df['wire_plane'].append('YA')
Noise_df = pd.DataFrame(Noise_df)

filename = f"RMS_histograms_{Run_num}"
Non_zero_mask = (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000);
print("Average noise of SBND: ",np.mean(Noise_df['Raw_rms'][Non_zero_mask]),"Median noise of SBND: ",np.median(Noise_df['Raw_rms'][Non_zero_mask]))


mask_UB = (Noise_df['wire_plane'] == 'UB') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_UB = np.mean(Noise_df['Raw_rms'][mask_UB])

mask_VB = (Noise_df['wire_plane'] == 'VB') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_VB = np.mean(Noise_df['Raw_rms'][mask_VB])

mask_YB = (Noise_df['wire_plane'] == 'YB') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_YB = np.mean(Noise_df['Raw_rms'][mask_YB])

mask_UA = (Noise_df['wire_plane'] == 'UA') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_UA = np.mean(Noise_df['Raw_rms'][mask_UA])

mask_VA = (Noise_df['wire_plane'] == 'VA') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_VA = np.mean(Noise_df['Raw_rms'][mask_VA])

mask_YA = (Noise_df['wire_plane'] == 'YA') & (Noise_df['Raw_rms'] > 0) & (Noise_df['Raw_rms'] <10000)
mean_YA = np.mean(Noise_df['Raw_rms'][mask_YA])


fig = make_subplots(rows=3,cols=2,column_widths = [0.5,0.5],subplot_titles = (f'<b><br>West 1<sup>st</sup> Induction RMS</b> <br> Mean RMS:{mean_UB:.2f} ADC',f'<b>East 1<sup>st</sup> Induction RMS </b> <br> Mean RMS:{mean_UA:.2f} ADC',f'<b>West 2<sup>nd</sup> Induction RMS </b> <br> Mean RMS:{mean_VB:.2f} ADC',f'<b>East 2<sup>nd</sup> Induction RMS </b> <br> Mean RMS:{mean_VA:.2f} ADC',f'<b> West Collection RMS </b> <br> Mean RMS:{mean_YB:.2f} ADC',f'<b>East Collection RMS</b> <br> Mean RMS:{mean_YA:.2f} ADC',))

fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_UB],marker_color = 'darkgreen',xbins=dict(start = mean_UB - 5,end = mean_UB+5,size=.05)),row = 1, col =2)
fig.update_layout(xaxis = dict(range = [0,4.3]))



fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_VB],marker_color = 'darkorange',xbins=dict(start = mean_VB - 5,end = mean_VB+5,size=.05)),row = 2, col =2)
fig.update_layout(xaxis3 = dict(range = [0,4.3]))



fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_YB],marker_color = 'purple',xbins=dict(start = mean_YB - 5,end = mean_YB + 5,size=.05)),row = 3, col =2)
fig.update_layout(xaxis5 = dict(range = [0,4.3]))




fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_UA],marker_color = 'darkgreen',xbins=dict(start = mean_UA-5,end = mean_UA+5,size=.05)),row = 1, col =1)
fig.update_layout(xaxis2 = dict(range = [0,4.3]))




fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_VA],marker_color = 'darkorange',xbins=dict(start = mean_VA-5,end = mean_VA+5,size=.05)),row = 2, col =1)
fig.update_layout(xaxis4 = dict(range = [0,4.3]))




fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_YA],marker_color = 'purple',xbins=dict(start = mean_YA-5,end = mean_YA+5,size=.05)),row = 3, col =1)
fig.update_layout(xaxis6 = dict(range = [0,4.3]))
if Run_num != 'sim':
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 1, col =1)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 1, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 2, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 2, col =1)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 3, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 115,yshift=157,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 3, col =1)
    #200,250


fig.update_xaxes(title_text = "RMS [ADC]",row = 1, col = 1)
fig.update_xaxes(title_text = "RMS [ADC]",row = 2, col = 1)
fig.update_xaxes(title_text = "RMS [ADC]",row = 3, col = 1)
fig.update_xaxes(title_text = "RMS [ADC]",row = 1, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 2, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 3, col = 2)

fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 1, col = 1)
fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 2, col = 1)
fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 3, col = 1)
fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 1, col = 2)
fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 2, col = 2)
fig.update_yaxes(title_text = "Channels/0.05 ADC",row = 3, col = 2)
#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
#fig.update_layout(yaxis = dict(range = [0,5]),yaxis3 = dict(range = [0,5]),yaxis5 = dict(range = [0,5]))
#fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = 64),xaxis3 = dict(tickmode = 'linear',dtick = 64),xaxis5 = dict(tickmode = 'linear',dtick = 64))
fig.update_layout(title_text='<span style="font-size: 24px;">TPC Noise per Plane</span>', title_x=0.5, title_y=0.98)
fig.update_layout(height = 1000, width = 1000,showlegend = False)

fig.write_image(directory+filename+".png")
fig.write_image(directory+filename+".pdf")
time.sleep(0.2)
fig.write_image(directory+filename+".pdf")
fig.show()

