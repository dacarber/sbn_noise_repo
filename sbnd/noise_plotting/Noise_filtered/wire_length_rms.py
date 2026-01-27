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
import os
from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
import sys

Run_num = input("Enter the Run Number: ")
noise = input("What Noise do you want (raw, orig, or int): ")

files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_run{Run_num}.root")
#files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/noise_output_coh_run{Run_num}.root")

files['tpc_noise;5'].keys()
bad_wires = [4800,4801,4802,4803,4804,4805,10438,10439,10440,10441,10442,104443]
bad_count = 0
raw_rms = files['tpc_noise;5'][f'{noise}_rms'].array().to_list()
Noise_df = {'Channel_id':[],'Raw_rms':[],'wire_plane':[]}
for r,rms in enumerate(raw_rms): 
    if r in bad_wires:
        print("Bad",r)
        bad_count +=1
        continue
    if r == 0:
        continue
        #print(r-bad_count,rms)
    Noise_df['Channel_id'].append(r-bad_count)
    print(rms)
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


Argon_df = pd.read_excel("/Users/danielcarber/Documents/SBND/Noise_Analysis/sbn_noise_repo/sbnd/datafiles/SBND_LAr_level.xlsx")
print(Argon_df['Lar_level'][27])
Argon_df['Lar_level'] = Argon_df['Lar_level']/Argon_df['Lar_level'][28]*100
#print(datetime.fromtimestamp(Argon_df['Time'][0]))


wire_plane_list = ['UB','VB','YB','UA','VA','YA']
wire_df = {'Channel_id':[],'cryo':[],'tpc':[],'tpc':[],'plane':[],'rel_wire':[],'x_0':[],'y_0':[],'z_0':[],'x_1':[],'y_1':[],'z_1':[],'r':[]}
wire_txt = '/Users/danielcarber/Documents/SBND/Noise_Analysis/sbn_noise_repo/sbnd/datafiles/Wire_lengths.txt'

with open(wire_txt) as f:
    for line in f:
        #print(line)
        currentline = line.split(" ")
        #print(currentline)
        wire_df['Channel_id'].append(int(currentline[0]))
        wire_df['cryo'].append(int(currentline[1]))
        wire_df['tpc'].append(int(currentline[2]))
        wire_df['plane'].append(int(currentline[3]))
        wire_df['rel_wire'].append(int(currentline[4]))
        wire_df['x_0'].append(float(currentline[5]))
        wire_df['y_0'].append(float(currentline[6]))
        wire_df['z_0'].append(float(currentline[7]))
        wire_df['x_1'].append(float(currentline[8]))
        wire_df['y_1'].append(float(currentline[9]))
        z_1 = currentline[10]
        #print(z_1[:-2])
        wire_df['z_1'].append(float(currentline[10][:-2]))
        length = np.sqrt(np.square(float(currentline[8])-float(currentline[5]))+np.square(float(currentline[9])-float(currentline[6]))+np.square(float(currentline[10][:-2])-float(currentline[7])))
        wire_df['r'].append(length)
wire_df = pd.DataFrame(wire_df)        

Noise_df =Noise_df.merge(wire_df, on='Channel_id')
print(Noise_df)
print(Noise_df[Noise_df['plane']==2])



mask_U = (Noise_df['r']>576.5) & (Noise_df['plane']==0) & (Noise_df['Raw_rms']>1.1)
mean_U = np.mean(Noise_df['Raw_rms'][mask_U])

mask_V = (Noise_df['r']>576.5) & (Noise_df['plane']==1) & (Noise_df['Raw_rms']>1.1)
mean_V = np.mean(Noise_df['Raw_rms'][mask_V])

mask_Y =  (Noise_df['plane']==2) & (Noise_df['Raw_rms']>1.4)
mean_Y = np.mean(Noise_df['Raw_rms'][mask_Y])

mask_UA = Noise_df['wire_plane'] == 'UA'
mean_UA = np.mean(Noise_df['Raw_rms'][mask_UA])

mask_VA = Noise_df['wire_plane'] == 'VA'
mean_VA = np.mean(Noise_df['Raw_rms'][mask_VA])

mask_YA = Noise_df['wire_plane'] == 'YA'
mean_YA = np.mean(Noise_df['Raw_rms'][mask_YA])


def Extract(lst,pos):
    return [item[pos] for item in lst]
#print(Times)


x = np.array(Noise_df['r'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1][Noise_df['r']<287])
y = np.array(Noise_df['Raw_rms'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1][Noise_df['r']<287])
#a_00, b_00 = np.polyfit(x[:830],y[:830], 1)

x = np.array(Noise_df['r'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1][Noise_df['r']>287][Noise_df['r']<576])
y = np.array(Noise_df['Raw_rms'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1][Noise_df['r']>287][Noise_df['r']<576])
#a_01, b_01 = np.polyfit(x[:700],y[:700], 1)
#print(a_01)
x = np.array(Noise_df['r'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1][Noise_df['r']<287])
y = np.array(Noise_df['Raw_rms'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1][Noise_df['r']<287])
#a_10, b_10 = np.polyfit(x[:830],y[:830], 1)

x = np.array(Noise_df['r'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1][Noise_df['r']>287][Noise_df['r']<576])
y = np.array(Noise_df['Raw_rms'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1][Noise_df['r']>287][Noise_df['r']<576])
#a_11, b_11 = np.polyfit(x[:830],y[:830], 1)



fig = make_subplots(rows=3,cols=2,specs=[[{"rowspan": 3}, {}],
           [None, {}],[None, {}]],subplot_titles =("<span style='font-size: 36px;''>TPC Noise Vs. Wire Length</span>","<span style='font-size: 20px;''>Longest Wires<br>RMS per Plane</span>"),column_widths=[0.8, 0.2],horizontal_spacing = 0.055,vertical_spacing = 0.085)

fig.add_trace(go.Scatter(x=Noise_df['r'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1],y = Noise_df['Raw_rms'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1],marker_color = 'DarkGreen',mode='markers',opacity=.8,text=Noise_df['Channel_id'][Noise_df['plane']==0][Noise_df['Raw_rms']>1.1],name='1<sup>st</sup> Induction'),row = 1, col = 1)
fig.add_trace(go.Scatter(x=Noise_df['r'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1],y = Noise_df['Raw_rms'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1],marker_color = 'darkorange',mode='markers',opacity=.8,text=Noise_df['Channel_id'][Noise_df['plane']==1][Noise_df['Raw_rms']>1.1],name='2<sup>nd</sup> Induction'),row = 1, col = 1)
fig.add_trace(go.Scatter(x=Noise_df['r'][Noise_df['plane']==2][Noise_df['Raw_rms']>1.4],y = Noise_df['Raw_rms'][Noise_df['plane']==2][Noise_df['Raw_rms']>1.4],marker_color = 'Purple',mode='markers',opacity=.8,text=Noise_df['Channel_id'][Noise_df['plane']==2][Noise_df['Raw_rms']>1.4],name='Collection'),row = 1, col = 1)


fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_U],marker_color = 'darkgreen',xbins=dict(start = mean_U - 5,end = mean_U+5,size=.05),showlegend = False),row = 1, col =2)


fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_V],marker_color = 'darkorange',xbins=dict(start = mean_V - 5,end = mean_V+5,size=.05),showlegend = False),row = 2, col =2)


fig.add_trace(go.Histogram(x=Noise_df['Raw_rms'][mask_Y],marker_color = 'purple',xbins=dict(start = mean_Y - 5,end = mean_Y + 5,size=.05),showlegend = False),row = 3, col =2)



fig.add_vrect(
    x0="287", x1="599",
     opacity=0.7,
    layer="below", line_width=2,row = 1,col =1
)

#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))
#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=120,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=110,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)

fig.update_layout(yaxis2 = dict(title="Channels/0.05 ADC",side="right",title_font=dict(size=14),title_standoff = 25),
                  yaxis3 = dict(title="Channels/0.05 ADC",side="right",title_font=dict(size=14),title_standoff = 25),
                  yaxis4 = dict(title="Channels/0.05 ADC",side="right",title_font=dict(size=14)),)
fig.update_layout(xaxis2 = dict(title="RMS [ADC]",side="left",title_font=dict(size=14),title_standoff= 5),
                  xaxis3 = dict(title="RMS [ADC]",side="left",title_font=dict(size=14),title_standoff= 5),
                  xaxis4 = dict(title="RMS [ADC]",side="left",title_font=dict(size=14),title_standoff= 5),)
fig.update_layout(xaxis2 = dict(range = [0,4.2]),xaxis4 = dict(range = [0,4.2]),xaxis3 = dict(range = [0,4.2]),)

fig.update_yaxes(title_text = "RMS [ADC]",title_font=dict(size=24),row = 1, col = 1)

fig.update_xaxes(title_text = "Wire length [cm]",title_font=dict(size=24),row = 1, col = 1)

#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
fig.update_layout(yaxis = dict(range = [1,3.5]))
fig.update_layout(xaxis = dict(range= [-10,600],tickmode = 'linear',dtick = 100))
fig.update_annotations(font_size=20)
fig.update_layout(font = dict(size=20),legend=dict(
    yanchor="top",
    y=0.99,
    xanchor="left",
    x=0.01,
    bgcolor="LightSteelBlue",
        bordercolor="Black",
))
#fig.add_annotation(dict(font = dict(size = 25,color="Black",)),xshift= -270,yshift=150,text = f"SBND<br>Preliminary Data",showarrow = False)
fig.add_annotation(dict(font = dict(size = 25,color="Black",)),opacity=0.7,xshift= 160,yshift=-270,text = f"Induction wires connected<br>across two wire frames",showarrow = False)

fig.update_layout(height = 800, width = 1200,showlegend = True)

fig.write_image(f'/Users/danielcarber/Documents/SBND/Noise_Analysis/Plots/RMS_vs_wire_length.png')
fig.write_image(f'/Users/danielcarber/Documents/SBND/Noise_Analysis/Plots/RMS_vs_wire_length.pdf')

fig.show()