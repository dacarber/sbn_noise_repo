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
anno = input("Do you want annotations? ")
noise = input("What Noise do you want (raw, coh, or int): ")
directory = f"/Users/danielcarber/Documents/SBND/Noise_Analysis/Plots/run{Run_num}/"
if not os.path.exists(directory):
    os.mkdir(directory)
files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise_Analysis/data/noise_output_run{Run_num}.root")
files['tpc_noise;1'].keys()

raw_rms = files['tpc_noise;1'][f'{noise}_rms'].array().to_list()
Noise_df = {'Channel_id':[],'Raw_rms':[],'wire_plane':[]}
for r,rms in enumerate(raw_rms):
    Noise_df['Channel_id'].append(r)
    if rms ==0:
        Noise_df['Raw_rms'].append(None)
    else:
        Noise_df['Raw_rms'].append(rms)
    if r <1984:
        Noise_df['wire_plane'].append('UB')
    elif r<3968:
        Noise_df['wire_plane'].append('VB')
    elif r<5638:#5632
        Noise_df['wire_plane'].append('YB')
    elif r<7616:
        Noise_df['wire_plane'].append('UA')
    elif r<9606:#9600
        Noise_df['wire_plane'].append('VA')
    elif r<11270:#11264
        Noise_df['wire_plane'].append('YA')
Noise_df = pd.DataFrame(Noise_df)

filename = f"RMS_channels_{Run_num}"
Non_zero_mask = Noise_df['Raw_rms'] > 0
print("Average noise of SBND: ",np.mean(Noise_df['Raw_rms'][Non_zero_mask]),"Median noise of SBND: ",np.median(Noise_df['Raw_rms'][Non_zero_mask]))
fig = make_subplots(rows=3,cols=2,column_widths = [0.5,0.5],subplot_titles = (r'West 1<sup>st</sup> Induction',f'East 1<sup>st</sup> Induction',f'West 2<sup>nd</sup> Induction',f'East 2<sup>nd</sup> Induction',f'West Collection',f'East Collection',),shared_yaxes=False,horizontal_spacing = 0.01)
mask = Noise_df['wire_plane'] == 'UB'
#median = np.median(Noise_df['Raw_rms'][mask])
#mean = np.mean(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Scatter(x=list(range(0,1984)),y = Noise_df['Raw_rms'][mask],marker_color = 'darkgreen'),row = 1, col = 2)
mask = Noise_df['wire_plane'] == 'UA'
#Noise_df['Channel_id'][mask] for LArSoft Channels, list(range(0,1984)) for local channel number
fig.add_trace(go.Scatter(x=list(range(0,1984)),y = Noise_df['Raw_rms'][mask],marker_color = 'darkgreen'),row = 1, col = 1)
#Noise_df['Channel_id'][mask] for LArSoft Channels, list(range(0,1984)) for local channel number

mask = Noise_df['wire_plane'] == 'VB'
fig.add_trace(go.Scatter(x=list(range(0,1984)),y = Noise_df['Raw_rms'][mask],marker_color = 'darkorange'),row = 2, col = 2)
#Noise_df['Channel_id'][mask] for LArSoft Channels, list(range(0,1984)) for local channel number
mask = Noise_df['wire_plane'] == 'VA'
fig.add_trace(go.Scatter(x=list(range(0,1984)),y = Noise_df['Raw_rms'][mask],marker_color = 'darkorange'),row = 2, col = 1)

#Noise_df['Channel_id'][mask] for LArSoft Channels, list(range(0,1984)) for local channel number
mask = Noise_df['wire_plane'] == 'YB'
fig.add_trace(go.Scatter(x=list(range(0,1664)),y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 3, col = 2)
#Noise_df['Channel_id'][mask] for LArSoft Channels, list(range(0,1664)) for local channel number
mask = Noise_df['wire_plane'] == 'YA'
fig.add_trace(go.Scatter(x=list(range(0,1664)),y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 3, col = 1)
fig.add_vrect(
    x0="1248", x1="1280",
    fillcolor="Grey", opacity=0.5,
    layer="below", line_width=0,row=2,col = 2
)
fig.add_vrect(
    x0="191", x1="224",
    fillcolor="Grey", opacity=0.5,
    layer="below", line_width=0,row=3,col = 2
)

#fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 1)
#fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 1)
#fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 1)
#fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 2)
#fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 2)
#fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 2)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 1, col = 1)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 2, col = 1)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 3, col = 1)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 1, col = 2)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 2, col = 2)
fig.update_xaxes(title_text = "TPC Plane Channel Number",row = 3, col = 2)
if Run_num != 'sim':
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 1, col =1)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 1, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 2, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 2, col =1)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 3, col =2)
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),xshift= 730,yshift=40,text = f"<b>SBND<br>Preliminary Data</b>",showarrow = False),row = 3, col =1)
if anno == "Yes" or anno == "yes":
    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 1528,y=3,text = f"Shorted Wire",showarrow = True,arrowhead=1,xanchor="right",arrowwidth=2,arrowcolor="Black"),row = 1, col =1)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 751,y=3,text = f"Shorted Wire",showarrow = True,arrowhead=1,xanchor="right",arrowwidth=2,arrowcolor="Black"),row = 2, col =1)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 952,y=1.5,text = f"Disconnected Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 2, col =1)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 406,y=1.25,text = f"No Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 3, col =1)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 1257,y=1.25,text = f"No Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 3, col =1)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 406,y=1.25,text = f"No Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 3, col =2)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 1257,y=1.25,text = f"No Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 3, col =2)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 200,y=3,text = f"Non-responisve channels",showarrow = True,arrowhead=1,xanchor="left",ax=30,arrowwidth=2,arrowcolor="Black"),row = 3, col =2)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 1260,y=4,text = f"Non-responisve channels",showarrow = True,arrowhead=1,xanchor="right",ay=0,ax=-20,arrowwidth=2,arrowcolor="Black"),row = 2, col =2)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 546,y=1.5,text = f"Disconnected Wire",showarrow = True,arrowhead=1,xanchor="right",yanchor="top",ay=10,arrowwidth=2,arrowcolor="Black"),row = 1, col =2)

    fig.add_annotation(dict(font = dict(size = 15,color="Black",),x= 606,y=3.2,text = f"Interconnected Wires",showarrow = True,arrowhead=1,xanchor="right",ay=-10,arrowwidth=2,arrowcolor="Black"),row = 1, col =2)




#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
fig.update_layout(yaxis = dict(range = [0,5],title="RMS [ADC]",side="left",),
                  yaxis2 = dict(range = [0,5],title="RMS [ADC]",side="right",),
                  yaxis3 = dict(range = [0,5],title="RMS [ADC]",side="left",),
                  yaxis4 = dict(range = [0,5],title="RMS [ADC]",side="right",),
                  yaxis5 = dict(range = [0,5],title="RMS [ADC]",side="left",),
                  yaxis6 = dict(range = [0,5],title="RMS [ADC]",side="right",))
tick_size = 128 #128 64
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = tick_size),
                  xaxis2 = dict(tickmode = 'linear',dtick = tick_size),
                  xaxis3 = dict(tickmode = 'linear',dtick = tick_size),
                  xaxis4 = dict(tickmode = 'linear',dtick = tick_size),
                  xaxis5 = dict(tickmode = 'linear',dtick = tick_size),
                  xaxis6 = dict(tickmode = 'linear',dtick = tick_size))
fig.update_layout(height = 800, width = 1800,showlegend = False)
fig.update_layout(
    title={
        'text': "<span style='font-size: 26px;'>TPC Noise per Channel</span>",
        'y':.95,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top',
        })
print(directory+filename)
fig.write_image(directory+filename+".png")
fig.write_image(directory+filename+".pdf")
time.sleep(0.1)
fig.write_image(directory+filename+".pdf")
fig.show()


