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


Run_num = input("Enter the Run Number: ")
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
    elif r<11264:
        Noise_df['wire_plane'].append('YA')
Noise_df = pd.DataFrame(Noise_df)

filename = f"RMS_channels_{Run_num}.png"
Non_zero_mask = Noise_df['Raw_rms'] > 0
print("Average noise of SBND: ",np.mean(Noise_df['Raw_rms'][Non_zero_mask]),"Median noise of SBND: ",np.median(Noise_df['Raw_rms'][Non_zero_mask]))
fig = make_subplots(rows=3,cols=2,column_widths = [0.5,0.5],subplot_titles = (r'West 1st Induction',f'East 1st Induction',f'West 2nd Induction',f'East 2nd Induction',f'West Collection',f'East Collection',),shared_yaxes=True,
                    horizontal_spacing=0.01)
mask = Noise_df['wire_plane'] == 'UB'
#median = np.median(Noise_df['Raw_rms'][mask])
#mean = np.mean(Noise_df['Raw_rms'][mask])
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'red'),row = 1, col = 2)
mask = Noise_df['wire_plane'] == 'UA'
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'red'),row = 1, col = 1)
#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))

mask = Noise_df['wire_plane'] == 'VB'
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 2, col = 2)

mask = Noise_df['wire_plane'] == 'VA'
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'purple'),row = 2, col = 1)


mask = Noise_df['wire_plane'] == 'YB'
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'blue'),row = 3, col = 2)

mask = Noise_df['wire_plane'] == 'YA'
fig.add_trace(go.Scatter(x=Noise_df['Channel_id'][mask],y = Noise_df['Raw_rms'][mask],marker_color = 'blue'),row = 3, col = 1)


fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 1, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 2, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 3, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 1, col = 2)
fig.update_xaxes(title_text = "Channel #",row = 2, col = 2)
fig.update_xaxes(title_text = "Channel #",row = 3, col = 2)
#fig.update_layout(xaxis2 = dict(range = [0,5]),xaxis4 = dict(range = [0,5]),xaxis6 = dict(range = [0,5]))
fig.update_layout(yaxis = dict(range = [0,5]),yaxis2 = dict(range = [0,5]),yaxis3 = dict(range = [0,5]),yaxis4 = dict(range = [0,5]),yaxis5 = dict(range = [0,5]),yaxis6 = dict(range = [0,5]))
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = 64),xaxis2 = dict(tickmode = 'linear',dtick = 64),xaxis3 = dict(tickmode = 'linear',dtick = 64),xaxis4 = dict(tickmode = 'linear',dtick = 64),xaxis5 = dict(tickmode = 'linear',dtick = 64),xaxis6 = dict(tickmode = 'linear',dtick = 64))
fig.update_layout(height = 800, width = 1800,showlegend = False,title_text=f"Channel RMS Run {Run_num}",title_x=0.5)

fig.write_image(directory+filename)
fig.show()


