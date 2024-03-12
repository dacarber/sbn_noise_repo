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
from plotly.subplots import make_subplots
from plotly import tools
import plotly.offline as pyo
import sys

file = uproot.open("/Users/danielcarber/Documents/ICARUS/Noise_Analysis/TPC_Noise_repo/Waveforms/run11524_50.root")
#file2 = uproot.open("/Users/danielcarber/Documents/ICARUS/Noise_Analysis/TPC_Noise_repo/Waveforms/run11676.root")

run_num = 11524
#run_num_2 = 11676
df = {}
for key in file["tpcnoisefull;1"]['tpcnoise;5'].keys():
    print(key)
    df[key] = file["tpcnoisefull;1"]['tpcnoise;5'][key].array().to_list()
    
#df_2 = {}
#for key in file2["tpcnoisefull;1"]['tpcnoise;9'].keys():
#    print(key)
#    df_2[key] = file2["tpcnoisefull;1"]['tpcnoise;9'][key].array().to_list()


df=pd.DataFrame.from_dict(df)
#df_2=pd.DataFrame.from_dict(df_2)
channel_map = pd.read_csv('/Users/danielcarber/Documents/ICARUS/Noise_Analysis/TPC_Noise_repo/icarus_channels.csv')
df = df.merge(channel_map, on='channel_id')
#df_2 = df_2.merge(channel_map, on='channel_id')


#plt.scatter(list(raw_rms_2.index.values),raw_rms_2['raw_rms'],alpha = 1, label = "RMS Run 11615",s = 10)
#plt.scatter(list(raw_rms.index.values),raw_rms['raw_rms'],alpha = .6, label = "RMS Run 11524",s = 10)
filename = f"RMS_plots_EE.pdf"
directory = "/Users/danielcarber/Documents/ICARUS/Noise_Analysis/"
fig = make_subplots(rows=3,cols=2,column_widths = [0.7,0.3],subplot_titles = (f'EE Induction 1','EE Ind 1 RMS',f'EE Induction 2','EE Ind 2 RMS',f'EE Collection','EE Coll RMS',))

mask = (df['plane_number'] ==0) & (df['tpc_number'] ==0)
#mask2 = (df_2['plane_number'] ==0) & (df_2['tpc_number'] ==0)
raw_rms = df[mask].groupby(['channel_id'])['coh_rms'].median().to_frame()
#raw_rms_2 = df_2[mask2].groupby(['channel_id'])['coh_rms'].median().to_frame()
median = np.median(raw_rms)
mean = np.mean(raw_rms)
#median2 = np.median(raw_rms_2)
#mean2 = np.mean(raw_rms_2)
#fig.add_trace(go.Histogram(x=raw_rms_2['coh_rms'],marker_color = 'red',xbins=dict(start = median - 10,end = median+10,size=.2),showlegend = False),row = 1, col =2)

#fig.add_trace(go.Scatter(x=list(raw_rms_2.index.values),y = raw_rms_2['coh_rms'],marker_color = 'red',mode='markers',name = f'Run {run_num_2}',opacity = .5),row = 1, col = 1)

fig.add_trace(go.Histogram(x=raw_rms['coh_rms'],marker_color = 'blue',xbins=dict(start = median - 10,end = median+10,size=.2),showlegend = False),row = 1, col =2)

fig.add_trace(go.Scatter(x=list(raw_rms.index.values),y = raw_rms['coh_rms'],marker_color = 'blue',mode='markers',name = f'Run {run_num}',opacity = .5),row = 1, col = 1)

#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=120,text = f"Run {run_num}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=106,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =1,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=92,text = f"Median RMS:{median:.2f}",showarrow = False),row =1,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=60,text = f"Run {run_num_2}",showarrow = False),row =1,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=46,text = f"Mean RMS:{mean2:.2f}",showarrow = False),row =1,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=32,text = f"Median RMS:{median2:.2f}",showarrow = False),row =1,col=2)

#Plane 1 EE TPC
mask = (df['plane_number'] ==1) & (df['tpc_number'] ==0)
#mask2 = (df_2['plane_number'] ==1) & (df_2['tpc_number'] ==0)
raw_rms = df[mask].groupby(['channel_id'])['coh_rms'].median().to_frame()
#raw_rms_2 = df_2[mask2].groupby(['channel_id'])['coh_rms'].median().to_frame()
median = np.median(raw_rms)
mean = np.mean(raw_rms)
#median2 = np.median(raw_rms_2)
#mean2 = np.mean(raw_rms_2)
#fig.add_trace(go.Histogram(x=raw_rms_2['coh_rms'],marker_color = 'red',xbins=dict(start = median - 10,end = median+10,size=.05),showlegend = False),row = 2, col =2)

#fig.add_trace(go.Scatter(x=list(raw_rms_2.index.values),y = raw_rms_2['coh_rms'],marker_color = 'red',mode='markers',name = f'Run {run_num_2}',opacity = .5,showlegend = False),row = 2, col = 1)

fig.add_trace(go.Histogram(x=raw_rms['coh_rms'],marker_color = 'blue',xbins=dict(start = median - 10,end = median+10,size=.05),showlegend = False),row = 2, col =2)

fig.add_trace(go.Scatter(x=list(raw_rms.index.values),y = raw_rms['coh_rms'],marker_color = 'blue',mode='markers',name = f'Run {run_num}',opacity = .5,showlegend = False),row = 2, col = 1)

#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=100,text = f"Run {run_num}",showarrow = False),row =2,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=86,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =2,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=72,text = f"Median RMS:{median:.2f}",showarrow = False),row =2,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=40,text = f"Run {run_num_2}",showarrow = False),row =2,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=26,text = f"Mean RMS:{mean2:.2f}",showarrow = False),row =2,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=12,text = f"Median RMS:{median2:.2f}",showarrow = False),row =2,col=2)

#Plane 2 EE TPC
mask = (df['plane_number'] ==2) & (df['tpc_number'] ==0)
#mask2 = (df_2['plane_number'] ==2) & (df_2['tpc_number'] ==0)
raw_rms = df[mask].groupby(['channel_id'])['coh_rms'].median().to_frame()
#raw_rms_2 = df_2[mask2].groupby(['channel_id'])['coh_rms'].median().to_frame()
median = np.median(raw_rms)
mean = np.mean(raw_rms)
#median2 = np.median(raw_rms_2)
#mean2 = np.mean(raw_rms_2)
#fig.add_trace(go.Histogram(x=raw_rms_2['coh_rms'],marker_color = 'red',xbins=dict(start = median - 10,end = median+10,size=.05),showlegend = False),row = 3, col =2)

#fig.add_trace(go.Scatter(x=list(raw_rms_2.index.values),y = raw_rms_2['coh_rms'],marker_color = 'red',mode='markers',name = f'Run {run_num_2}',opacity = .5,showlegend = False),row = 3, col = 1)

fig.add_trace(go.Histogram(x=raw_rms['coh_rms'],marker_color = 'blue',xbins=dict(start = median - 10,end = median+10,size=.05),showlegend = False),row = 3, col =2)

fig.add_trace(go.Scatter(x=list(raw_rms.index.values),y = raw_rms['coh_rms'],marker_color = 'blue',mode='markers',name = f'Run {run_num}',opacity = .5,showlegend = False),row = 3, col = 1)

#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=100,text = f"Run {run_num}",showarrow = False),row =3,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=86,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =3,col=2)
fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=72,text = f"Median RMS:{median:.2f}",showarrow = False),row =3,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=40,text = f"Run {run_num_2}",showarrow = False),row =3,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=26,text = f"Mean RMS:{mean2:.2f}",showarrow = False),row =3,col=2)
#fig.add_annotation(dict(font = dict(size = 14),xshift= 80,yshift=12,text = f"Median RMS:{median2:.2f}",showarrow = False),row =3,col=2)
#mask = (df['channel_id'] >=13824) & (df['channel_id'] < 27648)
#median = np.median(df['raw_rms'][mask])
#mean = np.mean(df['raw_rms'][mask])
#fig.add_trace(go.Histogram(x=df['raw_rms'],marker_color = 'red',xbins=dict(start = median - 5,end = median+5,size=.05)),row = 2, col =2)
#fig.add_trace(go.Scatter(x=df['channel_id'],y = df['raw_rms'],marker_color = 'red'),row = 2, col = 1)
#fig.update_layout(xaxis2 = dict(range = [0,median+5]))
#fig.update_layout(margin = dict(r=200))
#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=60,text = f"Mean RMS:{mean:.2f}",showarrow = False),row =2,col=2)
#fig.add_annotation(dict(font = dict(size = 10),xshift= 180,yshift=50,text = f"Median RMS:{median:.2f}",showarrow = False),row =2,col=2)


fig.update_yaxes(title_text = "RMS [ADC]",row = 1, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 2, col = 1)
fig.update_yaxes(title_text = "RMS [ADC]",row = 3, col = 1)

fig.update_xaxes(title_text = "Channel #",row = 1, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 2, col = 1)
fig.update_xaxes(title_text = "Channel #",row = 3, col = 1)

fig.update_xaxes(title_text = "RMS [ADC]",row = 1, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 2, col = 2)
fig.update_xaxes(title_text = "RMS [ADC]",row = 3, col = 2)

fig.update_layout(xaxis2 = dict(range = [0,10]),xaxis4 = dict(range = [0,7]),xaxis6 = dict(range = [0,7]))
fig.update_layout(yaxis = dict(range = [0,7]),yaxis3 = dict(range = [0,7]),yaxis5 = dict(range = [0,7]))
fig.update_layout(xaxis = dict(tickmode = 'linear',dtick = 288),xaxis3 = dict(tickmode = 'linear',dtick = 576),xaxis5 = dict(tickmode = 'linear',dtick = 576))
fig.update_layout(height = 800, width = 1500,showlegend = True)

fig.write_image(f"/Users/danielcarber/Documents/ICARUS/Noise_Analysis/Coh_RMS_EE_{run_num}.pdf")
fig.show()

