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
from tqdm import tqdm
class waveform_calc:
    def __init__(self, waveform_df, calc,grouping,wire_plane,APA):
        self.waveform_df = waveform_df
        self.calc = calc
        self.range = (0,5000)
        self.type = 'ADC'
        self.grouping = grouping
        self.wire_plane = wire_plane
        self.APA = APA
    def Run_calc(self):
        calculation = []
        bin_width = 10
        if self.calc == 'RMS':
            calculation = self.RMS_calc()
            self.range = (0,10)
            bin_width = .001
        elif self.calc == 'Max':
            calculation = self.Max_calc()
        elif self.calc == 'Min':
            calculation = self.Min_calc()
        elif self.calc == 'Range':
            calculation = self.Range_calc()
        elif self.calc == 'Max_min_time':
            self.range = (0,30)
            calculation = self.Max_min_time_calc()
            self.type = 'Tick'
        elif self.calc == 'Integral':
            calculation = self.Integral_calc()
            self.range = (-50000,50000)
            #bin_width = .01
        else:
            print("You didn't enter a correct calculation\n Please enter Max, Min, Range, or RMS")
        square = 0
        for tick in calculation:
            square += (tick - np.mean(calculation))*(tick - np.mean(calculation))
            mean = square / len(calculation)
            RMSE =np.sqrt(mean)
        
        nbins = int((max(calculation)-min(calculation))/bin_width)
        #print(nbins)
        fig = px.histogram(x=calculation,nbins=nbins)
        fig.update_layout(xaxis_title=self.calc)
        fig.update_layout(height = 600, width = 1000,showlegend = False)
        fig.add_annotation(dict(font = dict(size = 20),xshift=350,yshift=400,text = f"RMSE:{RMSE:.2f}",showarrow = False))
        fig.update_layout(xaxis_range=self.range)
        fig.show()
        calculation = self.group_sel(wires = calculation)
        return calculation
      
    def RMS_calc(self):
        RMS = []
        for keys in tqdm(self.waveform_df):
            waveform  = self.waveform_df[keys]-np.median(self.waveform_df[keys])
            square = 0
            for tick in waveform:
                #square += (tick - np.mean(waveform))*(tick - np.mean(waveform))
                square += (tick)*(tick)
            mean = square / len(waveform)
            RMS.append(np.sqrt(mean))
            #if np.sqrt(mean) < 1.35:
                #print(keys)
        return RMS
    def Max_calc(self):
        waveform_max = []
        for keys in self.waveform_df:
            waveform_max.append(np.max(self.waveform_df[keys]))
        return waveform_max
    def Min_calc(self):
        waveform_min = []
        for keys in self.waveform_df:
            waveform_min.append(np.min(self.waveform_df[keys]))
        return waveform_min
    def Range_calc(self):
        waveform_range = []
        for keys in self.waveform_df:
            waveform_range.append(np.max(self.waveform_df[keys])-np.min(self.waveform_df[keys]))
        return waveform_range
    def Max_min_time_calc(self):
        waveform_time_range = []
        for keys in self.waveform_df:
            waveform_time_range.append(abs(np.argmax(self.waveform_df[keys])-np.argmin(self.waveform_df[keys])))
        return waveform_time_range
    def Integral_calc(self):
        waveform_sum = []
        for keys in self.waveform_df:
            integral = 0
            #median = 0
            median = np.median(self.waveform_df[keys][:np.argmax(self.waveform_df[keys])-50])
            pulse = self.waveform_df[keys][np.argmax(self.waveform_df[keys])-200:np.argmax(self.waveform_df[keys])+200]-median 
            
            for k in pulse:
                integral += k

            waveform_sum.append(integral)
        return waveform_sum
    def group_sel(self,wires):
        channel_map = pd.read_csv("../datafiles/channel_mapping.txt", sep = " ")
        channel_map = channel_map[(channel_map["Wire_plane"] == self.wire_plane) & (channel_map["APA"] == self.APA)]
        groups = channel_map[self.grouping].unique()
        channel_map = channel_map.sort_values('LArSoft_ch')
        wires = np.array(wires)
        if self.grouping == None:
            return wires
        for group in groups:
            mask = list(channel_map[self.grouping] == group)
            
            #print(wires[mask])
            new_wires = wires
            mean = np.mean(wires[mask])
            for w,wire in enumerate(wires):
                if mask[w]== True:
                    new_wires[w] = mean
                else:
                    continue
            #new_wires.extend([np.mean(wires[mask])]*sum(mask))
        return new_wires
def threshold_info(waveform,value,group = "LArSoft_ch", threshold = None):
    channel_map = pd.read_csv("../datafiles/channel_mapping.txt", sep = " ")
    
    if threshold is None:
        return
    elif threshold == "Greater":
        for w,wire in enumerate(waveform):
            if wire >= float(value):
                print(channel_map[channel_map["LArSoft_ch"] == w])
    elif threshold == "Less":
        for w,wire in enumerate(waveform):
            if wire <= float(value):
                print(channel_map[channel_map["LArSoft_ch"] == w])
    elif threshold == "Equal":
        for w,wire in enumerate(waveform):
            if wire == float(value):
                print(channel_map[channel_map["LArSoft_ch"] == w])
    
def load_wire_info():
    wire_plane_list = ['UB','VB','YB','UA','VA','YA']
    wire_df = {'Channel_id':[],'cryo':[],'tpc':[],'tpc':[],'plane':[],'rel_wire':[],'x_0':[],'y_0':[],'z_0':[],'x_1':[],'y_1':[],'z_1':[],'r':[]}
    wire_txt = '/Users/danielcarber/Documents/SBND/Noise Analysis/sbn_noise_repo/sbnd/datafiles/Wire_lengths.txt'

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
    #wire_df = pd.DataFrame(wire_df) 
    ch_mask = wire_df['Channel_id'] == 1
    return wire_df

def plot_wireplanes(event_number,metric,value,threshold,plot_type):
    if not os.path.exists(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/"):
        os.mkdir(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/")
    directory = f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/{metric}/"
    if not os.path.exists(directory):
        os.mkdir(directory)
    files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/waveform_output_{event_number}.root")
    print(files['tpc_noise;1'].keys())
    
    
    wire_df = load_wire_info()
    Waveform_df = {}
    ch_id = 0
    raw_rms = files['tpc_noise;1']['UB_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'U',APA = 'East')
    waveform = calc.Run_calc()
    #waveform = calc.group_selection()
    threshold_info(waveform, value = value,group = plot_type, threshold = threshold)

    print(max(waveform),min(waveform),np.mean(waveform))
    
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(0,ch_id,1)):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            c.append(waveform[ch])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "ADC",range_color=calc.range,title=f"East TPC First Ind Wire {metric} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)
    
    fig.write_image(directory+f'UB_plane_diagram_{event_number}_{metric}_{plot_type}.png')
    print("Done")
    fig.show()
    
    
# VB Plane
    
    Waveform_df = {}
    ch_id = 1984
    raw_rms = files['tpc_noise;1']['VB_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)
    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'V',APA = 'East')
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value,group =plot_type, threshold = threshold)
    
    print(max(waveform),min(waveform),np.mean(waveform))
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(1984,ch_id,1)):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-1984])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color ="ADC",range_color=calc.range,title=f"East TPC Second Ind Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f"VB_plane_diagram_{event_number}_{metric}_{plot_type}.png")
    print("Done")
    fig.show()
    

#.      YB Plane    

    Waveform_df = {}
    ch_id = 3968
    raw_rms = files['tpc_noise;1']['YB_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'Y',APA = 'East')
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value,group =plot_type, threshold = threshold)

    print(max(waveform),min(waveform),np.mean(waveform))
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(3968,ch_id,1)):
        #a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        #b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['y_0'][ch]),int(wire_df['y_1'][ch]),1):
            x_i = wire_df['z_0'][ch]
            x.append(wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-3968])
            y.append(i)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}

    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "ADC",range_color=calc.range,title=f"East TPC Coll Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'YB_plane_diagram_{event_number}_{metric}_{plot_type}.png')
    print("Done")
    fig.show()
    
    
    Waveform_df = {}
    ch_id = 5632
    raw_rms = files['tpc_noise;1']['UA_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)
    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'U',APA = 'West')
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value,group =plot_type, threshold = threshold)

    print(max(waveform),min(waveform),np.mean(waveform))
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(5632,ch_id,1)):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-5632])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "ADC",range_color=calc.range,title=f"West TPC First Ind Wire {str(metric)} Signal")


    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'UA_plane_diagram_{event_number}_{metric}_{plot_type}.png')
    print("Done")
    fig.show()
    
    
    Waveform_df = {}
    ch_id = 7616
    raw_rms = files['tpc_noise;1']['VA_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'V',APA = 'West')
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value,group =plot_type, threshold = threshold)
    
    print(max(waveform),min(waveform),np.mean(waveform))
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(7616,ch_id,1)):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-7616])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "ADC",range_color=calc.range,title=f"West TPC Second Ind Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'VA_plane_diagram_{event_number}_{metric}_{plot_type}.png')
    print("Done")
    fig.show()
                    
        
    Waveform_df = {}
    ch_id = 9600
    raw_rms = files['tpc_noise;1']['YA_plane'].array().to_list()
    for r,rms in tqdm(enumerate(raw_rms)):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric,grouping = plot_type,wire_plane = 'Y',APA = 'West')
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value,group =plot_type, threshold = threshold)
    
    print(max(waveform),min(waveform),np.mean(waveform))
    print(ch_id)
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in tqdm(range(9600,ch_id,1)):
        #a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        #b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['y_0'][ch]),int(wire_df['y_1'][ch]),1):
            x_i = wire_df['z_0'][ch]
            x.append(wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-9600])
            y.append(i)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"ADC":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "ADC",range_color=calc.range,title=f"West TPC Coll Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'YA_plane_diagram_{event_number}_{metric}_{plot_type}.png')
    print("Done")
    fig.show()

def main():

    event_number = input("Please enter event number:")
    metric = input("Please enter calculation type (RMS, Max, Min,Integral or Range):")
    channel_map = pd.read_csv("../datafiles/channel_mapping.txt", sep = " ")
    threshold = input(f"If you want to print threshold enter type of threshold (Greater, Less, or Equal): ")
    
    if (threshold != "" or threshold != None):
        value = input("Please enter value for threshold: ")
    plot_type = input(f"Enter grouping you want ({channel_map.keys()}):")
    plot_wireplanes(event_number,str(metric),value,threshold,plot_type)
    
if __name__=="__main__": 
    main() 