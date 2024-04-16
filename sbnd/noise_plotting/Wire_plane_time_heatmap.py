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
    def __init__(self, waveform_df, calc):
        self.waveform_df = waveform_df
        self.calc = calc
        self.range = (0,200)
        self.type = 'Tick'
        
    def Run_calc(self):
        calculation = []
        bin_width = 1
        if self.calc == 'RMS':
            calculation = self.RMS_calc()
            self.range = (0,1000)
        elif self.calc == 'Max':
            calculation = self.Max_calc()
        elif self.calc == 'Min':
            calculation = self.Min_calc()
        elif self.calc == 'Range':
            calculation = self.Range_calc()
            self.range = (0,30)
        elif self.calc == 'Rise':
            calculation = self.Rise_calc()
            self.range = (0,100)
        elif self.calc == 'Slope':
            calculation = self.Slope_calc()
            self.range = (0,30)
        elif self.calc == 'Zero':
            calculation = self.Zero_calc()
            self.range = (0,20)
        elif self.calc == 'Start':
            calculation = self.Start_calc()
            self.range = (calculation[0]-10,calculation[0]+10)
            bin_width = .5
        elif self.calc == 'Lead_rise':
            calculation = self.Lead_rise_calc()
            self.range = (0,100)
            bin_width = .5
        elif self.calc == 'Ledge':
            calculation = self.Ledge_calc()
            self.range = (0,300)
            bin_width = 1
        else:
            print("You didn't enter a correct calculation\n Please enter Max, Min, Rise, or Range")
        square = 0
        for tick in calculation:
            square += (tick - np.mean(calculation))*(tick - np.mean(calculation))
            mean = square / len(calculation)
            RMSE =np.sqrt(mean)
        #print(calculation)
        nbins = int((max(calculation)-min(calculation))/bin_width)
        print(nbins)
        fig = px.histogram(x=calculation,nbins=nbins)
        fig.update_layout(xaxis_title=self.calc)
        fig.update_layout(height = 600, width = 1000,showlegend = False)
        fig.add_annotation(dict(font = dict(size = 20),xshift=350,yshift=400,text = f"RMSE:{RMSE:.2f}",showarrow = False))
        fig.update_layout(xaxis_range=self.range)
        fig.show()
        return calculation
    '''  
    def RMS_calc(self):
        RMS = []
        for keys in tqdm(self.waveform_df):
            waveform  = self.waveform_df[keys]
            square = 0
            for tick in range(waveform):
                square += (tick - np.mean(waveform))*(tick - np.mean(waveform))
            mean = square / len(waveform)
            RMS.append(np.sqrt(mean))
        return RMS
    '''
    def Max_calc(self):
        waveform_max = []
        for keys in self.waveform_df:
            waveform_max.append(np.argmax(self.waveform_df[keys]))
        return waveform_max
    def Min_calc(self):
        waveform_min = []
        for keys in self.waveform_df:
            waveform_min.append(np.argmin(self.waveform_df[keys]))
        return waveform_min
    def Range_calc(self):
        waveform_time_range = []
        for keys in self.waveform_df:
            number_0 = len(self.waveform_df[keys][self.waveform_df[keys] == 0])
            number_max = len(self.waveform_df[keys][self.waveform_df[keys] == np.max(self.waveform_df[keys])])
            waveform_time_range.append(abs(np.argmax(self.waveform_df[keys])+int(number_max/2)-(np.argmin(self.waveform_df[keys])+(int(number_0/2)))))
        return waveform_time_range
    def Rise_calc(self):
        rise_time = []
        for keys in self.waveform_df:
            waveform = self.waveform_df[keys]
            median = np.median(waveform[:np.argmax(waveform)-50])
            i = np.argmax(waveform)
            if i == 0:
                rise_time.append(-1)
                continue
            while waveform[i]>=median:
                i-=1
                if i == 0:
                    rise_time.append(-1)
                    continue
                #print(i)
            rise_time.append(np.argmax(waveform)-i)
        return rise_time

    def Slope_calc(self):
        waveform_time_slope = []
        for keys in self.waveform_df:
            pulse = self.waveform_df[keys][np.argmax(self.waveform_df[keys]):np.argmax(self.waveform_df[keys])+50]
            x=0
            for tick in pulse:
                x+=1
                if tick == 0:
                    break
            y = np.min(self.waveform_df[keys]) - np.max(self.waveform_df[keys])
            waveform_time_slope.append(y/x)
        return waveform_time_slope
    def Zero_calc(self):
        waveform_time_zero = []
        for keys in self.waveform_df:
            number_0 = len(self.waveform_df[keys][self.waveform_df[keys] ==np.min(self.waveform_df[keys])])
            waveform_time_zero.append(number_0)
        return waveform_time_zero
    def Start_calc(self):
        waveform_time_start = []
        for keys in self.waveform_df:
            waveform = self.waveform_df[keys]
            pulse_start = 0
            for tick in range(1,len(waveform),1):
                    if abs(waveform[tick]-waveform[tick-1]) > 20:
                        pulse_start = tick-1
                        break
            waveform_time_start.append(pulse_start)
        return waveform_time_start
    def Lead_rise_calc(self):
        waveform_time_lead = []
        k=0
        waveform = [0]*3415
        for keys in self.waveform_df:
            pulse = self.waveform_df[keys][:np.argmax(self.waveform_df[keys])-100]
            number_min = len(pulse[pulse < 500])
            number_max = len(pulse[pulse > 2500])
            if (number_min > 1 or number_max>1):
                continue
            waveform += self.waveform_df[keys]
            k+=1
        waveform = waveform/k
        rise_start = 0
        for tick in range(0,len(waveform),1):
                if abs(waveform[tick+20]-waveform[tick]) >=2:
                    rise_start = tick
                    break
        for keys in self.waveform_df:
            waveform = self.waveform_df[keys]
            pulse_start = 0
            for tick in range(1,len(waveform),1):
                    if abs(waveform[tick]-waveform[tick-1]) > 20:
                        pulse_start = tick-1
                        break
            waveform_time_lead.append(pulse_start-rise_start)
        return waveform_time_lead
    def Ledge_calc(self):
        waveform_time_zero = []
        long_ledge = 0
        ledge_time_end = 0
        ledge_time_beg = 0
        ledge_time = []
        for keys in self.waveform_df:
            number_0 = len(self.waveform_df[keys][(self.waveform_df[keys] ==np.min(self.waveform_df[keys]))&(np.min(self.waveform_df[keys]) < 100)])
            waveform = self.waveform_df[keys]
            waveform = list(waveform - np.mean(waveform[0:50]))
            beg_ped = np.mean(waveform[0:50])
            waveform.reverse()
            end_ped = np.mean(waveform[0:50])
            if abs(end_ped -beg_ped) >= 150:
                print("Waveform ledge goes all the way to the end")
                ledge_time.append(np.argmax(waveform))
                long_ledge +=1
                continue
            for t,tick in enumerate(waveform):
                if waveform[t] > 100 and waveform[t+1]-waveform[t] < 10 :
                    ledge_time_beg = t
                    break
            if abs(np.argmax(waveform)-ledge_time_beg) < 60:
                #print("No ledge")
                ledge_time.append(-10)
                continue
            if number_0 >1:
                #print("No ledge")
                ledge_time.append(-10)
                continue
            for k in range(ledge_time_beg, len(waveform),1):
                if waveform[k+1]-waveform[k] < -10:
                    ledge_time_end=0
                    ledge_time_bed = 10
                    #print("No ledge")
                    break
                if waveform[k+1]-waveform[k] > 10:
                    ledge_time_end = k
                    break
                if k == 3413:
                    break
            #print(ledge_time_end,ledge_time_beg)
            ledge_time.append(ledge_time_end-ledge_time_beg)
        print("Number of channels with ledges past event display: ",long_ledge)
        return ledge_time
        
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

def plot_wireplanes(event_number,metric):
    if not os.path.exists(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/"):
        os.mkdir(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/")
    directory = f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/{event_number}/time_{metric}/"
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()

    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
    
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"East TPC First Ind Wire {metric} Signal",color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(height = 800, width = 1200,showlegend = False)
    
    fig.write_image(directory+f'UB_plane_diagram_time_{event_number}_{metric}.png')
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()

    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"East TPC Second Ind Wire {str(metric)} Signal",color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f"VB_plane_diagram_time_{event_number}_{metric}.png")
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()

    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}

    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"East TPC Coll Wire {str(metric)} Signal",color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'YB_plane_diagram_time_{event_number}_{metric}.png')
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()

    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"West TPC First Ind Wire {str(metric)} Signal",color_continuous_scale=px.colors.sequential.Viridis)


    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'UA_plane_diagram_time_{event_number}_{metric}.png')
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    
    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"West TPC Second Ind Wire {str(metric)} Signal",color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'/VA_plane_diagram_time_{event_number}_{metric}.png')
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
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    
    print("Max: ",max(waveform)*.5,"Min: ",min(waveform)*.5,"Mean: ",np.mean(waveform)*.5)
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask])*.5,"Min within range: ",np.min(waveform[mask])*.5,"Mean within range: ",np.mean(waveform[mask])*.5)
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
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),"Tick":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = "Tick",range_color=calc.range,title=f"West TPC Coll Wire {str(metric)} Signal",color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'/YA_plane_diagram_time_{event_number}_{metric}.png')
    print("Done")
    fig.show()
                    

def main():

    event_number = input("Please enter event number:")
    metric = input("Please enter calculation type of time (Max, Min, Rise, Slope, Zero, Start, Lead_rise,Ledge or Range):")
    plot_wireplanes(event_number,str(metric))
    
if __name__=="__main__": 
    main() 