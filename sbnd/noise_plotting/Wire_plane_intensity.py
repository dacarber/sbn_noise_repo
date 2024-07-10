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
from datetime import datetime, timezone, timedelta

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
        self.range = (0,5000)
        self.type = 'ADC'
        self.ledge = 0
    def Run_calc(self):
        calculation = []
        bin_width = 10
        if self.calc == 'RMS':
            calculation = self.RMS_calc()
            self.range = (0,2)
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
            self.range = (-100000,100000)
            #bin_width = .01
        elif self.calc == 'Rise':
            calculation = self.Rise_calc()
            self.range = (-500,500)
        elif self.calc == 'Ledge':
            ledge_info = self.Ledge_calc()
            calculation = ledge_info[0]
            self.ledge = ledge_info[1]
            self.range = (0,2000)
            bin_width = 1
            self.type = 'Tick'
        elif self.calc == 'Rise_time':
            calculation = self.Rise_time_calc()
            self.range = (0,500)
        else:
            print("You didn't enter a correct calculation\n Please enter Max, Min, Range, or RMS")
        square = 0
        for tick in calculation:
            square += (tick - np.mean(calculation))*(tick - np.mean(calculation))
            mean = square / len(calculation)
            RMSE =np.sqrt(mean)
        '''
        nbins = int((max(calculation)-min(calculation))/bin_width)
        fig = px.histogram(x=calculation,nbins=nbins)
        fig.update_layout(xaxis_title=self.calc)
        fig.update_layout(height = 600, width = 1000,showlegend = False)
        fig.add_annotation(dict(font = dict(size = 20),xshift=350,yshift=400,text = f"RMSE:{RMSE:.2f}",showarrow = False))
        fig.update_layout(xaxis_range=self.range)
        fig.show()
        
        '''
        return calculation
      
    def RMS_calc(self):
        RMS = []
        for keys in self.waveform_df:
            waveform  = self.waveform_df[keys]-np.median(self.waveform_df[keys])
            square = 0
            for tick in waveform:
                #square += (tick - np.mean(waveform))*(tick - np.mean(waveform))
                square += (tick)*(tick)
            mean = square / len(waveform)
            RMS.append(np.sqrt(mean))
            #if np.sqrt(mean) < 1.35:
            #    print(keys)
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
    def Rise_calc(self):
        waveform_adc_lead = []
        for keys in self.waveform_df:
            waveform = self.waveform_df[keys]
            median = np.median(self.waveform_df[keys][:np.argmax(self.waveform_df[keys])-50])
            pulse_start = 0
            for tick in range(1,len(waveform),1):
                    if abs(waveform[tick]-waveform[tick-1]) > 20:
                        pulse_start = tick-1
                        break
            waveform_adc_lead.append(waveform[pulse_start]-median)
        return waveform_adc_lead
    def Rise_time_calc(self):
        rise_time = []
        for keys in self.waveform_df:
            waveform = self.waveform_df[keys]
            median = np.median(waveform[0:50])
            i = np.argmax(waveform)
            if i == 0:
                rise_time.append(-1)
                continue
            #print(waveform[i]," ",median, " ")
            while waveform[i]>=median:
                i-=1
                if i == 0:
                    break
                #print(i)
            rise_time.append(np.argmax(waveform)-i)
        return rise_time
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
                #print("Waveform ledge goes all the way to the end")
                ledge_time.append(np.argmax(waveform))
                long_ledge +=1
                continue
            for t,tick in enumerate(waveform):
                if waveform[t] > 100 and waveform[t+1]-waveform[t] < 10 :
                    ledge_time_beg = t
                    break
            if abs(np.argmax(waveform)-ledge_time_beg) < 60:
                #print("No ledge")
                ledge_time.append(0)
                continue
            if number_0 >1:
                #print("No ledge")
                ledge_time.append(0)
                continue
            for k in range(ledge_time_beg, len(waveform),1):
                if waveform[k+1]-waveform[k] < -10:
                    ledge_time_end=0
                    ledge_time_beg = 0
                    #print("No ledge")
                    break
                if waveform[k+1]-waveform[k] > 10:
                    ledge_time_end = k
                    break
                if k == 3413:
                    ledge_time_end=0
                    ledge_time_beg = 0
                    break
            #print(ledge_time_end,ledge_time_beg)
            ledge_time.append(ledge_time_end-ledge_time_beg)
        
        return ledge_time,long_ledge
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

def threshold_info(waveform,value, threshold = None):
    channel_map = pd.read_csv("../datafiles/channel_mapping.txt", sep = " ")
    
    if threshold is None:
        return
    elif threshold == "Greater":
        for w,wire in enumerate(waveform):
            if wire >= float(value):
                print(channel_map[channel_map['LArSoft_ch'] == w])
    elif threshold == "Less":
        for w,wire in enumerate(waveform):
            if wire <= float(value):
                print(channel_map[channel_map['LArSoft_ch'] == w])
    elif threshold == "Equal":
        for w,wire in enumerate(waveform):
            if wire == float(value):
                print(channel_map[channel_map['LArSoft_ch'] == w])
    

def plot_wireplanes(event_number,run_number,metric,):#value,threshold
    value = None
    threshold = None
    #"Event":[],"UA Blob Max":[],"UA Blob Mean":[],"Long Ledge":[]
    return_list = [0,0,0,0,0]
    return_list[0] = event_number
    print(f"Run {run_number}, Event: {event_number}")
    if not os.path.exists(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/run_{run_number}/"):
        os.mkdir(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/run_{run_number}/")
    if not os.path.exists(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/run_{run_number}/{event_number}/"):
        os.mkdir(f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/run_{run_number}/{event_number}/")
    directory = f"/Users/danielcarber/Documents/SBND/Noise Analysis/Plots/Event_diagnostics/run_{run_number}/{event_number}/{metric}/"
    if not os.path.exists(directory):
        os.mkdir(directory)
    files =uproot.open(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/run{run_number}/waveform_{event_number}.root")
    #print(files['tpc_noise;1'].keys())
    print(metric)
    wire_df = load_wire_info()
    '''
    Waveform_df = {}
    ch_id = 0

    raw_rms = files['tpc_noise;1']['UB_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: UB First East Induction")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask]),"Min within range: ",np.min(waveform[mask]),"Mean within range: ",np.mean(waveform[mask]))

    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(0,ch_id,1):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            c.append(waveform[ch])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = f"{calc.type}",range_color=calc.range,title=f"East TPC First Ind Wire {metric} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)
    
    fig.write_image(directory+f'UB_plane_diagram_{event_number}_{metric}.png')
    #fig.show()
    
    
# VB Plane
    
    Waveform_df = {}
    ch_id = 1984
    raw_rms = files['tpc_noise;1']['VB_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)
    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: VB Second East Induction")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask]),"Min within range: ",np.min(waveform[mask]),"Mean within range: ",np.mean(waveform[mask]))
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(1984,ch_id,1):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-1984])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color =f"{calc.type}",range_color=calc.range,title=f"East TPC Second Ind Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f"VB_plane_diagram_{event_number}_{metric}.png")
    #fig.show()
    

#.      YB Plane    

    Waveform_df = {}
    ch_id = 3968
    raw_rms = files['tpc_noise;1']['YB_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: YB East Collection")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if sum(mask) > 0:
        print("Max within range: ",np.max(waveform[mask]),"Min within range: ",np.min(waveform[mask]),"Mean within range: ",np.mean(waveform[mask]))
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(3968,ch_id,1):
        #a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        #b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['y_0'][ch]),int(wire_df['y_1'][ch]),1):
            x_i = wire_df['z_0'][ch]
            x.append(wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-3968])
            y.append(i)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}

    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = f"{calc.type}",range_color=calc.range,title=f"East TPC Coll Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'YB_plane_diagram_{event_number}_{metric}.png')
    #fig.show()
    '''
    
    Waveform_df = {}
    ch_id = 5632
    raw_rms = files['tpc_noise;1']['UA_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)
    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: UA First West Induction")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    if metric == 'Ledge':
        print("Max within range: ",np.max(waveform[132:832]),"Min within range: ",np.min(waveform[132:832]),"Mean within range: ",np.mean(waveform[132:832]))

        return_list[1] = np.max(waveform[132:832])
        return_list[2] = np.mean(waveform[132:832])
        return_list[3] = calc.ledge
        print("Number of channels with ledges past event display: ",calc.ledge)
    else:
        print("Max within range: ",np.max(waveform),"Min within range: ",np.min(waveform),"Mean within range: ",np.mean(waveform))

        return_list[1] = np.max(waveform)
        return_list[2] = np.mean(waveform)
    event_info = files['Event_info;1']['Time'].array().to_list()
    print(event_info[0])
    return_list[4] = datetime.fromtimestamp(event_info[0], timezone(timedelta(hours=-7)))
        #return_list[3] = calc.ledge
        #print("Number of channels with ledges past event display: ",calc.ledge)
    #else:
    #    return_list[1] = -10
    #    return_list[2] = -10
    #    return_list[3] = -10
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(5632,ch_id,1):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-5632])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = f"{calc.type}",range_color=calc.range,title=f"West TPC First Ind Wire {str(metric)} Signal")


    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'UA_plane_diagram_{event_number}_{metric}.png')
    #fig.show()

    
    Waveform_df = {}
    ch_id = 7616
    raw_rms = files['tpc_noise;1']['VA_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: VA Second West Induction")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    #if sum(mask) > 0:
    if metric == 'Ledge':
        print("Max within range: ",np.max(waveform[132:832]),"Min within range: ",np.min(waveform[132:832]),"Mean within range: ",np.mean(waveform[132:832]))

        print("Number of channels with ledges past event display: ",calc.ledge)
    else:
        print("Max within range: ",np.max(waveform),"Min within range: ",np.min(waveform),"Mean within range: ",np.mean(waveform))

    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(7616,ch_id,1):
        a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['r'][ch])):
            x_i = i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch]
            x.append(i*((wire_df['z_1'][ch]-wire_df['z_0'][ch])/int(wire_df['r'][ch]))+wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-7616])
            y.append(a*x_i+b)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = f"{calc.type}",range_color=calc.range,title=f"West TPC Second Ind Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'VA_plane_diagram_{event_number}_{metric}.png')
    #fig.show()
                    
        
    Waveform_df = {}
    ch_id = 9600
    raw_rms = files['tpc_noise;1']['YA_plane'].array().to_list()
    for r,rms in enumerate(raw_rms):
        if r%3415 == 0:
            Waveform_df[f'Channel_{ch_id}'] = []
            Waveform_df[f'Channel_{ch_id}'].append(rms)
            ch_id +=1
        else:
            Waveform_df[f'Channel_{ch_id-1}'].append(rms)

    Waveform_df = pd.DataFrame(Waveform_df)
    calc = waveform_calc(Waveform_df,metric)
    waveform = calc.Run_calc()
    threshold_info(waveform, value = value, threshold = threshold)
    print("\nWire plane: YA West Collection")
    print("Max: ",max(waveform),"Min: ",min(waveform),"Mean: ",np.mean(waveform))
    waveform = np.array(waveform)
    mask = (waveform > calc.range[0]) & (waveform < calc.range[1])
    #if sum(mask) > 0:
    if metric == 'Ledge':
        print("Max within range: ",np.max(waveform[132:832]),"Min within range: ",np.min(waveform[132:832]),"Mean within range: ",np.mean(waveform[132:832]))

        print("Number of channels with ledges past event display: ",calc.ledge)
    else:
        print("Max within range: ",np.max(waveform),"Min within range: ",np.min(waveform),"Mean within range:", np.mean(waveform))
    #ch_id = 10
    color = []
    x = []
    y=[]
    c = []
    #fig = make_subplots(rows=1,cols=1)
    for ch in range(9600,ch_id,1):
        #a = (wire_df['y_1'][ch]-wire_df['y_0'][ch])/(wire_df['z_1'][ch]-wire_df['z_0'][ch])
        #b = wire_df['y_1'][ch] -a*wire_df['z_1'][ch]

        color.append(ch)
        for i in range(int(wire_df['y_0'][ch]),int(wire_df['y_1'][ch]),1):
            x_i = wire_df['z_0'][ch]
            x.append(wire_df['z_0'][ch])
            #print(ch)
            c.append(waveform[ch-9600])
            y.append(i)
    df = {'Z [cm]':np.array(x),'Y [cm]':np.array(y),f"{calc.type}":np.array(c)}
    #df = pd.DataFrame(df)
    fig=px.scatter(df,x="Z [cm]",y ="Y [cm]",color = f"{calc.type}",range_color=calc.range,title=f"West TPC Coll Wire {str(metric)} Signal")
    fig.update_layout(height = 800, width = 1200,showlegend = False)

    fig.write_image(directory+f'YA_plane_diagram_{event_number}_{metric}.png')
    #fig.show()

    
    return return_list
def main():
    
    #event_number = input("Please enter event number:")
    run_number = input("Please enter run number:")
    dir_list = os.listdir(f"/Users/danielcarber/Documents/SBND/Noise Analysis/data/run{run_number}/")
    #metric = input("Please enter calculation type (RMS, Max, Min,Integral, Rise or Range):")
    metric_list = ['Rise']#'Rise','Integral','Rise_time','Ledge'
    channel_map = pd.read_csv("../datafiles/channel_mapping.txt", sep = " ")
    df_metric = {}
    #print(dir_list)
    for metric in metric_list:
        df_metric = {}
        df_metric[f'Event'] = []
        df_metric[f'{metric}_max'] = []
        df_metric[f'{metric}_mean'] = []
        df_metric[f'Event_time'] = []
        if metric == 'Ledge':
            df_metric[f'Long Ledge'] = []
        for event in tqdm(dir_list):
            if event[:4] != "wave":
                continue
            event_number = event[9:]
            
            print("Event:",event_number)
            event_number = int(event_number[:-5])
            
            #print(event_number)
        #threshold = input(f"If you want to print threshold for {metric} enter type of threshold (Greater, Less, or Equal): ")
        #if (threshold != "" or threshold != None):
        #    value = input("Please enter value for threshold: ")

            metric_list =plot_wireplanes(event_number,run_number,str(metric),)#value,threshold
            df_metric['Event'].append(metric_list[0])
            df_metric[f'{metric}_max'].append(metric_list[1])
            df_metric[f'{metric}_mean'].append(metric_list[2])
            if metric == 'Ledge':
                df_metric[f'Long Ledge'].append(metric_list[3])
            df_metric['Event_time'].append(metric_list[4])
        print(df_metric)
        df_metric = pd.DataFrame(df_metric)
        df_metric.to_csv(f'/Users/danielcarber/Documents/SBND/Noise Analysis/data/Run_{run_number}_{metric}.csv', index=False) 
    
if __name__=="__main__": 
    main() 