//#include "TCanvas.h"
//#include "TStyle.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVectorT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <typeinfo>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <TH1F.h>
#include <TVirtualFFT.h>
#include <fstream>
#include <TChain.h>
#include <string>
//#include <bits/stdc++.h>

//#include <bits/stdc++.h> 

using namespace std;

vector<short> Hit_removal(vector<short> channels,float Pedestal){	
	vector<short> noise_channels;
	//for (int i = 0; i <=channels.GetSize();i++){
	//cout<<"Start of Hit Removal"<<endl;
		vector<short> ADCs = channels;
		float pedestal = Pedestal;
		vector<short> noise;
		noise.clear();
		for (int j = 0; j < channels.size();j++){
			//cout<<"Start of searching for hits"<<endl;
			short ADC = abs(ADCs.at(j)-pedestal);
			if (ADC > 10){
				noise.push_back(ADC);
				//continue;
			}
			else{
				noise.push_back(ADC);
			}
		
		}
		//noise_channels.push_back(noise);
	//}
	return noise;

}

double Noise_levels(vector<float> noise_channels){
	double RMS;
	float square;
	float sum;
	for (int i = 0; i<noise_channels.size();i++){
		square = noise_channels[i] * noise_channels[i];
		sum = sum+square;
	}
	//cout<<sum<<endl;
	//float mean = sum/noise_channels.size();
	//float mean = accumulate(noise_channels.begin(),noise_channels.end(),0.0f)/noise_channels.size();
	RMS =sqrt(sum/noise_channels.size());
	//}
	cout<<"Size:"<<sum<<endl;
	return RMS;	
}
vector<float> Coherent_RMS(vector<vector<short>> noise_group){
	vector<float> waveform;
	vector<float> mean_waveform(3415,0.0);
	short tick;
	int group_size = 32;
	

	for (int j = 0; j<noise_group.size();j++){
		cout<<noise_group[j].size()<<endl;
		transform(mean_waveform.begin(),mean_waveform.end(),noise_group[j].begin(),mean_waveform.begin(),plus<float>());
	}
	
	//for (int i = 0; i<noise_group[0].size();i++){
	//	vector<short> ADCs;
	//	
	//	for (int j = 0; j<noise_group.size();j++){
	//		ADCs.push_back(noise_group[j][i]);
	//		//cout<<"ADC Size: "<<noise_group.size()<<endl;
	//	}
	//	sort(ADCs.begin(), ADCs.end());

    //	if (ADCs.size()+1 % 2 == 0) { // Even number of elements
    //    	tick = (ADCs[ADCs.size() / 2 - 1] + ADCs[ADCs.size()/ 2]) / 2.0;
	//	cout<<"Tick"<<tick<< endl;
    //	} else { // Odd number of elements
    //   	tick = ADCs[ADCs.size() / 2];
    //	}
		//cout<<"Tick"<<tick<< endl;
	//	waveform.push_back(tick);
	//}
	transform(mean_waveform.begin(),mean_waveform.end(),mean_waveform.begin(),[group_size](float &c){ return c/group_size; });
	//cout<<"Coh ADC "<<noise_group[0][0]<<endl;
	return mean_waveform;
}
vector<float> Coh_removal(vector<short> noise, vector<float> coh_noise,float kh_length,float mid_length){
	vector<float> int_waveform(noise.size(),0);
	//short tick;
	cout<<"Check size"<<endl;
	if (noise.size() != 3415){
		return int_waveform;
	}
	cout<<"Transforming vector"<<endl;
	//transform(coh_noise.begin(),coh_noise.end(),coh_noise.begin(),[kh_length,mid_length](float &c){ return 0.0012502364558886422*kh_length+c-0.0012502364558886422*mid_length; });
	transform(noise.begin(),noise.end(),coh_noise.begin(),int_waveform.begin(),minus<float>());
	//cout<<"Coh ADC "<<noise_group[0][0]<<endl;
	cout<<"Returning vector"<<endl;
	return int_waveform;
}
float get_wire_length(int wire_number,vector<float> wire_lengths){
	
	return wire_lengths.at(wire_number);

}
void LoadRawDigits(TFile *inFile)
{	
	ifstream wire_length_file;
	float x_0;
	float y_0;
	float z_0;
	float x_1;
	float y_1;
	float z_1;
	float entry_string;
	wire_length_file.open("/exp/sbnd/app/users/dcarber/noise_analysis/sbn_noise_repo/sbnd/datafiles/Wire_lengths.txt");
	int entry = 0;
	string line;
	string ent;
	vector<float> wire_lengths;
	while(getline(wire_length_file, line)){
		istringstream ss(line);
		while(ss >>ent){
			cout<<ent<<endl;
			entry_string = stof(ent);
			entry +=1;
			if (entry%6==0){
				x_0 = float(entry_string);
			}
			if (entry%7==0){
				y_0 = float(entry_string);
			}
			if (entry%8==0){
				z_0 = float(entry_string);
			}
			if (entry%9==0){
				x_1 = float(entry_string);
			}
			if (entry%10==0){
				y_1 = float(entry_string);
			}
			if (entry%11==0){
				z_1 = float(entry_string);
				wire_lengths.push_back(sqrt(pow(x_1-x_0,2)+pow(y_1-y_0,2)+pow(z_1-z_0,2)));
			}
			entry_string = 0;
		}
	}
	cout<<"Got Wire lengths"<<wire_lengths.size()<<endl;
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//Events.Print();
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	vector<float> Int_RMS_total(11264,0.0f);
	vector<float> RMS_total(11264,0.0f);
	vector<vector<float>> RMS_wave_total(352,vector<float>(3415,0));
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		//vector<short> ADC = myADC[1].ADCs();
		vector<vector<short>> channel_group;
		bool responsive_channel = true;
		vector<short> channels;
		short group_size = 32;
		for(int p=0; p<myADC.GetSize();p++){		//Puts all of the channel ids into a vector in the order the files have the events
            channels.push_back(myADC[p].Channel()); 
        }
		for(int ki=0; ki<11264;ki++){
			//cout<<ki<<endl;
			auto in = find(channels.begin(),channels.end(), ki); //finds the the location of the channel corresponding to ki
            int index = in-channels.begin();
			responsive_channel = true;
			short channel = myADC[index].Channel();
			cout<<"Channel index: "<<index<<" Channel: "<< myADC[index].Channel()<<endl;
			cout<<"Channel size: "<<myADC[index].NADC()<<endl;
			//Checks if the channel is dead
			if (myADC[index].NADC() < 3000 && (ki+1)%group_size != 0){
				responsive_channel = false;
				continue;
			}
			else if(myADC[index].NADC() < 3000 && (ki+1)%group_size == 0){
				if (channel_group.size() == 0){
					channel_group.clear();
					continue;				
				}
				vector<float> coherent_waveform = Coherent_RMS(channel_group);
				float Coh_RMS = Noise_levels(coherent_waveform);
				channel_group.clear();
				cout<<"Coh RMS:"<<Coh_RMS<<endl;
				for (int kh=0; kh < group_size; kh++){
					RMS_total[channel-kh] =  RMS_total.at(channel-kh)+Coh_RMS;
				}
				transform(RMS_wave_total[channel/31].begin(),RMS_wave_total[channel/31].end(),coherent_waveform.begin(),RMS_wave_total[channel/31].begin(),plus<float>());

				continue;
			}
			vector<short> x(myADC[index].Samples(),0); //Makes a vector the size of the uncompressed channel
			for (size_t itick=0; itick < myADC[index].Samples(); ++itick) x[itick] = myADC[index].ADC(itick);
			vector<short> noise_channels = Hit_removal(x,myADC[index].GetPedestal());
			cout<<"Completed noise  "<<ki<<endl;

			if ((ki+1)%group_size == 0 && responsive_channel == true){
				channel_group.push_back(noise_channels);
				vector<float> coherent_waveform = Coherent_RMS(channel_group);
				float Coh_RMS = Noise_levels(coherent_waveform);
				
				cout<<"Coh RMS:"<<Coh_RMS<<endl;
				for (int kh=0; kh < group_size; kh++){
					RMS_total[channel-kh] =  RMS_total.at(channel-kh)+Coh_RMS;
					cout<<"1 channel in group: "<<kh<<endl;
					float kh_length = get_wire_length(channel-kh,wire_lengths);
					cout<<"2"<<endl;
					float mid_length = get_wire_length(channel-16,wire_lengths);
					cout<<"3"<<channel_group[31-kh].size()<<endl;
					if (channel_group[31-kh].size() != 3415){
						continue;
					}
					vector<float> int_wave = Coh_removal(channel_group[31-kh],coherent_waveform,kh_length,mid_length);
					cout<<"4"<<endl;
					float Int_RMS = Noise_levels(int_wave);
					cout<<"5"<<endl;
					Int_RMS_total[channel-kh] =  Int_RMS_total.at(channel-kh)+Int_RMS;
				}
				cout<<"6"<<endl;
				transform(RMS_wave_total[channel/31].begin(),RMS_wave_total[channel/31].end(),coherent_waveform.begin(),RMS_wave_total[channel/31].begin(),plus<short>());
				cout<<"combine waveform"<<endl;
				channel_group.clear();
			}
			else{
				//cout<<"Adding another channel "<<noise_channels[0]<<endl; 
				channel_group.push_back(noise_channels);
			}

		}

		evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;

	}
	

	TFile* file = new TFile("noise_output_coh.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_rms;
	//vector<float> avg_FFT;
	tree->Branch("raw_rms", &avg_rms, "avg_rms/F");
	//tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	for(int ch = 0; ch<RMS_total.size(); ch++){
		avg_rms = RMS_total.at(ch)/evt;
		//for (size_t c = 0; c < FFT_total[ch].size(); ++c) {
                //	avg_FFT[c] = FFT_total[ch][c]/evt;
                //}
		tree->Fill();	
	}
	short coh_wave;
	tree->Branch("coh_wave", &coh_wave, "coh_wave/F");
	for(int ch = 0; ch<RMS_wave_total.size(); ch++){
		transform(RMS_wave_total[ch].begin(),RMS_wave_total[ch].end(),RMS_wave_total[ch].begin(),[evt](float &c){ return c/evt; });
		//coh_wave = RMS_wave_total[ch];
		for (size_t c = 0; c < RMS_wave_total[ch].size(); ++c) {
                	coh_wave = RMS_wave_total[ch][c];
			tree->Fill();
                }
	}
	float int_rms;
	//vector<float> avg_FFT;
	tree->Branch("int_rms", &int_rms, "int_rms/F");
	//tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	for(int ch = 0; ch<Int_RMS_total.size(); ch++){
		int_rms = Int_RMS_total.at(ch)/evt;
		//for (size_t c = 0; c < FFT_total[ch].size(); ++c) {
                //	avg_FFT[c] = FFT_total[ch][c]/evt;
                //}
		tree->Fill();	
	}
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;



}

void TPC_Coh_removal(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run11799/run_11799.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}
