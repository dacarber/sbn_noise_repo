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
	cout<<"Noise: "<<noise[0]<<endl;
	return noise;

}

double Noise_levels(vector<short> noise_channels){
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
vector<short> Coherent_RMS(vector<vector<short>> noise_group){
	vector<short> waveform;
	short tick;
	for (int i = 0; i<noise_group[0].size();i++){
		vector<short> ADCs;
		
		for (int j = 0; j<noise_group.size();j++){
			ADCs.push_back(noise_group[j][i]);
			//cout<<"ADC Size: "<<noise_group.size()<<endl;
		}
		sort(ADCs.begin(), ADCs.end());
    	if (ADCs.size()+1 % 2 == 0) { // Even number of elements
        	tick = (ADCs[ADCs.size() / 2 - 1] + ADCs[ADCs.size()/ 2]) / 2.0;
		cout<<"Tick"<<tick<< endl;
    	} else { // Odd number of elements
        	tick = ADCs[ADCs.size() / 2];
    	}
		//cout<<"Tick"<<tick<< endl;
		waveform.push_back(tick);
	}
	//cout<<"Coh ADC "<<noise_group[0][0]<<endl;
	return waveform;
}


void LoadRawDigits(TFile *inFile)
{	
	//TTree *Events = (TTree*)inFile->Get("Events;1");
	//TString rootfilename(filename.c_str());	
	//TFile *inFile = TFile::Open(rootfilename.Data());	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//Events.Print();
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	//TTreeReaderArray<int> myADC(Events, "raw::RawDigits_daq__DetSim.obj.fADC");

	//TTreeReaderArray<Float_t> myPedestal(Events, "raw::RawDigits_daq__DECODER.obj.fPedestal");
	//vector<short> ADC;
	//vector<uint32_t> Pedestal;
	//cout<<myADC.GetSize()<<endl;
	//size_t channel_size = 2000;
	vector<float> RMS_total(11264,0.0f);
	vector<vector<short>> RMS_wave_total(352,vector<double>(3415,0));
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		
		//for(int i = 0; i<myPedestal.GetSize();i++){
	//	cout<<myPedestal.GetSize()<<endl;
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<vector<short>> channel_group;
		vector<short> noise_channels(ADC.size(),0);
		vector<short> channels;
		bool responsive_channel = true;
		short group_size = 32;
		for(int p=0; p<myADC.GetSize();p++){		//Puts all of the channel ids into a vector in the order the files have the events
            channels.push_back(myADC[p].Channel()); 
        }
		for(int ki=0; ki<myADC.GetSize();ki++){
			//cout<<ki<<endl;
			auto in = find(channels.begin(),channels.end(), ki); //finds the the location of the channel corresponding to ki
            int index = in-channels.begin();
			vector<double> x(myADC[index].Samples(),0); //Makes a vector the size of the uncompressed channel

			for (size_t itick=0; itick < myADC[index].Samples(); ++itick) x[itick] = myADC[index].ADC(itick);

			//Checks if the channel is dead
			if (myADC[index].NADC() < 3000 && (ki+1)%group_size != 0){
				responsive_channel = false;
				continue;
			}
			else if(myADC[index].NADC() < 3000 && (ki+1)%group_size == 0){
				vector<short> coherent_waveform = Coherent_RMS(channel_group);
				float Coh_RMS = Noise_levels(coherent_waveform);
				channel_group.clear();
				cout<<"Coh RMS:"<<Coh_RMS<<endl;
				for (int kh=0; kh < group_size; kh++){
					RMS_total[ki-kh] =  RMS_total.at(ki-kh)+Coh_RMS;
				}
				transform(RMS_wave_total[ki].begin(),RMS_wave_total[ki].end(),coherent_waveform.begin(),RMS_wave_total[ki].begin(),plus<float>());

				continue;
			}
			noise_channels = Hit_removal(x,myADC[index].GetPedestal());
			cout<<"Completed noise removal"<<channel<<endl;

			if ((ki+1)%group_size == 0 && responsive_channel == true){
				vector<short> coherent_waveform = Coherent_RMS(channel_group);
				float Coh_RMS = Noise_levels(coherent_waveform);
				channel_group.clear();
				cout<<"Coh RMS:"<<Coh_RMS<<endl;
				for (int kh=0; kh < group_size; kh++){
					RMS_total[ki-kh] =  RMS_total.at(ki-kh)+Coh_RMS;
				}
				transform(RMS_wave_total[ki].begin(),RMS_wave_total[ki].end(),coherent_waveform.begin(),RMS_wave_total[ki].begin(),plus<short>());
				
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
	tree.Branch("raw_rms", &avg_rms, "avg_rms/F");
	//tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	for(int ch = 0; ch<RMS_total.size(); ch++){
		avg_rms = RMS_total.at(ch)/evt;
		//for (size_t c = 0; c < FFT_total[ch].size(); ++c) {
                //	avg_FFT[c] = FFT_total[ch][c]/evt;
                //}
		tree.Fill();	
	}
	vector<short> coh_wave;
	tree->Branch("coh_wave", &coh_wave, "coh_wave/F");
	for(int ch = 0; ch<RMS_wave_total.size(); ch++){
		transform(RMS_wave_total[ch].begin(),RMS_wave_total[ch].end(),RMS_wave_total[ch].begin(),[evt](double &c){ return c/evt; });
		coh_wave = RMS_wave_total[ch]
		//for (size_t c = 0; c < RMS_wave_total[ch].size(); ++c) {
         //       	avg_fft = RMS_wave_total[ch][c];
			tree.Fill();
                }
	}
	file.Write();
	file.Close();
	
	cout<<"Got ADC and Pedestal"<<endl;

	/*TBranch* RawDigits = Events->GetBranch("raw::RawDigits_daq__DetSim.obj");
	cout<<"Got Branch"<<endl;
	TLeaf* ADC = RawDigits->GetLeaf("fADC");
	int entries = Events->GetEntries("EventAuxiliary.id_.event_");

	cout<<entries<<endl;
	for (int i=0;i<entries;i++){
		cout<<RawDigits->GetRow(i)<<endl;
		cout<<ADC->GetValue(i)<<endl;
	}*/
	//TBranch* RawDigits = RawDigits_branch->GetSubBranch('raw::RawDigits_daq__DetSim.obj');


	//Double_t ADC = RawDigits->GetLeaf("fADC")->GetValue(1);
	

	//cout<<ADC<<endl;

	//TH1F* hist = new TH1F('fADC');
	//Events->Draw('raw::RawDigits_daq__DetSim.obj.fADC');

}

void TPC_coherent_noise(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run11552/run_11552.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}
/*void Hit_removal(auto channels)
{
	for (int i = 0; i <=channels.GetSize();i++){
		vector<short> noise_channels;
		int sum = accumulate(channels[i].begin(), channels[i].end(), 0);
  		double mean = double(sum) / channels[i].size();
		cout<<"Length of channel and mean before:"+channels[i].size()+mean<<endl;
		for (int j = 0; j =< 34;j++){
			double max_sig = *max_element(channels[i][j*100:(j+1)*100]);
			double min_sig = *min_element(channels[i][j*100:(j+1)*100]);
			double sig_diff = max_sig-min_sig;
			if sig_diff > 30{
				noise_channels.insert(i*100,vector<short> zeros(100,0.0));
			noise_channels.insert(i*100,channels[i][j*100:(j+1)*100]);
			}
		}
		int sum = accumulate(noise_channels.begin(), noise_channels.end(), 0);
  		mean = double(sum) / noise_channels.size();
		channels[i] = noise_channels;
		cout<<"Length of channel and mean after:"+noise_channels[i].size()+mean<<endl;
	}
	cout<<"Number of channels"+noise_channels[i].size()<<endl;

}*/
