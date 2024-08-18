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
	short tick;
	cout<<noise_group.size()<<endl;
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
vector<float> Coh_removal(vector<short> noise, vector<float> coh_noise){
	vector<float> int_waveform(noise.size(),0);
	int sum = accumulate(noise.begin(),noise.end(),0);
	cout<<"Check size"<<endl;
	if (noise.size() != 3415 || sum == 0 ){
		return int_waveform;
	}
	transform(noise.begin(),noise.end(),noise.begin(),[](float x) {return x * x});
	transform(coh_noise.begin(),coh_noise.end(),coh_noise.begin(),[](float x) {return x * x});
	transform(noise.begin(),noise.end(),coh_noise.begin(),int_waveform.begin(),minus<float>());
	transform(int_waveform.begin(),int_waveform.end(),int_waveform.begin(),[](float x) {return sqrt(x)});
	//cout<<"Coh ADC "<<noise_group[0][0]<<endl;
	cout<<"Returning vector"<<endl;
	return int_waveform;
}

void LoadRawDigits(TFile *inFile)
{	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);

	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj"); //For Data
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_simtpc2d_daq_DetSim.obj"); //For MC


	vector<float> RMS_total(11264,0.0f); //Stores the Coherent noise levels for entire TPC
	vector<float> INT_RMS_total(11264,0.0f); //Stores the Intrinsic noise levels for entire TPC
	vector<int> Entries(11264,0.0f);
	vector<int> Int_Entries(11264,0.0f);
	//vector<vector<float>> RMS_wave_total(352,vector<float>(3415,0));
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		

		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<vector<short>> channel_group;
		vector<vector<short>> int_channel_group;
		vector<short> noise_channels(ADC.size(),0);
		bool responsive_channel = true;
		vector<short> channels;
		short group_size = 8;

		//Puts all of the channel ids into a vector in the order the files have the events
		for(int p=0; p<myADC.GetSize();p++){		
            channels.push_back(myADC[p].Channel()); 
        }

        //Goes over all of the channels and does the analysis
		for(int ki=0; ki<11264;ki++){
			auto in = find(channels.begin(),channels.end(), ki); //finds the the location of the channel corresponding to ki
            int index = in-channels.begin();
			responsive_channel = true;
			short channel = myADC[index].Channel();
			cout<<"Channel index: "<<index<<" Channel: "<< myADC[index].Channel()<<endl;
			cout<<"Channel size: "<<myADC[index].NADC()<<endl;


			//Checks if the channel is dead
			if (myADC[index].Samples() != 3415 && (ki+1)%group_size != 0){
				responsive_channel = false;
				int_channel_group.push_back(vector<short>(3415,0))
				continue;
			}
			else if(myADC[index].Samples() != 3415 && (ki+1)%group_size == 0){
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
					Entries[channel-kh] = Entries.at(channel-kh)+1;
					if (accumulate(int_channel_group[kh].begin(),int_channel_group[kh].end(),0) == 0){
						continue;
					}
					vector<float> intrinsic_waveform = Coh_removal(int_channel_group[kh],coherent_waveform);
					float Int_RMS = Noise_levels(intrinsic_waveform);
					INT_RMS_total[channel-kh] = INT_RMS_total.at(channel-kh)+Int_RMS;
					Int_Entries[channel-kh] = Int_Entries.at(channel-kh)+1;
				}
				int_channel_group.clear();
				//transform(RMS_wave_total[channel/7].begin(),RMS_wave_total[channel/7].end(),coherent_waveform.begin(),RMS_wave_total[channel/31].begin(),plus<float>());

				continue;
			}

			//If channel is responsive the channel will grab the noise 
			bool skip_channel = false;
			vector<short> x(myADC[index].Samples(),0);
			for (size_t itick=0; itick < myADC[index].Samples(); ++itick){ 
				if (abs(myADC[index].ADC(itick)-myADC[index].GetPedestal()) >  20){
					skip_channel = true;
					int_channel_group.push_back(vector<short>(3415,0))
					break;
				}
				x[itick] = myADC[index].ADC(itick)-myADC[index].GetPedestal();

			}

			if ((ki+1)%group_size == 0 && responsive_channel == true){
				if (skip_channel == true){
					if (channel_group.size() == 0){
					channel_group.clear();
					int_channel_group.clear();
					continue;				
				}
					vector<float> coherent_waveform = Coherent_RMS(channel_group);
					float Coh_RMS = Noise_levels(coherent_waveform);
					channel_group.clear();
					cout<<"Coh RMS:"<<Coh_RMS<<endl;
					for (int kh=0; kh < group_size; kh++){
						RMS_total[channel-kh] =  RMS_total.at(channel-kh)+Coh_RMS;
						Entries[channel-kh] = Entries.at(channel-kh)+1;
						if (accumulate(int_channel_group[kh].begin(),int_channel_group[kh].end(),0) == 0){
							continue;
						}
						vector<float> intrinsic_waveform = Coh_removal(int_channel_group[kh],coherent_waveform);
						float Int_RMS = Noise_levels(intrinsic_waveform);
						INT_RMS_total[channel-kh] = INT_RMS_total.at(channel-kh)+Int_RMS;
						Int_Entries[channel-kh] = Int_Entries.at(channel-kh)+1;
					}
					int_channel_group.clear();
				}
				else{
					if (channel_group.size() == 0){
					channel_group.clear();
					int_channel_group.clear();
					continue;				
				}
					channel_group.push_back(x);
					int_channel_group.push_back(x);
					vector<float> coherent_waveform = Coherent_RMS(channel_group);
					float Coh_RMS = Noise_levels(coherent_waveform);
					channel_group.clear();
					cout<<"Coh RMS:"<<Coh_RMS<<endl;
					for (int kh=0; kh < group_size; kh++){
						RMS_total[channel-kh] =  RMS_total.at(channel-kh)+Coh_RMS;
						Entries[channel-kh] = Entries.at(channel-kh)+1;
						if (accumulate(int_channel_group[kh].begin(),int_channel_group[kh].end(),0) == 0){
							continue;
						}
						vector<float> intrinsic_waveform = Coh_removal(int_channel_group[kh],coherent_waveform);
						float Int_RMS = Noise_levels(intrinsic_waveform);
						INT_RMS_total[channel-kh] = INT_RMS_total.at(channel-kh)+Int_RMS;
						Int_Entries[channel-kh] = Int_Entries.at(channel-kh)+1;
					}
					int_channel_group.clear();
				}
				
				//transform(RMS_wave_total[channel/7].begin(),RMS_wave_total[channel/7].end(),coherent_waveform.begin(),RMS_wave_total[channel/31].begin(),plus<short>());
				//cout<<"combine waveform"<<endl;
			}
			else{
				if (skip_channel == true){ 
					continue;
				}
				//cout<<"Adding another channel "<<noise_channels[0]<<endl; 
				channel_group.push_back(x);
				int_channel_group.push_back(x);
			}

		}

		evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;
	}
	
	TFile* file = new TFile("noise_output_coh.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_rms;
	int entries;
	float int_rms;
	int int_entries;
	//vector<float> avg_FFT;
	tree->Branch("coh_rms", &avg_rms, "avg_rms/F");
	tree->Branch("entries", &entries, "entries/I");
	tree->Branch("int_rms", &int_rms, "int_rms/F");
	tree->Branch("int_entries", &int_entries, "int_entries/I");
	//tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	tree->SetBranchStatus("coh_rms", 1);
    tree->SetBranchStatus("entries", 0);
    tree->SetBranchStatus("int_rms", 0);
    tree->SetBranchStatus("int_entries", 0);
	for(int ch = 0; ch<RMS_total.size(); ch++){
		avg_rms = RMS_total.at(ch)/Entries.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("coh_rms", 0);
    tree->SetBranchStatus("entries", 1);
    tree->SetBranchStatus("int_rms", 0);
    tree->SetBranchStatus("int_entries", 0);
	for(int ch = 0; ch<Entries.size(); ch++){
		entries = Entries.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("int_rms", 1);
    tree->SetBranchStatus("int_entries", 0);
    tree->SetBranchStatus("coh_rms", 0);
    tree->SetBranchStatus("entries", 0);
	for(int ch = 0; ch<INT_RMS_total.size(); ch++){
		int_rms = INT_RMS_total.at(ch)/Int_Entries.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("coh_rms", 0);
    tree->SetBranchStatus("entries", 1);
    tree->SetBranchStatus("int_rms", 0);
    tree->SetBranchStatus("int_entries", 0);
	for(int ch = 0; ch<Int_Entries.size(); ch++){
		int_entries = Int_Entries.at(ch);

		tree->Fill();	
	}

	//short coh_wave;
	//tree->Branch("coh_wave", &coh_wave, "coh_wave/F");
	//for(int ch = 0; ch<RMS_wave_total.size(); ch++){
	//	transform(RMS_wave_total[ch].begin(),RMS_wave_total[ch].end(),RMS_wave_total[ch].begin(),[evt](float &c){ return c/evt; });
		//coh_wave = RMS_wave_total[ch];
	//	for (size_t c = 0; c < RMS_wave_total[ch].size(); ++c) {
    //           	coh_wave = RMS_wave_total[ch][c];
	//		tree->Fill();
    //            }
	//}
	file->Write();
	file->Close();
	
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

void TPC_coherent_noise(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run11799/run_11799.root")
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
