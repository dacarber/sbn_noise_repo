#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include "TMath.h"
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
#include "lardataobj/RawData/RawDigit.h"

//#include <bits/stdc++.h> 

using namespace std;

vector<float> Hit_removal(vector<float> channel,float Pedestal){	
	vector<float> ADCs = channel;

	float pedestal = Pedestal;
	vector<float> noise;
	int skips =0;
	for (int j = 0; j < channel.size();j++){
		//cout<<"ADC with pedestal: "<<ADCs.at(j)<<endl;
		float ADC = (float) ADCs.at(j)-pedestal;
		//float ADC = (float) TMath::Abs(ADCs.at(j)-TMath::Median(channel.size(),channel.data()));
		//cout<<"ADC without pedestal: "<<ADC<<endl;
		if (TMath::Abs(ADC) > 10.0 || skips > 0){
			//noise.push_back(ADC);
			if (skips > 0){
				skips-=1;
				continue;
			}
			skips = 50;
			for (int i=0; i<50;i++){
				if (noise.size() == 0) continue;
				noise.pop_back();
			}
			continue;
		}
		else{
			noise.push_back(ADC);
		}
	
	}
	return noise;

}

float Noise_levels(vector<float> noise_channels){
	float RMS;
	float square;
	float sum;
	float mean =TMath::Mean(noise_channels.begin(),noise_channels.end());
	//float mean = accumulate(noise_channels.begin(),noise_channels.end(),0.0f)/noise_channels.size();
	for (int i = 0; i<noise_channels.size();i++){
		square = (noise_channels[i]-mean) * (noise_channels[i]-mean);
		sum = sum+square;
	}
	//RMS = TMath::RMS(noise_channels.begin(),noise_channels.end());
	//cout<<sum<<endl;
	//float mean = sum/noise_channels.size();
	
	RMS =sqrt(sum/noise_channels.size());
	//}
	cout<<"Size:"<<mean<<endl;
	return RMS;	
}


void LoadRawDigits(TFile *inFile)
{	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__DECODE.obj"); //New files with full decode
	TTreeReaderValue<unsigned int> event_info(Events, "EventAuxiliary.id_.event_");
	//TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj"); //For Data
	//TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_simtpc2d_daq_DetSim.obj"); //For MC
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_sptpc2d_raw_WCLSNF.obj"); //For Data Noise filter

	vector<float> RMS_total(11264,0.0f);
	int event_len = 3427;
	cout<<"Running Events"<<endl;
	float evt = 0.0;
	int true_evt = 0;
	while (Events.Next())
	{
		//if (evt == 0.0){
		//	evt+=1;
		//	true_evt+=1;
		//	continue;
		//}
		//for(int i = 0; i<myPedestal.GetSize();i++){
	//	cout<<myPedestal.GetSize()<<endl;
		unsigned int *event_num = event_info.Get();
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<float> channels;
		for(int p=0; p<myADC.GetSize();p++){
			channels.push_back(myADC[p].Channel());
		}
		int NhighBurst = 0;
		int NLowBurst = 0;
		bool burst = false;
		for(int ki=0; ki<11264;ki++){
			bool burst_high = false;
			bool burst_low = false;
			int channel = myADC[ki].Channel();
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			if (myADC[in].Samples() != event_len){//5995 for long readout windows, 3415 for standard readout windows
				continue;
			} 
			vector<short> x(myADC[in].Samples(),0); //Makes a vector the size of the uncompressed channel
			for (size_t itick=0; itick < myADC[in].Samples(); ++itick){
				if (myADC[in].ADC(itick) - myADC[in].GetPedestal() > 1000 && burst_high == false){					
					NhighBurst++;
					burst_high = true;
				}
				if (myADC[in].ADC(itick) - myADC[in].GetPedestal() < -1000 && burst_low == false){
					NLowBurst++;
					burst_low = true;
				}
				if (burst_low == true && burst_high == true){ 
					break;
				}
			}
			if (NhighBurst > 1500 && NLowBurst > 100 && NhighBurst > NLowBurst){ //|| NLowBurst > 1500
				
				burst = true;
				//cout<< "Burst: "<< NhighBurst<<" "<<NLowBurst<<endl;
				//cout<< "Event: "<< *event_num<< "Local event: "<<evt<<endl;
			}
		}
		
		if (burst == true){
			true_evt+=1;
			cout<<"Skip event: "<<*event_num<<endl;
			continue;
		}
		for(int ki=0; ki<11264;ki++){
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			cout<<"Channel: "<<myADC[in].Channel()<<endl;
			int channel = myADC[in].Channel();
			cout<<"Number of ticks: "<< myADC[in].NADC()<<endl;
			if (myADC[in].NADC() != event_len){ //5995 Long readout windows, 3415 for normal readout windows 
				RMS_total[ki] =  0.0;
				continue;

			}

			cout<<"Index:"<<in<<", Channel:"<<myADC[in].Channel()<<", Loop index:"<<ki<<endl;
			vector<float> x(myADC[in].Samples(),0.0f);
			for (size_t itick=0; itick < myADC[in].Samples(); ++itick) x[itick] =myADC[in].ADC(itick);
			
			
			vector<float> noise_channels = Hit_removal(x,myADC[in].GetPedestal());
			
			cout<<"Completed hit removal"<<endl;
			float RMS = Noise_levels(noise_channels);
			
			cout<<"RMS:"<<RMS<<endl;
			RMS_total[ki] =  RMS_total.at(ki)+RMS;
		}
		evt+=1.0;
		true_evt+=1;
		cout<<"Event:"<<evt<<"True Event num"<<true_evt<<endl;
	}
	
	TFile* file = new TFile("noise_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_rms;
	tree->Branch("raw_rms", &avg_rms, "avg_rms/F");
	for(int ch = 0; ch<RMS_total.size(); ch++){
		avg_rms = RMS_total.at(ch)/(evt-1);
		tree->Fill();	
	}
	
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;

}

//void TPC_Noise_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run14275/run_14275.root")
void TPC_Noise_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run17470/data_evb01_EventBuilder1_art1_run17470_112_20241027T005545_tpcdecode_WCLSNF-20241206T030417.root")
//void TPC_Noise_analysis(TString inputFile="/pnfs/sbn/data_add/sbnd/commissioning/run14401_decoded/decode_data_evb03_EventBuilder3_art4_run14401_14_20240704T014829-d07546f2-c49c-47f0-bc39-17d3d2f4226a.root")

{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}

