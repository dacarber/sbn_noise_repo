//#include "TCanvas.h"
//#include "TStyle.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
//#include "TLegend.h"
//#include "TCollection.h"
#include <iostream>
//#include "THStack.h"
//#include "TAxis.h"
//#include "TGraphErrors.h"
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

//#include <bits/stdc++.h> 

using namespace std;

vector<float> Hit_removal(vector<float> channel,float Pedestal){	
	vector<short> noise_channels;
	//for (int i = 0; i <=channels.GetSize();i++){
	//cout<<"Start of Hit Removal"<<endl;
		vector<float> ADCs = channel;

		float pedestal = Pedestal;
		vector<float> noise;
		for (int j = 0; j < channel.size();j++){
			//cout<<"Start of searching for hits"<<endl;
			float ADC = abs(ADCs.at(j)-pedestal);
			if (ADC > 10){
				//noise.push_back(ADC);
				continue;
			}
			else{
				noise.push_back(ADC);
			}
		
		}
		//noise_channels.push_back(noise);
	//}
	return noise;

}

float Noise_levels(vector<float> noise_channels){
	float RMS;
	float square;
	float sum;
	//for (int i = 0; i<noise_channels.size();i++){
	//	square = noise_channels[i] * noise_channels[i];
	//	sum = sum+square;
	//}
	RMS = TMath::RMS(noise_channels.begin(),noise_channels.end());
	//cout<<sum<<endl;
	//float mean = sum/noise_channels.size();
	//float mean = accumulate(noise_channels.begin(),noise_channels.end(),0.0f)/noise_channels.size();
	//RMS =sqrt((sum/noise_channels.size()) - mean*mean);
	//}
	cout<<"Size:"<<sum<<endl;
	return RMS;	
}


void LoadRawDigits(TFile *inFile)
{	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	vector<float> RMS_total(11264,0.0f);
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
		vector<float> channels;
		for(int p=0; p<myADC.GetSize();p++){
			channels.push_back(myADC[p].Channel());
		}
		for(int ki=0; ki<myADC.GetSize();ki++){
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			cout<<"Channel: "<<myADC[in].Channel()<<endl;
			int channel = myADC[in].Channel();
			cout<<"Number of ticks: "<< myADC[in].NADC()<<endl;
			if (myADC[in].NADC() != 3415){
				RMS_total[ki] =  0;
				continue;

			}

			cout<<"Index:"<<in<<", Channel:"<<myADC[in].Channel()<<", Loop index:"<<ki<<endl;
			vector<float> x(myADC[in].Samples(),0);
			for (size_t itick=0; itick < myADC[in].Samples(); ++itick) x[itick] =myADC[in].ADC(itick);
			
			
			vector<float> noise_channels = Hit_removal(x,myADC[in].GetPedestal());
			
			cout<<"Completed hit removal"<<endl;
			float RMS = Noise_levels(noise_channels);
			
			cout<<"RMS:"<<RMS<<endl;
			RMS_total[ki] =  RMS_total.at(ki)+RMS;
		}
		evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;
	}
	
	TFile* file = new TFile("noise_output.root", "RECREATE");
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

void TPC_Noise_analysis(TString inputFile="/exp/sbnd/data/users/trj/run11505/tpcdecode_data_evb03_run11505_24_20240304T182922.root")
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
