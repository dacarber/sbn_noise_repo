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
	vector<float> ADCs = channel;

	float pedestal = Pedestal;
	vector<float> noise;
	for (int j = 0; j < channel.size();j++){
		float ADC = (float) TMath::Abs(ADCs.at(j)-pedestal);
		if (ADC > 10.0){
			//noise.push_back(ADC);
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
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	vector<float> RMS_total(11264,0.0f);
	cout<<"Running Events"<<endl;
	float evt = 0.0;
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
		for(int ki=0; ki<11264;ki++){
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			cout<<"Channel: "<<myADC[in].Channel()<<endl;
			int channel = myADC[in].Channel();
			cout<<"Number of ticks: "<< myADC[in].NADC()<<endl;
			if (myADC[in].NADC() != 3415){
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
		cout<<"Event:"<<evt<<endl;
	}
	
	TFile* file = new TFile("noise_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_rms;
	tree->Branch("raw_rms", &avg_rms, "avg_rms/F");
	for(int ch = 0; ch<RMS_total.size(); ch++){
		avg_rms = RMS_total.at(ch)/evt;
		tree->Fill();	
	}
	
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;

}

void TPC_Noise_analysis(TString inputFile="/exp/sbnd/data/users/trj/run11505/tpcdecode_data_evb03_run11505_24_20240304T182922.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}

