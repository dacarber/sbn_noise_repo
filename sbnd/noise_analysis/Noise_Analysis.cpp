#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCollection.h"
#include <iostream>
#include "THStack.h"
#include "TAxis.h"
#include "TGraphErrors.h"
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
//#include <bits/stdc++.h> 

using namespace std;

void Hit_removal(TTreeReaderArray<vector<short>> channels)
{
	for (int i = 0; i <=channels.GetSize();i++){
		int noise_channels[3400] = channels[i];
		int sum = accumulate(channels[i].begin(), channels[i].end(), 0);
  		double mean = double(sum) / channels[i].size();
		cout<<"Length of channel and mean before:"+channels[i].size()+mean<<endl;
		for (int j = 0; j <= 3400;j+100){
			double max_sig = *max_element(channels[i][j:(j+1)]);
			double min_sig = *min_element(channels[i][j:(j+1)]);
			double sig_diff = max_sig-min_sig;
			if sig_diff > 30{
				noise_channels.insert(j*100,vector<short> zeros(100,0.0));
			noise_channels.insert(j*100,channels[i][j:(j+1)]);
			}
		}
		int sum = accumulate(noise_channels.begin(), noise_channels.end(), 0);
  		mean = double(sum) / noise_channels.size();
		channels[i] = noise_channels;
		cout<<"Length of channel and mean after:"+noise_channels[i].size()+mean<<endl;
	}
	cout<<"Number of channels"+noise_channels[i].size()<<endl;

}

void LoadRawDigits(TFile *inFile)
{

	//TTree *Events = (TTree*)inFile->Get("Events;1");
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	Events.Print("raw::RawDigits_daq__DetSim.obj");
	TTreeReaderArray<vector<short>> myADC(Events, "raw::RawDigits_daq__DetSim.obj.fADC");
	//TTreeReaderArray<int> myADC(Events, "raw::RawDigits_daq__DetSim.obj.fADC");

	TTreeReaderArray<Float_t> myPedestal(Events, "raw::RawDigits_daq__DetSim.obj.fPedestal");
	vector<uint32_t> ADC;
	vector<uint32_t> Pedestal;
	while (Events.Next())
	{
		//for(int i = 0; i<myPedestal.GetSize();i++){
		cout<<myPedestal.GetSize()<<endl;
		cout<<myADC.GetSize()<<endl;
		cout<<myADC[1].size()<<endl;
		//cout<<typeid(myADC)<<endl;
		Hit_removal(myADC)
		
		//double ADC = (double)myADC[0];
		//TH1S *h1 = new TH1S("h1", "h1 title", 100, 1.0, 10000.0);
		//for (int i = 0; i<3400;i++){
		//cout<<i<<endl;
		//h1->Fill(myADC[0].at(i));
	//}
		//h1->Draw();

	//}
	}
	
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

void Noise_Analysis(TString inputFile="/exp/sbnd/data/users/dcarber/gen_g4_detsim_reco1-61193afb-d94b-460d-a307-3dd6417a9b80_HitDumper-20230523T200522.root")
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