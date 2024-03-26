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


		

void LoadRawDigits(TFile *inFile,int sel_evt)
{	
	//TTree *Events = (TTree*)inFile->Get("Events;1");
	//TString rootfilename(filename.c_str());	
	//TFile *inFile = TFile::Open(rootfilename.Data());	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//Events.Print();
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	TFile* file = new TFile("waveform_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		if (evt != sel_evt) {
			cout<<"Skip event: "<<evt<<endl;
			evt +=1;
			continue;
		}
		//for(int i = 0; i<myPedestal.GetSize();i++){
	//	cout<<myPedestal.GetSize()<<endl;
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<short> channels;
		
		float tick;
		
		tree->Branch("UB_plane", &tick,"tick/F");
		tree->Branch("VB_plane", &tick,"tick/F");
		tree->Branch("YB_plane", &tick,"tick/F");
		tree->Branch("UA_plane", &tick,"tick/F");
		tree->Branch("VA_plane", &tick,"tick/F");
		tree->Branch("YA_plane", &tick,"tick/F");
		for(int p=0; p<myADC.GetSize();p++){
			channels.push_back(myADC[p].Channel());
		}
		for(int ki=0; ki<11264;ki++){
        	if (ki < 1984) {
        		tree->SetBranchStatus("UB_plane", 1);
        		tree->SetBranchStatus("VB_plane", 0);
        		tree->SetBranchStatus("YB_plane", 0);
        		tree->SetBranchStatus("UA_plane", 0);
        		tree->SetBranchStatus("VA_plane", 0);
        		tree->SetBranchStatus("YA_plane", 0);
        	}
        	else if (ki < 3968){
        		tree->SetBranchStatus("UB_plane", 0);
        		tree->SetBranchStatus("VB_plane", 1);
        		tree->SetBranchStatus("YB_plane", 0);
        		tree->SetBranchStatus("UA_plane", 0);
        		tree->SetBranchStatus("VA_plane", 0);
        		tree->SetBranchStatus("YA_plane", 0);
        	}
        	else if (ki < 5632){
        		tree->SetBranchStatus("UB_plane", 0);
        		tree->SetBranchStatus("VB_plane", 0);
        		tree->SetBranchStatus("YB_plane", 1);
        		tree->SetBranchStatus("UA_plane", 0);
        		tree->SetBranchStatus("VA_plane", 0);
        		tree->SetBranchStatus("YA_plane", 0);
        	}
        	else if (ki < 7616){
        		tree->SetBranchStatus("UB_plane", 0);
        		tree->SetBranchStatus("VB_plane", 0);
        		tree->SetBranchStatus("YB_plane", 0);
        		tree->SetBranchStatus("UA_plane", 1);
        		tree->SetBranchStatus("VA_plane", 0);
        		tree->SetBranchStatus("YA_plane", 0);
        	}
        	else if (ki < 9600){
        		tree->SetBranchStatus("UB_plane", 0);
        		tree->SetBranchStatus("VB_plane", 0);
        		tree->SetBranchStatus("YB_plane", 0);
        		tree->SetBranchStatus("UA_plane", 0);
        		tree->SetBranchStatus("VA_plane", 1);
        		tree->SetBranchStatus("YA_plane", 0);
        	}
        	else if (ki < 11264){
        		tree->SetBranchStatus("UB_plane", 0);
        		tree->SetBranchStatus("VB_plane", 0);
        		tree->SetBranchStatus("YB_plane", 0);
        		tree->SetBranchStatus("UA_plane", 0);
        		tree->SetBranchStatus("VA_plane", 0);
        		tree->SetBranchStatus("YA_plane", 1);
        	}
			vector<float> x;
			int channel = myADC[ki].Channel();
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			//int in = index;
			cout<<"Index:"<<in<<", Channel:"<<myADC[in].Channel()<<", Loop index:"<<ki<<endl;
			if (myADC[in].NADC() != 3415){
				cout<<"Dead Channel"<<endl;
				for (size_t itick=0; itick < 3415; ++itick){
					tick = -1000;

					x.push_back(-1000);
					tree->Fill();
				}
				//continue;
			}
			else { 
			for (size_t itick=0; itick < 3415; ++itick){
				tick = myADC[in].ADC(itick);
				x.push_back(myADC[in].ADC(itick));
				tree->Fill();

			}
			
			//x.clear();
			} 
			cout<<x.size()<<endl;
			

		}
		evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;
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

void load_waveforms(TString inputFile="/pnfs/sbnd/scratch/users/jaz8600/Decoded/decoded_data_evb02_run12007_14_20240319T153634.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	int sel_evt = 18;
	LoadRawDigits(inFile,sel_evt);
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
