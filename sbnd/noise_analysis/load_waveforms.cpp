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
#include "TTreeReaderValue.h"
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

using namespace std;


		

void LoadRawDigits(TFile *inFile,int sel_evt)
{	
	//TTree *Events = (TTree*)inFile->Get("Events;1");
	//TString rootfilename(filename.c_str());	
	//TFile *inFile = TFile::Open(rootfilename.Data());	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//Events.Print();
	
	TTreeReaderValue<unsigned int> event_info(Events, "EventAuxiliary.id_.event_");
	TTreeReaderValue<unsigned int> time(Events, "EventAuxiliary.time_.timeHigh_");
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		unsigned int *event_num = event_info.Get();
		unsigned int *event_time = time.Get();
		cout<<"Event id: "<<*event_num<<endl;
		evt +=1;
		if (*event_num == 0 || *event_num >50){
            continue;
        }
		//for(int i = 0; i<myPedestal.GetSize();i++){
	//	cout<<myPedestal.GetSize()<<endl;
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<short> channels;
		
		
		for(int p=0; p<myADC.GetSize();p++){
			channels.push_back(myADC[p].Channel());
		}
		cout<<1<<endl;
		int NhighBurst = 0;
		int NLowBurst = 0;
		bool burst = false;
		for(int ki=0; ki<11264;ki++){
			bool burst_high = false;
			bool burst_low = false;
			int channel = myADC[ki].Channel();
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			if (myADC[in].Samples() != 3415){
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
				cout<< "Burst: "<< NhighBurst<<" "<<NLowBurst<<endl;
				cout<< "Event: "<< *event_num<< "Local event: "<<evt<<endl;
			}
		}
		
		if (burst == false){
			cout<<"Skip event: "<<*event_num<<endl;
			continue;
		}
		TString filename = "waveform_";
		filename += Form("%u", *event_num);
		filename += ".root";
		TFile* file = new TFile(filename, "RECREATE");
		TTree* tree = new TTree("tpc_noise", "tpc_noise");
		TTree* event_tree = new TTree("Event_info","Event_info")
		short tick;
		event_tree->Branch("Event",&event_num,"event_num/I")
		event_tree->Branch("Time",&event_time,"event_time/I")
		tree->Branch("UB_plane", &tick,"tick/S");
		tree->Branch("VB_plane", &tick,"tick/S");
		tree->Branch("YB_plane", &tick,"tick/S");
		tree->Branch("UA_plane", &tick,"tick/S");
		tree->Branch("VA_plane", &tick,"tick/S");
		tree->Branch("YA_plane", &tick,"tick/S");
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
			int total_tick = 0;
			if (myADC[in].NADC() != 3415){
				cout<<"Dead Channel"<<endl;
				for (size_t itick=0; itick < 3415; ++itick){
					tick = -1000;
					total_tick+=1;
					x.push_back(-1000);
					tree->Fill();
				}
				//continue;
			}
			else { 
			for (size_t itick=0; itick < 3415; ++itick){
				tick = myADC[in].ADC(itick);
				total_tick+=1;
				x.push_back(tick);
				tree->Fill();

			}
			cout<<tick<<" : "<<total_tick<<endl;
			
			//x.clear();
			} 
			cout<<x.size()<<endl;
			

		}
		//evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;
		file->Write();
		file->Close();
	}
	
	
	
	
	cout<<"Got ADC and Pedestal"<<endl;

}

void load_waveforms(TString inputFile="/pnfs/sbnd/scratch/users/jaz8600/Decoded/decoded_data_evb02_run12007_14_20240319T153634.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	int sel_evt = 18;
	LoadRawDigits(inFile,sel_evt);
}

