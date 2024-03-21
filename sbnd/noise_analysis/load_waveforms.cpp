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
/*float Coherent_RMS(vector<short> noise_channels){
	float coherent_RMS;
	int channel_index = 128
	for (int i = 0; i<noise_channels.size()/channel_index;i++){
		vector<short> channel_group(noise_channels.begin()+(i*channel_index),noise_channels.end()+((i+1)*channel_index));
		sort(channel_group.begin(),channel_group.end());
		coherent_RMS = (channel_group[64] + channel_group[63])/2.0;
	}
}*/
vector<float> FFT(vector<short> noise_channel){
	//vector<float> FFT;
	int vec_size = noise_channel.size();
	Int_t size = vec_size;
	//vector<double> noise_vector(noise_channel.begin(), noise_channel.end());
	std::vector<Double_t> inputSignalDouble(vec_size);
	cout<<"Turn vector into Double_t"<<endl;
    	for (size_t i = 0; i < vec_size; ++i) {
        	inputSignalDouble[i] = static_cast<Double_t>(noise_channel[i]);
    	}
	noise_channel.clear();
	//size_t* vectorSize = &size;
	//Int_t intVectorSize = static_cast<Int_t>(vectorSize);
	cout<<"Transforming"<<endl;
	TVirtualFFT::SetTransform(0);
   	TVirtualFFT* fft = TVirtualFFT::FFT(1, &size, "R2C");
	if (!fft) {
        	std::cerr << "Error: Failed to initialize FFT." << std::endl;
        	return vector<float>();
    	}
    	fft->SetPoints(inputSignalDouble.data());
    	fft->Transform();
	cout<<"Grabbing real and imag"<<endl;
	std::vector<Double_t> fftReal(vec_size / 2 + 1);
    	std::vector<Double_t> fftImag(vec_size / 2 + 1);
	fft->GetPoints(fftReal.data(), fftImag.data());
	delete fft;
	std::vector<float> fftMagnitude(vec_size / 2 + 1);
	cout<<"Getting magnitude"<<vec_size<<endl;
   	for (auto i = 0; i <fftReal.size(); ++i) {
		//cout<<"Find Mag"<<endl;
		//int j = static_cast<int>(i);
        	fftMagnitude[i] = static_cast<float>(sqrt(fftReal[0] * fftReal[0] + fftImag[0] * fftImag[0]));
		cout<<"Mag:"<<fftMagnitude[i]<<endl;
		fftReal.erase(fftReal.begin());
		fftImag.erase(fftImag.begin());
    	}
	cout<<"Finished mad"<<endl;
	//fftReal.clear();
	//fftImag.clear();
	return fftMagnitude;
		
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
	TFile* file = new TFile("waveform_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
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
		vector<short> channels;
		short tick;
		string branch_name = "waveform_" + to_string(evt); 
		tree->Branch(branch_name, &tick);
		for(int p=0; p<myADC.GetSize();p++){
			channels.push_back(myADC[p].Channel());
		}
		for(int ki=0; ki<myADC.GetSize();ki++){
			cout<<myADC[ki].Channel()<<endl;
			int channel = myADC[ki].Channel();
			auto index = find(channels.begin(),channels.end(), ki);
			int in = index-channels.begin();
			cout<<"Index:"<<in<<", Channel:"<<myADC[in].Channel()<<", Loop index:"<<ki<<endl;
			if (myADC[in].NADC() != 3415){
				continue;
			}
			for (size_t itick=0; itick < myADC[in].Samples(); ++itick){
				tick = myADC[in].ADC(itick);
				tree->Fill();

			} 
		file->Write();
			

		}
		evt+=1;
		cout<<"Event:"<<evt<<endl;
		//break;
	}
	
	
	
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

void load_waveforms(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run11665/run_11665.root")
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
