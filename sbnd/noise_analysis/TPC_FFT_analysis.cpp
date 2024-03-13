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
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
//#include <bits/stdc++.h> 

using namespace std;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

vector<double> Hit_removal(vector<double> channels,float Pedestal){	
	vector<short> noise_channels;
	//for (int i = 0; i <=channels.GetSize();i++){
	//cout<<"Start of Hit Removal"<<endl;
		vector<double> ADCs = channels;
		float pedestal = Pedestal;
		vector<double> noise;
		noise.clear();
		for (int j = 0; j < channels.size();j++){
			//cout<<"Start of searching for hits"<<endl;
			double ADC = abs(ADCs.at(j)-pedestal);
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

double Noise_levels(vector<double> noise_channels){
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
vector<double> FFT(vector<double> noise_channel){
	//vector<float> FFT;
	int vec_size = noise_channel.size();
	Int_t size = vec_size;
	//vector<double> noise_vector(noise_channel.begin(), noise_channel.end());
	double inputSignalDouble[vec_size];
	cout<<"Turn vector into Double_t"<<endl;
    	for (size_t i = 0; i < vec_size; ++i) {
        	inputSignalDouble[i] = noise_channel[i];
    	}
	noise_channel.clear();
	//size_t* vectorSize = &size;
	//Int_t intVectorSize = static_cast<Int_t>(vectorSize);
	cout<<"Transforming"<<endl;
	//TVirtualFFT::SetTransform(0);
	cout<<"Size of waveform: "<<size<<endl;
   	TVirtualFFT* fft = TVirtualFFT::FFT(1, &size, "R2C ES K");
	if (!fft) {
        std::cerr << "Error: Failed to initialize FFT." << std::endl;
        return vector<double>();
    }
    
    fft->SetPoints(inputSignalDouble);
    fft->Transform();
	cout<<"Grabbing real and imag"<<endl;
	double fftReal=0;
        double fftImag=0;
	vector<double> fftMag(vec_size / 2 + 1);
	for(size_t k=1;k<vec_size / 2 + 1;k++){
	fft->GetPointComplex(k,fftReal, fftImag);
	//delete fft;
		fftMag[k] = TMath::Sqrt(fftReal*fftReal + fftImag*fftImag);
	}
	//float fftMag[vec_size / 2 + 1];
	cout<<"Getting magnitude"<<vec_size<<endl;
	//transform(fftReal.begin(),fftReal.end(),fftReal.begin(),fftReal.begin(),multiplies<Double_t>());
	//transform(fftImag.begin(),fftImag.end(),fftImag.begin(),fftImag.begin(),multiplies<Double_t>());
	//transform(fftReal.begin(),fftReal.end(),fftImag.begin(),fftMag.begin(),plus<Double_t>());
	//delete fft;
	cout<<"Finished"<<fftMag[0]<<endl;
	//fftReal.clear();
	//fftImag.clear();
	return fftMag;
		
}

void LoadRawDigits(TFile *inFile)
{	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	//Events.Print();
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	//vector<float> RMS_total(11264,0.0f);
	vector<vector<double>> FFT_total(11264,vector<double>(3415/2+1,0));
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		evt+=1;
		//for(int i = 0; i<myPedestal.GetSize();i++){
	//	cout<<myPedestal.GetSize()<<endl;
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<"Grabbed ADCs"<<endl;
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<short> channels;
		//vector<double> x(myADC.GetSize());
		for(int p=0; p<myADC.GetSize();p++){
                	channels.push_back(myADC[p].Channel());
        	}
		for(int ki=0; ki<11264;ki++){
			//cout<<myADC[ki].Channel()<<endl;

			auto in = find(channels.begin(),channels.end(), ki);
            		int index = in-channels.begin();
			vector<double> x(myADC[index].Samples(),0);
			cout<<getValue()<<endl;
			//cout<<"Starting noise removal"<<index<<endl;
			if (myADC[index].NADC() < 3000){
				continue;
			}
			cout<<getValue()<<endl;
			for (size_t itick=0; itick < myADC[index].Samples(); ++itick) x[itick] = myADC[index].ADC(itick);
			vector<double> channel_fft = FFT(x); 
			//cout<<"FFT is calculated"<<endl;
			cout<<getValue()<<endl;
			//transform(FFT_total[ki].begin(),FFT_total[ki].end(),channel_fft.begin(),FFT_total[ki].begin(),plus<double>());
			for (size_t k = 0; k < FFT_total[ki].size(); ++k) {
    			FFT_total[ki][k] += channel_fft[k]-FFT_total[ki][k]/evt;
    		}
			cout<<"FFT is calculated"<<FFT_total[ki].size()<<endl;
			cout<<getValue()<<endl;

			cout<<"Channel "<<ki<<endl;

		}

		
		//if (evt==4){
		//		for(int wire = 0; wire<FFT_total.size(); wire++){
		//			short mov_avg = 4;
		//			transform(FFT_total[wire].begin(),FFT_total[wire].end(),FFT_total[ch].begin(),[mov_avg](float &c){ return c/mov_avg; });
		//		}
		//	}
		//else if (evt%4==0){
		//
		//}
		cout<<"Event:"<<evt<<endl;
		//break;
	}

	
	TFile* file = new TFile("fft_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_rms;
	double avg_fft;
	tree->Branch("fft_mag", &avg_fft, "avg_fft/F");
	//tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	for(int ch = 0; ch<FFT_total.size(); ch++){
		transform(FFT_total[ch].begin(),FFT_total[ch].end(),FFT_total[ch].begin(),[evt](double &c){ return c/evt; });
		//for_each(FFT_total[ch].begin(), FFT_total[ch].end(), [evt](float &c){ c /= evt; });
		//avg_fft = FFT_total[ch];
		//cout<<FFT_total[ch][100]<<endl;
		TString bid = "fft_mag";
		bid +=ch;
		//tree->Branch(bid, &avg_fft, "avg_fft/F");
		for (size_t c = 0; c < FFT_total[ch].size(); ++c) {
                	avg_fft = FFT_total[ch][c];
			tree->Fill();
                }
		//tree->Fill();	
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

void TPC_FFT_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run11527/run_11527.root")
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
