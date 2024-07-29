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



vector<double> FFT(vector<double> noise_channel){
	int vec_size = noise_channel.size();
	Int_t size = vec_size;
	double* inputSignalDouble = new double[vec_size];
    	for (size_t i = 0; i < vec_size; ++i) {
        	inputSignalDouble[i] = noise_channel[i];
    	}
	noise_channel.clear();
	//size_t* vectorSize = &size;
	//Int_t intVectorSize = static_cast<Int_t>(vectorSize);
	//cout<<"Transforming"<<endl;
	//TVirtualFFT::SetTransform(0);
	//cout<<"Size of waveform: "<<size<<endl;
   	TVirtualFFT* fft = TVirtualFFT::FFT(1, &size, "R2C ES K");
	if (!fft) {
        std::cerr << "Error: Failed to initialize FFT." << std::endl;
        return vector<double>();
    }
    
    fft->SetPoints(inputSignalDouble);
    fft->Transform();

	double fftReal=0;
        double fftImag=0;
	vector<double> fftMag(vec_size / 2 + 2);
	for(size_t k=1;k<vec_size / 2 + 2;k++){
		fft->GetPointComplex(k,fftReal, fftImag);
	//delete fft;
		if (k == vec_size / 2 + 1){
			fftMag[k] = 1;
			continue;
		}
		fftMag[k] = TMath::Sqrt(fftReal*fftReal + fftImag*fftImag);
	}


	delete[] inputSignalDouble;
	delete fft;
	//cout<<"Finished"<<fftMag[100]<<endl;
	return fftMag;
		
}

void LoadRawDigits(TFile *inFile)
{	
	cout<<"Got Events"<<endl;
	TTreeReader Events("Events;1", inFile);
	TTreeReaderArray<raw::RawDigit> myADC(Events, "raw::RawDigits_daq__TPCDECODER.obj");
	vector<vector<double>> FFT_total(11264,vector<double>(3415/2+2,0));
	cout<<"Running Events"<<endl;
	int evt = 0;
	while (Events.Next())
	{
		if (evt == 0){
			evt+=1;
			continue;
		}
		cout<<myADC.GetSize()<<endl; //Grabs the number of channels
		vector<short> ADC = myADC[1].ADCs();
		cout<<ADC.size()<<endl; //Grabs the number of time ticks
		vector<short> channels;
		for(int p=0; p<myADC.GetSize();p++){
                	channels.push_back(myADC[p].Channel());
        	}
		for(size_t ki=0; ki<11264;ki++){

			auto in = find(channels.begin(),channels.end(), ki);
            int index = find(channels.begin(),channels.end(), ki)-channels.begin();
			cout<<myADC[index].NADC()<<endl;
			if (myADC[index].NADC() != 3415){
				continue;
			}
			bool skip_channel = false;
			vector<double> x(myADC[index].Samples(),0);
			for (size_t itick=0; itick < myADC[index].Samples(); ++itick){ 
				if (abs(myADC[index].ADC(itick)-myADC[index].GetPedestal()) >  10){
					skip_channel = true;
					break;
				}
				x[itick] =myADC[index].ADC(itick);

			}
			if (skip_channel == true){
				continue;
			}
			vector<double> channel_fft = FFT(x);
			x.clear(); 
			transform(FFT_total[ki].begin(),FFT_total[ki].end(),channel_fft.begin(),FFT_total[ki].begin(),plus<double>());
	
			cout<<"FFT is calculated "<<FFT_total[ki][100]<<" "<<FFT_total[ki].capacity()<<" "<<x.capacity()<<endl;
			cout<<getValue()<<endl;

			//cout<<"Channel "<<ki<<endl;

		}

		evt +=1;		
		cout<<"Event:"<<evt<<endl;
	}

	
	TFile* file = new TFile("fft_output.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float avg_FFT;
	tree->Branch("avg_FFT", &avg_FFT, "avg_FFT/F");
	for(int ch = 0; ch<FFT_total.size(); ch++){
		//transform(FFT_total[ch].begin(),FFT_total[ch].end(),FFT_total[ch].begin(),[evt](double &c){ return c/evt; });
		cout<<FFT_total[ch][100]<<endl;

		for (size_t c = 0; c < FFT_total[ch].size(); ++c) {
			cout<<c<<" "<<FFT_total[ch][c]<<endl;
			avg_FFT = FFT_total[ch][c];
			tree->Fill();
                }
	}
	
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;
}

void TPC_FFT_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run14784/run_14784.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}

