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


double Noise_levels(vector<double> noise_channels){
	double RMS;
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
vector<double> Coherent_RMS(vector<vector<short>> noise_group){
	vector<double> waveform;
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
vector<double> Coh_removal(vector<short> noise, vector<double> coh_noise){
	vector<double> int_waveform(noise.size(),0);
	float coh;
	float raw;
	int sum = accumulate(noise.begin(),noise.end(),0);
	cout<<"Check size"<<noise.size()<<endl;
	if (noise.size() != 3427 || sum == 0 ){
		return int_waveform;
	}

	transform(noise.begin(),noise.end(),coh_noise.begin(),int_waveform.begin(),minus<float>());
	cout<<"Returning vector"<<int_waveform[100]<<endl;

	return int_waveform;
}
vector<double> FFT(vector<double> noise_channel){
	cout<<"Starting FFT"<<endl;
	int vec_size = noise_channel.size();
	Int_t size = vec_size;
	double* inputSignalDouble = new double[vec_size];
    	for (size_t i = 0; i < vec_size; ++i) {
        	inputSignalDouble[i] = noise_channel[i];
    	}
	noise_channel.clear();

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
	cout<<"Finished"<<fftMag[100]<<endl;
	return fftMag;
}
float median(vector<float> &vec) {
    int n = vec.size();

    // Sort the vector
    sort(vec.begin(), vec.end());

    if (n % 2 == 0) { 
        // If the vector has an even number of elements, return the average of the middle two elements
        return (vec[n / 2 - 1] + vec[n / 2]) / 2.0;
    } else { 
        // If the vector has an odd number of elements, return the middle element
        return vec[n / 2];
    }
}

void LoadRawDigits(TFile *inFile)
{	
	int TOTAL_EVT=300;
	int event_len = 3427;
	vector<double> RMS_total(11264,0.0f); //Stores the Coherent noise levels for entire TPC
	vector<int> Entries(11264,0.0f);
	vector<vector<double>> FFT_total(11264,vector<double>(event_len/2+2,0));


	//Grabs the histograms and merges the wire info into a 2D vector for all the wire info of an event
	TIter next(infile->GetListOfKeys());
    TKey* key;
    for (int e = 1; e <= TOTAL_EVT; ++e){
	    vector<vector<double>> TPC_wires(11264,vector<double>(event_len,0));
	    while ((key = (TKey*)next())) {
	        // Check if the object is a 2D histogram
	        if (TH2* hist2D = dynamic_cast<TH2*>(key->ReadObj())) {
	            std::cout << "2D Histogram: " << hist2D->GetName() << std::endl;
	        `	string hist_name =hist2D->GetName();
	        	event = stoi(hist_name.substr(9));
	        	if (e != event) continue;
	            if (hist_name[3] != 'r') continue;
	            // Access 2D histogram data (e.g., print bin contents)
	            int nBinsX = hist2D->GetNbinsX();
	            int nBinsY = hist2D->GetNbinsY();

	            vector<vector<double>> u0_wires();
	            vector<vector<double>> v0_wires();
	            vector<vector<double>> w0_wires();
	            vector<vector<double>> u1_wires();
	            vector<vector<double>> v1_wires();
	            vector<vector<double>> w1_wires();

	            if (hist_name[1] == 'u' && hist_name[7] == '0'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "U0: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	u0_wires.push_back(wire);
	            	}
	            }
	            else if (hist_name[1] == 'v' && hist_name[7] == '0'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "V0: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	v0_wires.push_back(wire);
	            	}
	            }
	            else if (hist_name[1] == 'w' && hist_name[7] == '0'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "W0: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	w0_wires.push_back(wire);
	            	}
	            }
	            else if (hist_name[1] == 'u' && hist_name[7] == '1'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "U1: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	u1_wires.push_back(wire);
	            	}
	            }
	            else if (hist_name[1] == 'v' && hist_name[7] == '1'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "V1: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	v1_wires.push_back(wire);
	            	}
	            }
	            else if (hist_name[1] == 'w' && hist_name[7] == '1'){
	            	for (int i = 1; i <= nBinsX; ++i) {
	            		vector<float> wire();
	                	for (int j = 1; j <= nBinsY; ++j) {
	                    	double binContent = hist2D->GetBinContent(i, j);
	                    	wire.push_back(j);
	                    	std::cout << "W1: (" << i << ", " << j << "): " << binContent << std::endl;
	                	}
	                	w1_wires.push_back(wire);
	            	}
	            }   
	    	}
    	}
    	TPC_wires.insert(TPC_wires.end(), u0_wires.begin(), u0_wires.end());
    	TPC_wires.insert(TPC_wires.end(), v0_wires.begin(), v0_wires.end());
    	TPC_wires.insert(TPC_wires.end(), w0_wires.begin(), w0_wires.end());
    	TPC_wires.insert(TPC_wires.end(), u1_wires.begin(), u1_wires.end());
    	TPC_wires.insert(TPC_wires.end(), v1_wires.begin(), v1_wires.end());
    	TPC_wires.insert(TPC_wires.end(), w1_wires.begin(), w1_wires.end());
    
	
	
		cout<<"Running Events"<<endl;
		int evt = 0;

        //Goes over all of the channels and does the analysis
		for(int ki=0; ki<11264;ki++){
			cout<<" Channel: "<<ki<<endl;
			cout<<"Channel size: "<<myADC[index].NADC()<<endl;

			//If channel is responsive the channel will grab the noise 
			bool skip_channel = false;
			vector<short> x(myADC[index].Samples(),0);
			vector<double> y(myADC[index].Samples(),0);
			for (size_t itick=0; itick < TPC_wires[ki].size(); ++itick){ 
				float pedestal = median(TPC_wires[ki])
				if (abs(TPC_wires[ki][itick]-pedestal) >  20){
					skip_channel = true;
					break;
				}
				x[itick] = TPC_wires[ki][itick];//-myADC[index].GetPedestal();//
				y[itick] = TPC_wires[ki][itick];//-myADC[index].GetPedestal();
;//
			}
			vector<double> raw_channel_fft = FFT(y);
			transform(Raw_FFT_total[ki].begin(),Raw_FFT_total[ki].end(),raw_channel_fft.begin(),Raw_FFT_total[ki].begin(),plus<double>());
			double RMS = Noise_levels(x);
			Raw_RMS_total[ki] = Raw_RMS_total.at(ki)+RMS;
			entries[ki] = entries.at(ki)+1;


		}

		cout<<"Event:"<<e<<endl;
		//break;
	}
	
	TFile* file = new TFile("noise_output_fft.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float raw_rms;
	int entries;
	float int_rms;
	int int_entries;
	float int_FFT;
	float coh_FFT;
	float raw_FFT;

	//vector<float> avg_FFT;
	//tree->Branch("coh_rms", &avg_rms, "avg_rms/F");
	tree->Branch("entries", &entries, "entries/I");
	tree->Branch("raw_rms", &raw_rms, "raw_rms/F");
	tree->Branch("raw_FFT", &raw_FFT, "raw_FFT/F");
	tree->SetBranchStatus("raw_rms", 1);
    tree->SetBranchStatus("entries", 0);
    tree->SetBranchStatus("raw_FFT", 0);
	for(int ch = 0; ch<RMS_total.size(); ch++){
		raw_rms = RMS_total.at(ch)/Entries.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("raw_rms", 0);
    tree->SetBranchStatus("entries", 1);
	for(int ch = 0; ch<Entries.size(); ch++){
		entries = Entries.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("raw_FFT", 1);
    tree->SetBranchStatus("entries", 0);
	for(int ch = 0; ch<Raw_FFT_total.size(); ch++){
		for (size_t c = 0; c < Raw_FFT_total[ch].size(); ++c) {
			raw_FFT = Raw_FFT_total[ch][c];
			tree->Fill();
        }
    }
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;


}

//void TPC_coherent_noise(TString inputFile="/exp/sbnd/data/users/dcarber/prodgenie_cosmic_rockbox_sbnd_GenieGen-20241015T135353_G4-20241015T140624_DetSim-20241017T143903.root")
void TPC_WC_noise_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run17470/sbnd-data-check.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}

