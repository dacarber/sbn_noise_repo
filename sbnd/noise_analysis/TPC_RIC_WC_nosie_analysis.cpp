
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
#include <sys/resource.h>
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
vector<double> Int_removal(vector<double> noise, vector<double> int_noise){
	vector<double> coh_waveform(noise.size(),0);
	float intrinsic;
	float raw;
	int sum = accumulate(noise.begin(),noise.end(),0);
	cout<<"Check size"<<noise.size()<<endl;
	if (noise.size() != 3427 || sum == 0 ){
		return coh_waveform;
	}
	transform(noise.begin(), noise.end(), noise.begin(), [](double x) { return x < 0 ? - (x * x) : (x * x);  });
	transform(int_noise.begin(), int_noise.end(), int_noise.begin(), [](double x) { return x < 0 ? - (x * x) : (x * x);  });
	
	transform(noise.begin(),noise.end(),int_noise.begin(),coh_waveform.begin(),minus<float>());
	cout<<"Returning vector"<<coh_waveform[100]<<endl;

	return coh_waveform;
}
vector<double> FFT(vector<double> noise_channel){
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
	return fftMag;
}
float Median(vector<double> &vec) {
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
void getMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    long memoryUsage = usage.ru_maxrss; // in kilobytes
    std::cout << "Memory usage: " << memoryUsage << " KB" << std::endl;
}
void LoadRawDigits(TFile *inFile)
{	
	int event_len = 3427;
	vector<double> RMS_orig_total(12276,0.0f); //Stores the Full noise levels for entire TPC
	vector<int> Entries_orig(12276,0.0f);
	vector<vector<double>> FFT_orig_total(12276,vector<double>(event_len/2+3,0));
	vector<double> RMS_raw_total(12276,0.0f); //Stores the NF noise levels for entire TPC
	vector<int> Entries_raw(12276,0.0f);
	vector<vector<double>> FFT_raw_total(12276,vector<double>(event_len/2+3,0));

	//Work in progress
	vector<double> RMS_coh_total(12276,0.0f); //Stores the Coherent noise levels for entire TPC
	vector<int> Entries_coh(12276,0.0f);
	vector<vector<double>> FFT_coh_total(12276,vector<double>(event_len/2+3,0));

	int channel_base = 0;
	int event = 0;
	//Grabs the histograms and merges the wire info into a 2D vector for all the wire info of an event
    TIter next(inFile->GetListOfKeys());
    int numKeys = inFile->GetNkeys();
	TKey* key;
    int i = 0;
    while ((key = (TKey*)next())) {
        // Check if the object is a 2D histogram
        getMemoryUsage();
        std::cout << "Key Name: " << key->GetName() << std::endl;
        TObject* obj = key->ReadObj();
        TH2D* hist2D = (TH2D*)key->ReadObj();
            string hist_name = hist2D->GetName();
        	std::cout <<hist_name[1]<<hist_name[3] <<hist_name.back() << std::endl;
        	if (hist_name[3] == 'r'){
        		delete key;
        		delete hist2D;
        		delete obj;
        		continue;
        	}
        	int cycle = key->GetCycle();  // Get cycle number
            string name = "h"+ string(1,hist_name[1]) +"_raw"+ string(1,hist_name.back()) +";"+to_string(cycle);
            TH2D* hist2D_raw = (TH2D*)inFile->Get(name.c_str());
    		std::cout << "Found Histogram: " << hist_name
            << ", Entries: " << hist2D->GetNbinsX() << std::endl;
            std::cout << "Key Name: " << key->GetName() << std::endl;
            int nBinsX = hist2D->GetNbinsX();
        	int nBinsY = hist2D->GetNbinsY();
        	int nBinsX_raw = hist2D_raw->GetNbinsX();
        	int nBinsY_raw = hist2D_raw->GetNbinsY();
        	event+=1;
        	
        	vector<double> wire_ped;
        	vector<double> wire;
        	vector<double> wire_ped_raw;
        	vector<double> wire_raw;


        	for (int x = 0; x <= nBinsX; ++x) {
        		wire_ped.clear();
        		wire.clear();
        		for (int y = 0; y <= nBinsY; ++y) {
            		double binContent = hist2D->GetBinContent(x, y);
            		if (binContent == 0){continue;};
            		wire_ped.push_back(binContent);
            		wire.push_back(binContent);
            	}
            for (int x = 0; x <= nBinsX_raw; ++x) {
        		wire_ped_raw.clear();
        		wire_raw.clear();
        		for (int y = 0; y <= nBinsY_raw; ++y) {
            		double binContent = hist2D_raw->GetBinContent(x, y);
            		if (binContent == 0){continue;};
            		wire_ped_raw.push_back(binContent);
            		wire_raw.push_back(binContent);
            	}


            	//Cleaning up the channels
            	if (wire.size()!=3427 or wire_raw.size()!=3427){continue;};
            	int pedestal = Median(wire_ped);
    			transform(wire.begin(), wire.end(), wire.begin(),[pedestal](double elem) { return elem - pedestal; });
    			auto max_el = max_element(wire.begin(), wire.end(), [](double a, double b) {return std::abs(a) < std::abs(b);});
    			if (abs(max_el[0]) > 20.0){
    				continue;
    			}
    			if (accumulate(wire.begin(), wire.end(), 0) == 0){
    				continue;
    			}


    			if (hist_name[3] == 'o'){
            		if (hist_name[1] == 'u' and hist_name.back() == '0'){
            			channel_base = 0;
            		}
            		else if (hist_name[1] == 'v' and hist_name.back() == '0'){
            			channel_base = 1984;
            		}
            		else if (hist_name[1] == 'w' and hist_name.back() == '0'){
            			channel_base = 1984*2;
            		}
            		else if (hist_name[1] == 'u' and hist_name.back() == '1'){
            			channel_base = 1984*2+1670;
            		}
            		else if (hist_name[1] == 'v' and hist_name.back() == '1'){
            			channel_base = 1984*3+1670;
            		}
            		else if (hist_name[1] == 'w' and hist_name.back() == '1'){
            			channel_base = 1984*4+1670;
            		}
            		int channel = x+channel_base-1;

        			vector<double> orig_channel_fft = FFT(wire);
					transform(FFT_orig_total[channel].begin(),FFT_orig_total[channel].end(),orig_channel_fft.begin(),FFT_orig_total[channel].begin(),plus<double>());
					double RMS_orig = Noise_levels(wire);
					RMS_orig_total[channel] = RMS_orig_total.at(channel)+RMS_orig;
					Entries_orig[channel] = Entries_orig.at(channel)+1;

					vector<double> raw_channel_fft = FFT(wire_raw);
					transform(FFT_raw_total[channel].begin(),FFT_raw_total[channel].end(),raw_channel_fft.begin(),FFT_raw_total[channel].begin(),plus<double>());
					double RMS_raw = Noise_levels(wire_raw);
					RMS_raw_total[channel] = RMS_raw_total.at(channel)+RMS_raw;
					Entries_raw[channel] = Entries_raw.at(channel)+1;

					vector<double> coh_wave = Int_removal(wire,wire_raw);
					vector<double> coh_channel_fft = FFT(coh_wave);
					transform(FFT_coh_total[channel].begin(),FFT_coh_total[channel].end(),coh_channel_fft.begin(),FFT_coh_total[channel].begin(),plus<double>());
					double RMS_coh = Noise_levels(coh_wave);
					RMS_coh_total[channel] = RMS_coh_total.at(channel)+RMS_coh;
					Entries_coh[channel] = Entries_coh.at(channel)+1;

            	}
        	}	
    	delete key;
        delete hist2D;
        delete obj;
        delete hist2D_raw;
	}
	
	TFile* file = new TFile("noise_output_fft.root", "RECREATE");
	TTree* tree = new TTree("tpc_noise", "tpc_noise");
	float raw_rms;
	int raw_entries;
	float orig_rms;
	int orig_entries;
	float orig_FFT;
	float raw_FFT;
	float coh_rms;
	int coh_entries;
	float coh_FFT;

	//vector<float> avg_FFT;
	//tree->Branch("coh_rms", &avg_rms, "avg_rms/F");
	tree->Branch("raw_entries", &raw_entries, "raw_entries/I");
	tree->Branch("raw_rms", &raw_rms, "raw_rms/F");
	tree->Branch("raw_FFT", &raw_FFT, "raw_FFT/F");

	tree->Branch("orig_entries", &orig_entries, "orig_entries/I");
	tree->Branch("orig_rms", &orig_rms, "orig_rms/F");
	tree->Branch("orig_FFT", &orig_FFT, "orig_FFT/F");

	tree->Branch("coh_entries", &coh_entries, "coh_entries/I");
	tree->Branch("coh_rms", &coh_rms, "coh_rms/F");
	tree->Branch("coh_FFT", &coh_FFT, "coh_FFT/F");

	tree->SetBranchStatus("raw_rms", 1);
    tree->SetBranchStatus("raw_entries", 0);
    tree->SetBranchStatus("raw_FFT", 0);
    tree->SetBranchStatus("orig_rms", 0);
    tree->SetBranchStatus("orig_entries", 0);
    tree->SetBranchStatus("orig_FFT", 0);
    tree->SetBranchStatus("coh_rms", 0);
    tree->SetBranchStatus("coh_entries", 0);
    tree->SetBranchStatus("coh_FFT", 0);


	for(int ch = 0; ch<RMS_raw_total.size(); ch++){
		raw_rms = RMS_raw_total.at(ch)/Entries_raw.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("raw_rms", 0);
    tree->SetBranchStatus("raw_entries", 1);
	for(int ch = 0; ch<Entries_raw.size(); ch++){
		raw_entries = Entries_raw.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("raw_FFT", 1);
    tree->SetBranchStatus("raw_entries", 0);
	for(int ch = 0; ch<FFT_raw_total.size(); ch++){
		for (size_t c = 0; c < FFT_raw_total[ch].size(); ++c) {
			raw_FFT = FFT_raw_total[ch][c];
			tree->Fill();
        }
    }
    tree->SetBranchStatus("orig_rms", 1);
    tree->SetBranchStatus("raw_FFT", 0);
    for(int ch = 0; ch<RMS_orig_total.size(); ch++){
		orig_rms = RMS_orig_total.at(ch)/Entries_orig.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("orig_rms", 0);
    tree->SetBranchStatus("orig_entries", 1);
	for(int ch = 0; ch<Entries_orig.size(); ch++){
		orig_entries = Entries_orig.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("orig_FFT", 1);
    tree->SetBranchStatus("orig_entries", 0);
	for(int ch = 0; ch<FFT_orig_total.size(); ch++){
		for (size_t c = 0; c < FFT_orig_total[ch].size(); ++c) {
			orig_FFT = FFT_orig_total[ch][c];
			tree->Fill();
        }
    }
    tree->SetBranchStatus("coh_rms", 1);
    tree->SetBranchStatus("orig_FFT", 0);
    for(int ch = 0; ch<RMS_coh_total.size(); ch++){
		coh_rms = RMS_coh_total.at(ch)/Entries_coh.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("coh_rms", 0);
    tree->SetBranchStatus("coh_entries", 1);
	for(int ch = 0; ch<Entries_coh.size(); ch++){
		coh_entries = Entries_coh.at(ch);
		tree->Fill();	
	}
	tree->SetBranchStatus("coh_FFT", 1);
    tree->SetBranchStatus("coh_entries", 0);
	for(int ch = 0; ch<FFT_coh_total.size(); ch++){
		for (size_t c = 0; c < FFT_coh_total[ch].size(); ++c) {
			coh_FFT = FFT_coh_total[ch][c];
			tree->Fill();
        }
    }
	file->Write();
	file->Close();
	
	cout<<"Got ADC and Pedestal"<<endl;


}

//void TPC_coherent_noise(TString inputFile="/exp/sbnd/data/users/dcarber/prodgenie_cosmic_rockbox_sbnd_GenieGen-20241015T135353_G4-20241015T140624_DetSim-20241017T143903.root")
void TPC_RIC_WC_noise_analysis(TString inputFile="/exp/sbnd/data/users/dcarber/tpcnoise/run17470/sbnd-data-check.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}

