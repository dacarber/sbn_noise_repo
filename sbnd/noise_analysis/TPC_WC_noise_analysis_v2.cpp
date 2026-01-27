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
#include <sys/resource.h>
//#include <bits/stdc++.h> 

using namespace std;
double Noise_levels(const vector<double>& noise_channels){
    if (noise_channels.empty()) return 0.0;
    double sum = 0.0;
    double sum_sq = 0.0;
    for (double val : noise_channels) {
        sum += val;
        sum_sq += val * val;
    }
    double n = (double)noise_channels.size();
    double mean = sum / n;
    double variance = (sum_sq / n) - (mean * mean);
    return sqrt(abs(variance));     
}

float Median(vector<double> vec) { 
    size_t n = vec.size();
    if (n == 0) return 0.0;
    
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    double median = *target;

    // Handle even-sized vectors (average of two middle elements)
    if (n % 2 == 0) {
        auto target_neighbor = std::max_element(vec.begin(), target);
        median = (*target + *target_neighbor) / 2.0;
    }
    return median;
}

// Change arguments to accept the pre-made planner
void FFT(const vector<double>& noise_channel, vector<double>& output_mag, TVirtualFFT* fft_planner, double* input_buffer) {
    int vec_size = noise_channel.size();
    
    // Copy data into the pre-allocated buffer (No 'new' memory allocation!)
    for (size_t i = 0; i < vec_size; ++i) {
        input_buffer[i] = noise_channel[i];
    }

    // Use the existing planner
    fft_planner->SetPoints(input_buffer);
    fft_planner->Transform();

    double fftReal, fftImag;
    
    // Fill the output vector (passed by reference to avoid copying return values)
    // We assume output_mag is already size (vec_size/2 + 2)
    for(size_t k=1; k < vec_size / 2 + 2; k++){
        fft_planner->GetPointComplex(k, fftReal, fftImag);
        if (k == vec_size / 2 + 1){
            output_mag[k] = 1.0;
        } else {
            output_mag[k] = TMath::Sqrt(fftReal*fftReal + fftImag*fftImag);
        }
    }
}

void LoadRawDigits(TFile *inFile)
{   
    // 1. MEMORY OPTIMIZATION: Reserve vectors outside the loop
    // This prevents creating/destroying millions of vectors
    int event_len = 3427; // Constant size
    
    // --- OPTIMIZATION 1: PRE-ALLOCATE FFT RESOURCES ---
    // Create the planner ONCE.
    Int_t n_size = event_len;
    TVirtualFFT* fft_planner = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
    // Pre-allocate the buffer arrays used by FFT
    double* fft_input_buffer = new double[event_len];
    vector<double> fft_output_buffer(event_len/2 + 2);

    // Standard pre-allocations
    vector<double> wire_ped; 
    wire_ped.reserve(event_len);
    vector<double> wire; 
    wire.reserve(event_len);

    // Data containers
    vector<double> RMS_orig_total(12276, 0.0);
    vector<int> Entries_orig(12276, 0);
    vector<vector<double>> FFT_orig_total(12276, vector<double>(event_len/2+3, 0));

    vector<double> RMS_raw_total(12276, 0.0);
    vector<int> Entries_raw(12276, 0);
    vector<vector<double>> FFT_raw_total(12276, vector<double>(event_len/2+3, 0));

    int channel_base = 0;

    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int total_keys = inFile->GetListOfKeys()->GetSize();
	int counter = 0;
	TStopwatch timer;
	timer.Start();
    // --- READING LOOP ---
    while ((key = (TKey*)next())) {
        // Fix: Check inheritance first, read once
        TObject* obj = key->ReadObj();
        //if (!obj->InheritsFrom("TH2D")) {
        //    delete obj;
        //    continue; 
        //}
        TH2D* hist2D = (TH2D*)obj;
        
        string hist_name = hist2D->GetName();
        // std::cout << "Processing: " << hist_name << std::endl;

        int nBinsX = hist2D->GetNbinsX();
        int nBinsY = hist2D->GetNbinsY();
        
        // Loop over wires in this histogram
        for (int x = 0; x <= nBinsX; ++x) {
            // OPTIMIZATION: Clear instead of re-declaring
            wire_ped.clear();
            wire.clear();
            
            // Extract waveform
            for (int y = 0; y <= nBinsY; ++y) {
                double binContent = hist2D->GetBinContent(x, y);
                if (binContent == 0) continue;
                wire_ped.push_back(binContent);
                wire.push_back(binContent);
            }
            
            if (wire.size() != event_len) continue;

            // Calculate Pedestal
            // Note: Median takes 'wire_ped' by copy, so it's safe
            int pedestal = Median(wire_ped); 
            
            // Subtract Pedestal
            transform(wire.begin(), wire.end(), wire.begin(), 
                      [pedestal](double elem) { return elem - pedestal; });

            // Check for Signal (Pulse detection)
            // Using a lambda to find max absolute value
            double max_val = 0;
            for(double w : wire) if(std::abs(w) > max_val) max_val = std::abs(w);
            
            // Note: "median" variable was undefined in your code, assuming 0 or pedestal logic?
            // I'll assume you meant comparing max_val to a threshold
            if (max_val > 20.0) {
                 // cout << "Wire " << x << " has signal" << endl;
                 continue;
            }
            
            // Determine Channel Base
            // (Condensed logic for readability)
            char plane = hist_name[1]; // 'u', 'v', 'w'
            char tpc = hist_name.back(); // '0', '1'
            
            if (tpc == '0') {
                if (plane == 'u') channel_base = 0;
                else if (plane == 'v') channel_base = 1984;
                else if (plane == 'w') channel_base = 1984*2;
            } else if (tpc == '1') {
                if (plane == 'u') channel_base = 1984*2 + 1670;
                else if (plane == 'v') channel_base = 1984*3 + 1670;
                else if (plane == 'w') channel_base = 1984*4 + 1670;
            }

            int channel = x + channel_base - 1;
            if (channel < 0 || channel >= 12276) continue; // Safety check

            // Process based on histogram type (orig vs raw)
            if (hist_name[3] == 'o') { // "orig"
                //vector<double> orig_channel_fft = FFT(wire);
                FFT(wire, fft_output_buffer, fft_planner, fft_input_buffer);
                // Add FFT to total
                for(size_t k=0; k < fft_output_buffer.size(); k++) {
                     if (k < FFT_orig_total[channel].size())
                        FFT_orig_total[channel][k] += fft_output_buffer[k];
                }
                //for(size_t k=0; k<orig_channel_fft.size(); k++) {
                //    if (k < FFT_orig_total[channel].size())
                //        FFT_orig_total[channel][k] += orig_channel_fft[k];
                //}
                
                RMS_orig_total[channel] += Noise_levels(wire);
                Entries_orig[channel]++;
            } 
            else if (hist_name[3] == 'r') { // "raw"
                //vector<double> raw_channel_fft = FFT(wire);
                FFT(wire, fft_output_buffer, fft_planner, fft_input_buffer);
                for(size_t k=0; k < fft_output_buffer.size(); k++) {
                     if (k < FFT_raw_total[channel].size())
                        FFT_raw_total[channel][k] += fft_output_buffer[k];
                }
                //for(size_t k=0; k<raw_channel_fft.size(); k++) {
                //     if (k < FFT_raw_total[channel].size())
                //        FFT_raw_total[channel][k] += raw_channel_fft[k];
                //}

                RMS_raw_total[channel] += Noise_levels(wire);
                Entries_raw[channel]++;
            }
        }
         
        delete hist2D; // Only delete the object we casted
        counter++;
        if (counter % 10 == 0 || counter == total_keys) { 
        Double_t time_elapsed = timer.RealTime();
        timer.Continue(); // Resume the timer
        
        Double_t rate = counter / time_elapsed; // keys per second
        Double_t time_remaining = (total_keys - counter) / rate;
        
        std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% "
                  << "| Rate: " << std::fixed << std::setprecision(1) << rate << " keys/s "
                  << "| ETA: " << int(time_remaining) << "s   " << std::flush;
    }
    // ---------------------------------------

std::cout << std::endl; // Move to next line when done
    }
    delete[] fft_input_buffer;
	delete fft_planner; 
    // --- SAVING LOOP (Completely Rewritten) ---
    cout << "Writing to TTree..." << endl;
    
    TFile* file = new TFile("noise_output_fft.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "tpc_noise");

    // Branch Variables
    int ch_id;
    float raw_rms_val, orig_rms_val;
    int raw_entries_val, orig_entries_val;
    vector<float> raw_fft_vec;  // Use vectors for array data
    vector<float> orig_fft_vec; 

    tree->Branch("channel", &ch_id, "channel/I");
    tree->Branch("raw_entries", &raw_entries_val, "raw_entries/I");
    tree->Branch("raw_rms", &raw_rms_val, "raw_rms/F");
    tree->Branch("raw_fft", &raw_fft_vec); // Automatic std::vector branch

    tree->Branch("orig_entries", &orig_entries_val, "orig_entries/I");
    tree->Branch("orig_rms", &orig_rms_val, "orig_rms/F");
    tree->Branch("orig_fft", &orig_fft_vec);

    // Single loop over all channels to fill the tree
    for(int ch = 0; ch < 12276; ch++){
        ch_id = ch;

        // Process RAW
        if (Entries_raw[ch] > 0) {
            raw_entries_val = Entries_raw[ch];
            raw_rms_val = RMS_raw_total[ch] / raw_entries_val;
            
            // Average the FFT
            raw_fft_vec.clear();
            for(double val : FFT_raw_total[ch]) {
                raw_fft_vec.push_back( val / raw_entries_val );
            }
        } else {
            raw_entries_val = 0;
            raw_rms_val = 0;
            raw_fft_vec.clear();
        }

        // Process ORIG
        if (Entries_orig[ch] > 0) {
            orig_entries_val = Entries_orig[ch];
            orig_rms_val = RMS_orig_total[ch] / orig_entries_val;
            
            // Average the FFT
            orig_fft_vec.clear();
            for(double val : FFT_orig_total[ch]) {
                orig_fft_vec.push_back( val / orig_entries_val );
            }
        } else {
            orig_entries_val = 0;
            orig_rms_val = 0;
            orig_fft_vec.clear();
        }

        tree->Fill(); // Saves one row containing ALL data for this channel
    }

    file->Write();
    file->Close();
    cout << "Done." << endl;
}



void TPC_WC_noise_analysis_v2(TString inputFile="sbnd-data-check.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}
