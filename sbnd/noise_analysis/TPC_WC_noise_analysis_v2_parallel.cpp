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
#include <omp.h> // The multithreading library
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

void LoadRawDigits(TFile *inFile) {
    ROOT::EnableThreadSafety(); // Important for ROOT 6+
    int event_len = 3427;

    // --- 1. DEFINE GLOBAL ACCUMULATORS (Move these to the top) ---
    // These hold the final sum of everything.
    vector<double> global_RMS_orig(12276, 0.0);
    vector<int>    global_Entries_orig(12276, 0);
    vector<vector<double>> global_FFT_orig(12276, vector<double>(event_len/2+3, 0));

    vector<double> global_RMS_raw(12276, 0.0);
    vector<int>    global_Entries_raw(12276, 0);
    vector<vector<double>> global_FFT_raw(12276, vector<double>(event_len/2+3, 0));

    // --- 2. CONVERT TLIST TO VECTOR (New Logic) ---
    // OpenMP can't handle "while(next())". It needs a simple vector to loop over.
    TList* keyList = inFile->GetListOfKeys();
    int nKeys = keyList->GetSize();
    std::vector<TKey*> all_keys;
    all_keys.reserve(nKeys);
    
    TIter next(keyList);
    TKey* k;
    while ((k = (TKey*)next())) { all_keys.push_back(k); }
    
    cout << "Processing " << nKeys << " events in parallel..." << endl;

    #pragma omp parallel 
    {
        // --- 4. SETUP PRIVATE VARIABLES (Per Core) ---
        // Every core creates its OWN temporary vectors here.
        // If we didn't do this, cores would fight over the same memory and crash.
        
        // Math resources (FFT plan)
        Int_t n_size = event_len;
        TVirtualFFT* fft_planner = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
        double* fft_in = new double[event_len];
        vector<double> fft_out(event_len/2 + 2);
        
        // Data vectors
        vector<double> wire_ped; wire_ped.reserve(event_len);
        vector<double> wire;     wire.reserve(event_len);

        // Local Accumulators (Results for just this core)
        vector<double> loc_RMS_orig(12276, 0.0);
        vector<int>    loc_Entries_orig(12276, 0);
        vector<vector<double>> loc_FFT_orig(12276, vector<double>(event_len/2+3, 0));
        
        vector<double> loc_RMS_raw(12276, 0.0);
        vector<int>    loc_Entries_raw(12276, 0);
        vector<vector<double>> loc_FFT_raw(12276, vector<double>(event_len/2+3, 0));

        // --- 5. THE PARALLEL LOOP ---
        // This splits the vector of keys: Core 1 gets keys 0-100, Core 2 gets 101-200, etc.
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < nKeys; ++i) {
            
            TKey* key = all_keys[i];
            TH2D* hist2D = nullptr;

            // --- 6. CRITICAL READ (Safety Lock) ---
            // Only one thread can touch the hard drive/TFile at a time.
            #pragma omp critical (FileIO) 
            {
                TObject* obj = key->ReadObj();
                if (obj->InheritsFrom("TH2D")) hist2D = (TH2D*)obj;
                else delete obj;
            } 
            
            if (!hist2D) continue; 

            // --- 7. YOUR ANALYSIS LOGIC ---
            // (Paste your loop over bins, FFT calls, and Median logic here)
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
                // 1. Run FFT
                FFT(wire, fft_output_buffer, fft_planner, fft_input_buffer);
                
                // 2. Add to LOCAL Accumulator
                for(size_t k=0; k < fft_output_buffer.size(); k++) {
                     if (k < loc_FFT_orig[channel].size()) {
                        loc_FFT_orig[channel][k] += fft_output_buffer[k]; 
                     }
                }
                
                // 3. Add RMS/Entries to LOCAL Accumulators
                loc_RMS_orig[channel] += Noise_levels(wire); 
                loc_Entries_orig[channel]++;                 
            } 
            else if (hist_name[3] == 'r') { // "raw"
                FFT(wire, fft_output_buffer, fft_planner, fft_input_buffer);

                for(size_t k=0; k < fft_output_buffer.size(); k++) {
                     if (k < loc_FFT_raw[channel].size()) {
                        loc_FFT_raw[channel][k] += fft_output_buffer[k]; 
                     }
                }

                loc_RMS_raw[channel] += Noise_levels(wire); 
                loc_Entries_raw[channel]++;                 
            }
        }
            // IMPORTANT: Save results to 'loc_RMS_orig' (local), NOT 'global_RMS_orig'.
            
            delete hist2D; // Clean up immediately
        } 

        // --- 8. CLEANUP LOCAL RESOURCES ---
        delete[] fft_in;
        delete fft_planner;

        // --- 9. CRITICAL MERGE (Combine Results) ---
        // The cores wait in line here to add their local bucket to the global bucket.
        #pragma omp critical (MergeData)
        {
            for(int c=0; c<12276; ++c) {
                global_RMS_orig[c]     += loc_RMS_orig[c];
                global_Entries_orig[c] += loc_Entries_orig[c];
                global_RMS_raw[c]     += loc_RMS_raw[c];
                global_Entries_raw[c] += loc_Entries_raw[c];
                // ... same for raw ...
                
                // Merge FFT vectors
                for(size_t k=0; k<loc_FFT_orig[c].size(); k++) {
                    global_FFT_orig[c][k] += loc_FFT_orig[c][k];
                }
                for(size_t k=0; k<loc_FFT_raw[c].size(); k++) {
                    global_FFT_raw[c][k] += loc_FFT_raw[c][k];
                }
                // ... same for raw FFT ...
            }
        }
    } // <--- End of Parallel Region
    // =======================================================
    // STEP 4: WRITE OUTPUT TO TTREE
    // =======================================================
    cout << "Writing output to noise_output_fft.root..." << endl;

    TFile* file = new TFile("noise_output_fft.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "tpc_noise");

    // --- Branch Variables ---
    // These hold the values for the "Current Row" being written
    int ch_id;
    
    // Raw Data Variables
    float raw_rms_val;
    int raw_entries_val;
    vector<float> raw_fft_vec; 
    
    // Orig Data Variables
    float orig_rms_val;
    int orig_entries_val;
    vector<float> orig_fft_vec; 

    // --- Define Branches ---
    tree->Branch("channel",      &ch_id,            "channel/I");
    
    tree->Branch("raw_entries",  &raw_entries_val,  "raw_entries/I");
    tree->Branch("raw_rms",      &raw_rms_val,      "raw_rms/F");
    tree->Branch("raw_fft",      &raw_fft_vec);     // Vector branch

    tree->Branch("orig_entries", &orig_entries_val, "orig_entries/I");
    tree->Branch("orig_rms",     &orig_rms_val,     "orig_rms/F");
    tree->Branch("orig_fft",     &orig_fft_vec);    // Vector branch

    // --- Loop Over All Channels and Fill ---
    // The previous code had 12276 channels total
    for(int ch = 0; ch < 12276; ch++){
        ch_id = ch;

        // 1. Process RAW Data
        if (global_Entries_raw[ch] > 0) {
            raw_entries_val = global_Entries_raw[ch];
            // Calculate Average RMS
            raw_rms_val = global_RMS_raw[ch] / raw_entries_val;
            
            // Calculate Average FFT
            // We copy the Double vector to Float vector and divide by N
            raw_fft_vec.clear();
            raw_fft_vec.reserve(global_FFT_raw[ch].size());
            for(double val : global_FFT_raw[ch]) {
                raw_fft_vec.push_back( (float)(val / raw_entries_val) );
            }
        } else {
            // If no data for this channel, fill with zeros/empty
            raw_entries_val = 0;
            raw_rms_val = 0.0;
            raw_fft_vec.clear();
        }

        // 2. Process ORIG Data
        if (global_Entries_orig[ch] > 0) {
            orig_entries_val = global_Entries_orig[ch];
            // Calculate Average RMS
            orig_rms_val = global_RMS_orig[ch] / orig_entries_val;
            
            // Calculate Average FFT
            orig_fft_vec.clear();
            orig_fft_vec.reserve(global_FFT_orig[ch].size());
            for(double val : global_FFT_orig[ch]) {
                orig_fft_vec.push_back( (float)(val / orig_entries_val) );
            }
        } else {
            orig_entries_val = 0;
            orig_rms_val = 0.0;
            orig_fft_vec.clear();
        }

        // Fill the row for this channel
        tree->Fill(); 
    }

    // Save and Close
    file->Write();
    file->Close();
    
    cout << "Done. File closed." << endl;
}
    
    



void TPC_WC_noise_analysis_v2_parallel(TString inputFile="sbnd-data-check.root")
{	
	cout<<"Get ready for the rollercoaster of me learning Root and C++"<<endl;
	
	TFile *inFile = TFile::Open(inputFile.Data());
	cout<<"Got File"<<endl;
	LoadRawDigits(inFile);
}
