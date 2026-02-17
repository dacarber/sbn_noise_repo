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
#include <TH2D.h>
#include <TVirtualFFT.h>
#include <fstream>
#include <TChain.h>
#include <TKey.h>
#include <TStopwatch.h>
#include <sys/resource.h>

using namespace std;

// --- Helper Functions ---

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

    if (n % 2 == 0) {
        auto target_neighbor = std::max_element(vec.begin(), target);
        median = (*target + *target_neighbor) / 2.0;
    }
    return median;
}

// FFT Function (Kept for reference, but not called in loop to save speed/space)
void FFT(const vector<double>& noise_channel, vector<double>& output_mag, TVirtualFFT* fft_planner, double* input_buffer) {
    int vec_size = noise_channel.size();
    for (size_t i = 0; i < vec_size; ++i) {
        input_buffer[i] = noise_channel[i];
    }
    fft_planner->SetPoints(input_buffer);
    fft_planner->Transform();
    double fftReal, fftImag;
    for(size_t k=1; k < vec_size / 2 + 2; k++){
        fft_planner->GetPointComplex(k, fftReal, fftImag);
        if (k == vec_size / 2 + 1){
            output_mag[k] = 1.0;
        } else {
            output_mag[k] = TMath::Sqrt(fftReal*fftReal + fftImag*fftImag);
        }
    }
}

// --- Main Loader ---

void LoadRawDigits(TFile *inFile)
{   
    // 1. SETUP OUTPUT TREE
    // We create the file and tree BEFORE looping, so we can fill it event-by-event.
    cout << "Setting up output TTree..." << endl;
    TFile* outFile = new TFile("noise_rms_per_event.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "RMS per wire per event");

    // Branch Variables
    std::string b_hist_name;           // Name of the histogram (contains info on plane/raw/orig)
    int b_event_idx;                   // Counter index
    std::vector<int> b_channels;       // Vector of Channel IDs for this event
    std::vector<float> b_rms;          // Vector of RMS values for this event

    // Create Branches
    tree->Branch("hist_name", &b_hist_name);
    tree->Branch("event_idx", &b_event_idx, "event_idx/I");
    tree->Branch("channels", &b_channels);
    tree->Branch("rms", &b_rms);

    // 2. PRE-ALLOCATION
    int event_len = 3427; 
    
    // Vectors for processing waveforms
    vector<double> wire_ped; 
    wire_ped.reserve(event_len);
    vector<double> wire; 
    wire.reserve(event_len);

    // FFT Setup (Initialized but not used to save space/time unless requested)
    //Int_t n_size = event_len;
    //TVirtualFFT* fft_planner = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
    //double* fft_input_buffer = new double[event_len];
    //vector<double> fft_output_buffer(event_len/2 + 2);

    // Iteration Setup
    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    int counter = 0;
    
    TStopwatch timer;
    timer.Start();

    // 3. MAIN EVENT LOOP
    while ((key = (TKey*)next())) {
        
        // Read Object
        TObject* obj = key->ReadObj();
        //if (!obj->InheritsFrom("TH2D")) {
        //     delete obj;
        //     continue; 
        //}
        TH2D* hist2D = (TH2D*)obj;
        
        // Reset Branch Vectors for this new event
        b_hist_name = hist2D->GetName();
        b_event_idx = counter;
        b_channels.clear();
        b_rms.clear();

        int nBinsX = hist2D->GetNbinsX();
        int nBinsY = hist2D->GetNbinsY();
        
        // Determine Channel Base from Histogram Name
        int channel_base = 0;
        string h_name_str = b_hist_name;
        
        // Safety check on string length
        if (h_name_str.length() > 5) {
            char plane = h_name_str[1]; // 'u', 'v', 'w'
            char tpc = h_name_str.back(); // '0', '1' or others
            
            if (tpc == '0') {
                if (plane == 'u') channel_base = 0;
                else if (plane == 'v') channel_base = 1984;
                else if (plane == 'w') channel_base = 1984*2;
            } else if (tpc == '1') {
                if (plane == 'u') channel_base = 1984*2 + 1670;
                else if (plane == 'v') channel_base = 1984*3 + 1670;
                else if (plane == 'w') channel_base = 1984*4 + 1670;
            }
        }

        // Loop over wires in this event
        for (int x = 0; x <= nBinsX; ++x) {
            wire_ped.clear();
            wire.clear();
            
            // Extract waveform
            for (int y = 0; y <= nBinsY; ++y) {
                double binContent = hist2D->GetBinContent(x, y);
                // Simple compression: don't store pure zeros if that's how the TH2D is packed
                if (binContent != 0) {
                    wire_ped.push_back(binContent);
                    wire.push_back(binContent);
                }
            }
            
            // Skip empty or partial wires
            if (wire.size() < (size_t)(event_len - 10)) continue; 

            // Calculate Pedestal & Subtract
            double pedestal = Median(wire_ped); 
            
            double max_val = 0;
            for(size_t k=0; k<wire.size(); k++) {
                wire[k] -= pedestal;
                if(std::abs(wire[k]) > max_val) max_val = std::abs(wire[k]);
            }

            // Signal Rejection (Skip wires that have actual physics signals)
            if (max_val > 20.0) continue;
            
            // Calculate Channel ID
            int channel = x + channel_base - 1;
            if (channel < 0 || channel >= 12276) continue;

            // Calculate RMS
            float rms_val = (float)Noise_levels(wire);

            // Save to vectors
            b_channels.push_back(channel);
            b_rms.push_back(rms_val);

            // Note: I have skipped FFT here. 
            // If you want FFT per wire, you would need a vector<vector<float>> branch.
        }
         
        // Fill the Tree for this Histogram/Event
        if (!b_channels.empty()) {
            tree->Fill();
        }

        delete hist2D; 
        counter++;

        // Progress Bar
        if (counter % 10 == 0 || counter == total_keys) { 
            Double_t time_elapsed = timer.RealTime();
            timer.Continue(); 
            Double_t rate = counter / time_elapsed; 
            Double_t time_remaining = (total_keys - counter) / rate;
            
            std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% "
                      << "| Rate: " << std::fixed << std::setprecision(1) << rate << " keys/s "
                      << "| ETA: " << int(time_remaining) << "s    " << std::flush;
        }
    }

    std::cout << std::endl;
    //delete[] fft_input_buffer;
    //delete fft_planner; 

    // Save and Close
    cout << "Writing to file..." << endl;
    outFile->Write();
    outFile->Close();
    cout << "Done." << endl;
}

void TPC_WC_noise_analysis_full_wire(TString inputFile="sbnd-data-check.root")
{   
    cout<<"Starting Per-Event RMS Analysis..."<<endl;
    TFile *inFile = TFile::Open(inputFile.Data());
    if (!inFile || inFile->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }
    cout<<"Got File"<<endl;
    LoadRawDigits(inFile);
}