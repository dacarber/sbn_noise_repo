#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include "TString.h"
#include "TTree.h"
#include <algorithm>
#include <TH2D.h>
#include "TStopwatch.h"
#include "TKey.h"

using namespace std;

// --- OPTIMIZED MEDIAN ---
// This sorts the vector IN PLACE. 
// It destroys the time-ordering of the waveform, but preserves the values.
// This is fine for RMS/Max calculation and saves copying memory.
double GetMedianInPlace(vector<double>& vec) { 
    size_t n = vec.size();
    if (n == 0) return 0.0;
    
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    
    if (n % 2 != 0) return *target; 
    
    // For even numbers, we need the average of the two middle elements
    auto target_neighbor = std::max_element(vec.begin(), target);
    return (*target + *target_neighbor) / 2.0;
}

void ProcessFast(TFile *inFile)
{   
    // 1. Setup Output
    TFile* outFile = new TFile("noise_rms_fast_serial.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "RMS per wire per event");

    // Variables
    std::string b_hist_name;
    int b_event_idx;
    std::vector<int> b_channels;
    std::vector<float> b_rms;

    tree->Branch("hist_name", &b_hist_name);
    tree->Branch("event_idx", &b_event_idx, "event_idx/I");
    tree->Branch("channels", &b_channels);
    tree->Branch("rms", &b_rms);

    // 2. Loop Prep
    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    int counter = 0;
    
    // PRE-ALLOCATE MEMORY to avoid creating/destroying vectors millions of times
    vector<double> waveform_buffer;
    waveform_buffer.reserve(4000); 

    TStopwatch timer;
    timer.Start();

    while ((key = (TKey*)next())) {
        
        TObject* obj = key->ReadObj();
        //if (!obj->InheritsFrom("TH2D")) { delete obj; continue; }
        TH2D* hist2D = (TH2D*)obj;
        
        b_hist_name = hist2D->GetName();
        b_event_idx = counter;
        b_channels.clear();
        b_rms.clear();

        int nWires = hist2D->GetNbinsX();
        int nTicks = hist2D->GetNbinsY();
        
        // --- THE SECRET WEAPON: RAW ARRAY ACCESS ---
        // Instead of GetBinContent(), we get the pointer to the raw data block.
        // ROOT stores TH2D data linearly: index = y * (nWires + 2) + x
        std::cout<<"1"<<std::endl;
        double* raw_data = hist2D->GetArray();
        std::cout<<"2"<<std::endl;
        int stride = nWires + 2; // +2 includes overflow/underflow bins

        // Parse Channel Base
        int channel_base = 0;
        string h_name = b_hist_name;
        if (h_name.length() > 5) {
            char plane = h_name[1]; 
            char tpc = h_name.back();
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

        // Loop over wires (Columns)
        for (int x = 1; x <= nWires; ++x) {
            
            // Clear the buffer, but keep the memory reserved (very fast)
            waveform_buffer.clear(); 

            // Copy raw data from the array stride
            // This is much faster than function calls
            for (int y = 1; y <= nTicks; ++y) {
                double val = raw_data[y * stride + x];
                if (val != 0) waveform_buffer.push_back(val);
            }

            if (waveform_buffer.size() < (size_t)(nTicks - 50)) continue;

            // 1. Calculate Pedestal (Sorts buffer in place)
            // Note: After this, the waveform is SCRAMBLED (not time-ordered).
            double pedestal = GetMedianInPlace(waveform_buffer);

            // 2. Calculate RMS & Max on the scrambled buffer
            // (RMS doesn't care about time order, so this is valid!)
            double sum_sq = 0.0;
            double max_val = 0.0;
            
            for (double val : waveform_buffer) {
                val -= pedestal;
                double abs_val = std::abs(val);
                if (abs_val > max_val) max_val = abs_val;
                sum_sq += val * val;
            }

            // Signal Rejection
            if (max_val > 20.0) continue;

            // Save
            int ch_id = (x - 1) + channel_base;
            if (ch_id >= 0) {
                b_channels.push_back(ch_id);
                b_rms.push_back((float)std::sqrt(sum_sq / waveform_buffer.size()));
            }
        } 

        if (!b_channels.empty()) tree->Fill();
        delete hist2D; 
        counter++;

        // Simple Progress Bar
        if (counter % 10 == 0) {
             std::cout << "\rProcessed " << counter << " / " << total_keys << " events" << std::flush;
        }
    }

    std::cout << std::endl << "Saving..." << std::endl;
    outFile->Write();
    outFile->Close();
}

void TPC_WC_noise_analysis_full_wire(TString inputFile="sbnd-data-check.root")
{   
    TFile *inFile = TFile::Open(inputFile.Data());
    if (!inFile || inFile->IsZombie()) return;
    ProcessFast(inFile);
}