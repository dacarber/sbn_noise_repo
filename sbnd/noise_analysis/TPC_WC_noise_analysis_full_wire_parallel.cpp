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

// --- OPTIMIZED HELPER FUNCTIONS ---

// Fast Median Calculation (O(N) performance)
double FastMedian(vector<double> vec) { 
    size_t n = vec.size();
    if (n == 0) return 0.0;
    
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    
    if (n % 2 != 0) return *target; 
    
    auto target_neighbor = std::max_element(vec.begin(), target);
    return (*target + *target_neighbor) / 2.0;
}

void LoadRawDigitsOptimized(TFile *inFile)
{   
    // --- 1. SETUP OUTPUT ---
    TFile* outFile = new TFile("noise_rms_optimized.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "RMS per wire per event");

    // Branch Variables
    std::string b_hist_name;
    int b_event_idx;
    std::vector<int> b_channels;
    std::vector<float> b_rms;

    tree->Branch("hist_name", &b_hist_name);
    tree->Branch("event_idx", &b_event_idx, "event_idx/I");
    tree->Branch("channels", &b_channels);
    tree->Branch("rms", &b_rms);

    // --- 2. LOOP PREP ---
    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    int counter = 0;
    
    TStopwatch timer;
    timer.Start();

    // Iterate over histograms (Events)
    while ((key = (TKey*)next())) {
        
        TObject* obj = key->ReadObj();
        if (!obj->InheritsFrom("TH2D")) { delete obj; continue; }
        TH2D* hist2D = (TH2D*)obj;
        
        // Reset Variables
        b_hist_name = hist2D->GetName();
        b_event_idx = counter;
        b_channels.clear();
        b_rms.clear();

        int nWires = hist2D->GetNbinsX();
        int nTicks = hist2D->GetNbinsY();
        
        // Parse Channel Base
        int channel_base = 0;
        string h_name_str = b_hist_name;
        if (h_name_str.length() > 5) {
            char plane = h_name_str[1]; 
            char tpc = h_name_str.back();
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

        // --- 3. OPTIMIZED SERIAL LOOP ---
        for (int x = 1; x <= nWires; ++x) {
            
            vector<double> wire;
            wire.reserve(nTicks);
            
            for (int y = 1; y <= nTicks; ++y) {
                double val = hist2D->GetBinContent(x, y);
                if (val != 0) wire.push_back(val);
            }

            if (wire.size() < (size_t)(nTicks - 50)) continue;

            // Fast Calculation
            double pedestal = FastMedian(wire);
            double sum_sq = 0.0;
            double max_val = 0.0;
            
            for(double &val : wire) {
                val -= pedestal; 
                double abs_val = std::abs(val);
                if(abs_val > max_val) max_val = abs_val;
                sum_sq += val * val;
            }

            if (max_val > 20.0) continue; // Skip signal

            // Save Result
            int ch_id = (x - 1) + channel_base;
            if (ch_id >= 0) {
                b_channels.push_back(ch_id);
                b_rms.push_back((float)std::sqrt(sum_sq / wire.size()));
            }
        } 

        if (!b_channels.empty()) tree->Fill();
        
        delete hist2D; 
        counter++;

        // Progress Update
        if (counter % 10 == 0 || counter == total_keys) {
            Double_t time_elapsed = timer.RealTime();
            timer.Continue();
            if (time_elapsed < 0.001) time_elapsed = 0.001;
            std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% "
                      << "| Rate: " << std::fixed << std::setprecision(1) << counter/time_elapsed << " ev/s " 
                      << std::flush;
        }
    }

    std::cout << std::endl;
    outFile->Write();
    outFile->Close();
    std::cout << "Done. Output saved to noise_rms_optimized.root" << std::endl;
}

void TPC_WC_noise_analysis_full_wire_parallel(TString inputFile="sbnd-data-check.root")
{   
    TFile *inFile = TFile::Open(inputFile.Data());
    if (!inFile || inFile->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }
    LoadRawDigitsOptimized(inFile);
}