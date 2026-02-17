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
#include <omp.h> // Multithreading

using namespace std;

// --- OPTIMIZED HELPER FUNCTIONS ---

// Fast Median Calculation
// Uses nth_element for O(N) performance. 
// We pass 'vec' by value intentionally so we don't modify the original waveform.
double FastMedian(vector<double> vec) { 
    size_t n = vec.size();
    if (n == 0) return 0.0;
    
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    
    if (n % 2 != 0) return *target; // Odd length
    
    // Even length handling
    auto target_neighbor = std::max_element(vec.begin(), target);
    return (*target + *target_neighbor) / 2.0;
}

// Struct to hold results safely during parallel execution
struct WireResult {
    int channel_id;
    float rms;
    bool keep;
};

void LoadRawDigitsOptimized(TFile *inFile)
{   
    // --- 1. SETUP OUTPUT ---
    TFile* outFile = new TFile("noise_rms_fast.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "RMS per wire per event");

    // Branch Variables
    std::string b_hist_name;
    int b_event_idx;
    std::vector<int> b_channels;
    std::vector<float> b_rms;

    // Set branch addresses
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
        
        // Read Object
        TObject* obj = key->ReadObj();
        if (!obj->InheritsFrom("TH2D")) { 
            delete obj; 
            continue; 
        }
        TH2D* hist2D = (TH2D*)obj;
        
        // Reset Variables for this event
        b_hist_name = hist2D->GetName();
        b_event_idx = counter;
        b_channels.clear();
        b_rms.clear();

        int nWires = hist2D->GetNbinsX();
        int nTicks = hist2D->GetNbinsY();
        
        // Parse Channel Base (Once per event)
        int channel_base = 0;
        string h_name_str = b_hist_name;
        
        // Determine channel offset based on plane name (u/v/w)
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

        // --- 3. PARALLELIZED WIRE LOOP ---
        // Temporary storage vector sized to holding results for all wires
        // We use this because writing to TTree vectors inside parallel loop is unsafe.
        std::vector<WireResult> temp_results(nWires + 1);

        // OpenMP Parallel Loop: Splits wires among CPU cores
        #pragma omp parallel for schedule(dynamic)
        for (int x = 1; x <= nWires; ++x) {
            
            // Local storage for waveform
            vector<double> wire;
            wire.reserve(nTicks);
            
            // Extract waveform from TH2D
            for (int y = 1; y <= nTicks; ++y) {
                double val = hist2D->GetBinContent(x, y);
                // Simple compression: only store non-zeros
                if (val != 0) wire.push_back(val);
            }

            // Skip bad/empty wires (too few ticks)
            if (wire.size() < (size_t)(nTicks - 50)) {
                temp_results[x].keep = false;
                continue;
            }

            // Calculate Pedestal (Median)
            double pedestal = FastMedian(wire);

            // COMBINED LOOP: Subtract Pedestal, Find Max, Calculate SumSq
            double sum_sq = 0.0;
            double max_val = 0.0;
            
            for(double &val : wire) {
                val -= pedestal; // Subtract in place
                double abs_val = std::abs(val);
                if(abs_val > max_val) max_val = abs_val;
                sum_sq += val * val;
            }

            // Signal Rejection (Skip wires with big physics pulses)
            if (max_val > 20.0) {
                temp_results[x].keep = false;
                continue;
            }

            // Calculate Final RMS
            double rms = 0.0;
            if (wire.size() > 0) rms = std::sqrt(sum_sq / wire.size());

            // Store Result in Thread-Safe Array
            temp_results[x].channel_id = (x - 1) + channel_base;
            temp_results[x].rms = (float)rms;
            temp_results[x].keep = true;
        } 
        // --- END PARALLEL SECTION ---

        // 4. AGGREGATE RESULTS (Serial)
        // Move valid results from temporary storage to TTree vectors
        for(int x = 1; x <= nWires; ++x) {
            if (temp_results[x].keep && temp_results[x].channel_id >= 0) {
                b_channels.push_back(temp_results[x].channel_id);
                b_rms.push_back(temp_results[x].rms);
            }
        }

        // Only save event if we found valid wires
        if (!b_channels.empty()) tree->Fill();
        
        delete hist2D; 
        counter++;

        // Progress Update
        if (counter % 10 == 0 || counter == total_keys) {
            Double_t time_elapsed = timer.RealTime();
            timer.Continue();
            
            // Avoid division by zero
            if (time_elapsed < 0.001) time_elapsed = 0.001;
            
            Double_t rate = counter / time_elapsed;
            Double_t time_remaining = (total_keys - counter) / rate;

            std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% "
                      << "| Rate: " << std::fixed << std::setprecision(1) << rate << " ev/s " 
                      << "| ETA: " << (int)time_remaining << "s   " << std::flush;
        }
    }

    std::cout << std::endl;
    outFile->Write();
    outFile->Close();
    std::cout << "Done. Output saved to noise_rms_fast.root" << std::endl;
}

void TPC_WC_noise_analysis_full_wire_parallel(TString inputFile="sbnd-data-check.root")
{   
    TFile *inFile = TFile::Open(inputFile.Data());
    if (!inFile || inFile->IsZombie()) {
        cout << "Error opening file: " << inputFile << endl;
        return;
    }
    LoadRawDigitsOptimized(inFile);
}