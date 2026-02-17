#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TKey.h"
#include "TStopwatch.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

// --- Helper: Fast Median ---
template <typename T>
double GetMedianInPlace(vector<T>& vec) { 
    size_t n = vec.size();
    if (n == 0) return 0.0;
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    if (n % 2 != 0) return (double)*target; 
    auto target_neighbor = std::max_element(vec.begin(), target);
    return ((double)*target + (double)*target_neighbor) / 2.0;
}

// --- Templated Processor ---
template <typename HistType, typename DataType>
void ProcessHistogram(TObject* obj, TTree* tree, 
                      int& b_event_idx, string& b_hist_name, 
                      vector<int>& b_channels, vector<float>& b_rms) 
{
    HistType* hist = (HistType*)obj;
    
    // Metadata
    b_hist_name = hist->GetName();
    b_channels.clear();
    b_rms.clear();

    int nWires = hist->GetNbinsX();
    int nTicks = hist->GetNbinsY();
    
    // --- SAFETY CHECKS ---
    // 1. Get the raw array pointer
    DataType* raw_data = hist->GetArray();
    
    // 2. Calculate the required size
    // ROOT TH2 structure: (XBins + 2) * (YBins + 2)
    Long64_t required_size = (Long64_t)(nWires + 2) * (Long64_t)(nTicks + 2);
    int stride = nWires + 2;

    // 3. Determine if we can use Fast Access
    bool use_fast_access = true;
    
    if (!raw_data) {
        // Pointer is null -> Empty histogram or weird allocation
        use_fast_access = false; 
    } 
    else if (hist->GetSize() < required_size) {
        // The array is smaller than the bins say it should be.
        // This indicates a corrupted header.
        // We will fallback to standard GetBinContent which handles bounds checking.
        use_fast_access = false;
    }

    // --- DECODING CHANNEL MAPPING ---
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

    // Pre-allocate buffer
    vector<DataType> waveform_buffer;
    waveform_buffer.reserve(nTicks);

    // --- LOOP ---
    for (int x = 1; x <= nWires; ++x) {
        waveform_buffer.clear();

        if (use_fast_access) {
            // FAST PATH: Pointer Arithmetic
            // We iterate Y (ticks) for a fixed X (wire)
            for (int y = 1; y <= nTicks; ++y) {
                // Determine index
                Long64_t idx = (Long64_t)y * stride + x;
                DataType val = raw_data[idx];
                if (val != 0) waveform_buffer.push_back(val);
            }
        } else {
            // SLOW (SAFE) PATH: Standard ROOT function calls
            // Use this if the array pointer looked suspicious
            for (int y = 1; y <= nTicks; ++y) {
                DataType val = (DataType)hist->GetBinContent(x, y);
                if (val != 0) waveform_buffer.push_back(val);
            }
        }

        // Quality Cuts
        if (waveform_buffer.size() < (size_t)(nTicks - 50)) continue;

        // Calc Stats
        double pedestal = GetMedianInPlace(waveform_buffer);
        
        double sum_sq = 0.0;
        double max_val = 0.0;

        for (DataType val : waveform_buffer) {
            double v = (double)val - pedestal;
            if (std::abs(v) > max_val) max_val = std::abs(v);
            sum_sq += v * v;
        }

        if (max_val > 20.0) continue; 

        // Save
        int ch_id = (x - 1) + channel_base;
        if (ch_id >= 0) {
            b_channels.push_back(ch_id);
            b_rms.push_back((float)std::sqrt(sum_sq / waveform_buffer.size()));
        }
    }
    
    if (!b_channels.empty()) tree->Fill();
}

void ProcessSafe(TFile *inFile)
{   
    TFile* outFile = new TFile("noise_rms_robust.root", "RECREATE");
    TTree* tree = new TTree("tpc_noise", "RMS per wire per event");

    std::string b_hist_name;
    int b_event_idx;
    std::vector<int> b_channels;
    std::vector<float> b_rms;

    tree->Branch("hist_name", &b_hist_name);
    tree->Branch("event_idx", &b_event_idx, "event_idx/I");
    tree->Branch("channels", &b_channels);
    tree->Branch("rms", &b_rms);

    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int counter = 0;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    
    TStopwatch timer;
    timer.Start();

    while ((key = (TKey*)next())) {
        
        TObject* obj = key->ReadObj();
        
        // Handle TH2D (Double Precision)
        if (obj->InheritsFrom("TH2D")) {
            b_event_idx = counter;
            ProcessHistogram<TH2D, double>(obj, tree, b_event_idx, b_hist_name, b_channels, b_rms);
        } 
        // Handle TH2F (Float Precision)
        else if (obj->InheritsFrom("TH2F")) {
            b_event_idx = counter;
            ProcessHistogram<TH2F, float>(obj, tree, b_event_idx, b_hist_name, b_channels, b_rms);
        }
        
        delete obj;
        counter++;

        if (counter % 10 == 0) {
             double rate = counter / (timer.RealTime() + 0.001);
             timer.Continue();
             std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% " 
                       << "| Rate: " << (int)rate << " ev/s" << std::flush;
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
    ProcessSafe(inFile);
}