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
#include <sys/resource.h>
#include <iomanip>
#include <TKey.h>

using namespace std;

// --- HELPER FUNCTIONS ---

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

void FFT(const vector<double>& noise_channel, vector<double>& output_mag, TVirtualFFT* fft_planner, double* input_buffer) {
    int vec_size = noise_channel.size();
    
    for (size_t i = 0; i < vec_size; ++i) {
        input_buffer[i] = noise_channel[i];
    }

    fft_planner->SetPoints(input_buffer);
    fft_planner->Transform();

    double fftReal, fftImag;
    
    // Fill the output vector
    for(size_t k=1; k < vec_size / 2 + 2; k++){
        fft_planner->GetPointComplex(k, fftReal, fftImag);
        if (k == vec_size / 2 + 1){
            output_mag[k] = 1.0;
        } else {
            output_mag[k] = TMath::Sqrt(fftReal*fftReal + fftImag*fftImag);
        }
    }
}

// --- MAIN ANALYSIS FUNCTION ---

void LoadRawDigits(TFile *inFile)
{   
    cout << "Setting up Output Tree..." << endl;

    // Output File
    TFile* outFile = new TFile("noise_output_event_avg.root", "RECREATE");
    TTree* tree = new TTree("event_noise", "Average Noise per Event/Plane");

    // --- TTree Branch Variables ---
    std::string b_hist_name;
    std::string b_type; // "raw" or "orig"
    int b_tpc;          // 0 or 1
    int b_plane;        // 0=u, 1=v, 2=w
    float b_avg_rms;    // The average RMS of all wires in this event
    int b_num_wires;    // How many wires contributed to this average
    vector<float> b_avg_fft; // The average FFT spectrum for this event

    // --- Branch Definitions ---
    tree->Branch("hist_name", &b_hist_name);
    tree->Branch("type", &b_type);
    tree->Branch("tpc", &b_tpc, "tpc/I");
    tree->Branch("plane", &b_plane, "plane/I");
    tree->Branch("avg_rms", &b_avg_rms, "avg_rms/F");
    tree->Branch("num_wires", &b_num_wires, "num_wires/I");
    tree->Branch("avg_fft", &b_avg_fft);

    // --- FFT Setup ---
    int event_len = 3427; 
    Int_t n_size = event_len;
    TVirtualFFT* fft_planner = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
    double* fft_input_buffer = new double[event_len];
    vector<double> fft_output_buffer(event_len/2 + 2);

    // Temp vectors for processing
    vector<double> wire_ped; 
    wire_ped.reserve(event_len);
    vector<double> wire; 
    wire.reserve(event_len);

    // Iterators
    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    int counter = 0;
    TStopwatch timer;
    timer.Start();

    // --- MAIN READING LOOP ---
    while ((key = (TKey*)next())) {
        
        // 1. Read Object
        TObject* obj = key->ReadObj();
        // Safety check to ensure it's a histogram
        //if (!obj->InheritsFrom("TH2D")) {
        //    delete obj;
        //    continue; 
        //}
        TH2D* hist2D = (TH2D*)obj;
        b_hist_name = hist2D->GetName();

        // 2. Parse Metadata (TPC, Plane, Type) from Name
        // Assuming format like: "h_orig_..._plane_u_tpc_0"
        
        // Determine Type
        if (b_hist_name.find("orig") != std::string::npos) b_type = "orig";
        else if (b_hist_name.find("raw") != std::string::npos) b_type = "raw";
        else b_type = "unknown";

        // Determine Plane
        char plane_char = b_hist_name[1]; // Based on your old logic
        if (plane_char == 'u') b_plane = 0;
        else if (plane_char == 'v') b_plane = 1;
        else if (plane_char == 'w') b_plane = 2;
        else b_plane = -1;

        // Determine TPC
        char tpc_char = b_hist_name.back();
        if (tpc_char == '0') b_tpc = 0;
        else if (tpc_char == '1') b_tpc = 1;
        else b_tpc = -1;

        // 3. Reset Event Accumulators
        double event_sum_rms = 0.0;
        b_num_wires = 0;
        
        // We need a vector to sum the FFTs. Initialize with 0s.
        // Size depends on output of FFT function
        vector<double> event_sum_fft(event_len/2 + 2, 0.0);

        int nBinsX = hist2D->GetNbinsX();
        int nBinsY = hist2D->GetNbinsY();

        // 4. Loop over Wires (X-axis) in this Histogram
        for (int x = 0; x <= nBinsX; ++x) {
            wire_ped.clear();
            wire.clear();

            // Extract Waveform
            for (int y = 0; y <= nBinsY; ++y) {
                double binContent = hist2D->GetBinContent(x, y);
                // Zero-suppression check (optional, kept from your code)
                if (binContent == 0) continue; 
                wire_ped.push_back(binContent);
                wire.push_back(binContent);
            }

            if (wire.size() != event_len) continue;

            // Pedestal Subtraction
            int pedestal = Median(wire_ped);
            transform(wire.begin(), wire.end(), wire.begin(), 
                      [pedestal](double elem) { return elem - pedestal; });

            // Signal Rejection (Skip wires with large pulses)
            double max_val = 0;
            for(double w : wire) if(std::abs(w) > max_val) max_val = std::abs(w);
            if (max_val > 20.0) continue; 

            // --- ACCUMULATE STATS ---
            
            // RMS
            double this_rms = Noise_levels(wire);
            event_sum_rms += this_rms;

            // FFT
            FFT(wire, fft_output_buffer, fft_planner, fft_input_buffer);
            for(size_t k=0; k < fft_output_buffer.size(); k++) {
                if(k < event_sum_fft.size()) {
                    event_sum_fft[k] += fft_output_buffer[k];
                }
            }

            b_num_wires++;
        } // End Wire Loop

        // 5. Finalize and Fill Tree
        if (b_num_wires > 0) {
            b_avg_rms = event_sum_rms / b_num_wires;

            // Average the FFT vector
            b_avg_fft.clear();
            for(double val : event_sum_fft) {
                b_avg_fft.push_back( val / b_num_wires );
            }

            tree->Fill();
        }

        // Cleanup current object
        delete hist2D;

        // Progress Bar
        counter++;
        if (counter % 10 == 0 || counter == total_keys) { 
            Double_t time_elapsed = timer.RealTime();
            timer.Continue();
            Double_t rate = counter / time_elapsed;
            Double_t time_remaining = (total_keys - counter) / rate;
            std::cout << "\rProgress: " << int((float)counter/total_keys * 100) << "% "
                      << "| Rate: " << std::fixed << std::setprecision(1) << rate << " keys/s "
                      << "| ETA: " << int(time_remaining) << "s    " << std::flush;
        }
    } // End Key Loop

    std::cout << std::endl << "Writing file..." << std::endl;
    delete[] fft_input_buffer;
    delete fft_planner; 
    
    outFile->Write();
    outFile->Close();
    cout << "Done. Saved to noise_output_event_avg.root" << endl;
}

void TPC_WC_noise_analysis_v2_avg_det(TString inputFile="sbnd-data-check.root")
{   
    cout<<"Starting Event-Wise Analysis..."<<endl;
    TFile *inFile = TFile::Open(inputFile.Data());
    if(!inFile || inFile->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }
    cout<<"File Opened."<<endl;
    LoadRawDigits(inFile);
}