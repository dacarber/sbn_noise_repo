#include "TFile.h"
#include "TH2D.h"
#include "TMatrixDSym.h"
#include "TString.h"
#include "TKey.h"
#include "TStopwatch.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

// --- Helper: Median Calculation for Pedestal ---
double GetMedian(vector<double> vec) { 
    if (vec.empty()) return 0.0;
    size_t n = vec.size();
    auto target = vec.begin() + n / 2;
    std::nth_element(vec.begin(), target, vec.end());
    double median = *target;
    if (n % 2 == 0) {
        auto target_neighbor = std::max_element(vec.begin(), target);
        median = (*target + *target_neighbor) / 2.0;
    }
    return median;
}

void ProcessCovariance(TString inputFile)
{
    // Open Input
    TFile *inFile = TFile::Open(inputFile);
    if (!inFile || inFile->IsZombie()) {
        cout << "Error opening input file!" << endl;
        return;
    }

    // Open Output
    TFile *outFile = new TFile("covariance_matrices.root", "RECREATE");
    
    // Iteration Setup
    TIter next(inFile->GetListOfKeys());
    TKey* key;
    int counter = 0;
    int total_keys = inFile->GetListOfKeys()->GetSize();
    
    TStopwatch timer;
    timer.Start();

    cout << "Starting Covariance Calculation..." << endl;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        
        // Ensure we are looking at a 2D Histogram (Event)
        //if (!obj->InheritsFrom("TH2D")) {
        //     delete obj;
        //     continue; 
        //}
        TH2D* hist2D = (TH2D*)obj;
        TString histName = hist2D->GetName();

        // --- Data Extraction ---
        int nWires = hist2D->GetNbinsX(); // Number of wires in this plane/hist
        int nTicks = hist2D->GetNbinsY(); // Number of time ticks

        // We use a vector of vectors to store the "cleaned" waveforms
        // Only storing valid wires to save memory before matrix calc
        // Format: data_cache[wire_index][tick_index]
        vector<vector<double>> data_cache(nWires, vector<double>(nTicks));
        vector<bool> is_wire_good(nWires, false);

        // 1. Load and Clean Data (Pedestal Subtraction)
        for (int w = 0; w < nWires; ++w) {
            vector<double> raw_wave;
            raw_wave.reserve(nTicks);
            
            // Extract single wire
            for (int t = 1; t <= nTicks; ++t) {
                // TH2 binning starts at 1
                raw_wave.push_back(hist2D->GetBinContent(w+1, t));
            }

            // Calculate Pedestal (Median)
            double pedestal = GetMedian(raw_wave);
            
            // Subtract Pedestal & Check for Signal
            double max_amp = 0.0;
            for (int t = 0; t < nTicks; ++t) {
                double val = raw_wave[t] - pedestal;
                data_cache[w][t] = val;
                if (std::abs(val) > max_amp) max_amp = std::abs(val);
            }

            // FILTER: If signal is too high (physics), mark as bad so we don't skew noise correlation
            // If you WANT physics correlations, comment out the check below.
            if (max_amp < 20.0 && max_amp > 0.001) {
                is_wire_good[w] = true;
            } else {
                // Fill with zeros so it doesn't affect sums
                std::fill(data_cache[w].begin(), data_cache[w].end(), 0.0);
                is_wire_good[w] = false; 
            }
        }

        // --- Covariance Calculation ---
        // We use TMatrixDSym because Covariance matrices are Symmetric (Cov[i][j] == Cov[j][i])
        // This saves half the memory.
        TMatrixDSym covMatrix(nWires);

        // Loop over all pairs of wires (Upper triangle only)
        for (int i = 0; i < nWires; ++i) {
            
            // Optimization: If wire i is empty/bad, the whole row is 0
            if (!is_wire_good[i]) continue; 

            for (int j = i; j < nWires; ++j) {
                
                if (!is_wire_good[j]) continue;

                double dot_product = 0.0;
                
                // Vectorized loop (compiler usually optimizes this well)
                for (int t = 0; t < nTicks; ++t) {
                    dot_product += data_cache[i][t] * data_cache[j][t];
                }

                double covariance = dot_product / (double)nTicks;
                
                // Set matrix element (indexing is valid for TMatrix)
                covMatrix(i, j) = covariance;
                covMatrix(j, i) = covariance; // Symmetric fill
            }
        }

        // --- Saving ---
        // We create a directory for each histogram to keep things organized
        // or just append the name. 
        outFile->cd();
        covMatrix.Write(Form("cov_%s", histName.Data()));

        // Cleanup
        delete hist2D;
        counter++;

        // --- Progress Bar ---
        if (counter % 1 == 0) { // Update every event because this is slow
            Double_t time_elapsed = timer.RealTime();
            timer.Continue();
            Double_t rate = counter / time_elapsed;
            Double_t time_remaining = (total_keys - counter) / (rate + 0.0001);
            
            cout << "\rProcessing: " << histName 
                 << " (" << nWires << "x" << nWires << ")"
                 << " | Rate: " << std::fixed << std::setprecision(2) << rate << " ev/s"
                 << " | ETA: " << (int)time_remaining << "s   " << flush;
        }
    }

    cout << endl << "Writing and Closing..." << endl;
    outFile->Close();
    inFile->Close();
    cout << "Done. Output saved to covariance_matrices.root" << endl;
}