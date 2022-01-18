#include <iostream>
#include <string>
#include "TTree.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TGraph.h"
#include <map>
#include <vector>
#include <ROOT/TProcessExecutor.hxx>
#include <tuple>
#include "TRandom3.h"


std::map<Double_t, Double_t> Beam_MaxEnergy = {{1.0, 100}, {2.0, 100}, {3.0, 150}, {5.0, 200}, {10.0, 400}, {20.0, 700}, {30.0, 1000},
     {40.0, 1500}, {50.0, 2000}, {60., 2500}, {70., 3000}, {80., 3500}, {90., 4000}, {100., 4500}}; // Determined arbitrarily; just an energy that includes entire distribution in the histogram
const Int_t num_total_layers = 51;
const Int_t num_towers = 36;
const Int_t num_tail = 3;
Bool_t EnableTailCatcher = kTRUE; // Tail Catcher used for hadrons


TH1D* ECalWeighting(std::vector<Double_t> ECalEdep_event, std::vector<Double_t> HCalEdep_event, std::vector<Double_t> TailCatcherEdep_event, Double_t energy, Double_t weight)
{
    const Int_t num_events = ECalEdep_event.size();
    char hist_name[100];
    sprintf(hist_name, "h_TotalEdep%f", weight);

    TH1D* h_TotalEdep_weighted = new TH1D(hist_name, "", 600, 0, Beam_MaxEnergy[energy]);

    for(Int_t i = 0; i < num_events; i++)
    {
        Double_t total_edep = ECalEdep_event[i]/weight + HCalEdep_event[i];
        Double_t fraction = 1.;
        if( (total_edep) != 0.) fraction = TailCatcherEdep_event[i] / total_edep; // Fraction of energy deposited in tail catcher region
        if(EnableTailCatcher && fraction < 0.01) // 0.01 threshold might be adjusted
        {
            h_TotalEdep_weighted->Fill(total_edep);
        }
        else if(!EnableTailCatcher) h_TotalEdep_weighted->Fill(total_edep);
    }

    return h_TotalEdep_weighted;
}
TH1D* ECalWeightingProcess(TTree* TotalTree, TTree* HCalTree, Double_t energy)
{

    Double_t ECalEdep, HCalTileEdep;
    Int_t ECal_EventID, HCal_EventID, HCal_LayerID;

    TotalTree->SetBranchAddress("ECal_Edep_Total", &ECalEdep); // Total ECal energy per event
    TotalTree->SetBranchAddress("eventID", &ECal_EventID);
    
    HCalTree->SetBranchAddress("HCal_Edep_Tile", &HCalTileEdep); // Energy per tile in each event
    HCalTree->SetBranchAddress("HCal_Layerid", &HCal_LayerID); // Layer number that the tile belongs to
    HCalTree->SetBranchAddress("eventID", &HCal_EventID);

    Int_t num_events = (Int_t) TotalTree->GetEntries();

    std::vector<Double_t> ECalEdep_event;
    std::vector<Double_t> HCalEdep_event;
    std::vector<Double_t> TailCatcherEdep_event; // Tail-Catcher layers are last num_tail (currently 3) layers of HCal

    for(Int_t i = 0; i < num_events; i++)
    {
        ECalEdep_event.push_back(0.);
        HCalEdep_event.push_back(0.);
        TailCatcherEdep_event.push_back(0.);
    }

    for(Int_t i = 0; i < num_events; i++)
    {
        TotalTree->GetEntry(i);
        ECalEdep_event[ECal_EventID] += ECalEdep;
        
        for(Int_t itile = 0; itile < num_towers*num_total_layers; itile++) 
        {
            HCalTree->GetEntry(itile + i*num_towers*num_total_layers); // 51 layers in each tower, 1 tile per layer in each tower. 36 towers. -> 36*51 tiles to read
            HCalTileEdep *= gRandom->Gaus(1., 0.2); // Smearing
            if(HCalTileEdep < 0.5) HCalTileEdep = 0.; // Tile cut
            HCalEdep_event[HCal_EventID] += HCalTileEdep;

            if( (HCal_LayerID + 1) > num_total_layers - num_tail ) // Tail Catcher (last 3 layers)
            {
                TailCatcherEdep_event[HCal_EventID] += HCalTileEdep;
            }
        }
    }

    std::vector<Double_t> weights;
    Double_t w0 = 0.6;
    for(Int_t i = 0; i < 20; i++)
    {
        weights.push_back(w0 + i*0.05); // Weights used for ECal weighting
    }

    std::vector<TH1D*> WeightedHists;

    for(auto& weight : weights)
    {
        WeightedHists.push_back(ECalWeighting(ECalEdep_event, HCalEdep_event, TailCatcherEdep_event, energy, weight));
    }

    TCanvas* c_WeightedHists = new TCanvas("c_WeightedHists", "", 1000, 1000);
    c_WeightedHists->Divide(4,5);

    Double_t minimum_resolution = 99999.;
    Double_t optimal_weight = -1.;
    Int_t optimal_index = -1;
    std::vector<Double_t> resolutions;

    // Looping through all weights and obtaining weight that minimizes the resolution
    for(Int_t i = 0; i < WeightedHists.size(); i++)
    {
        c_WeightedHists->cd(i+1);
        WeightedHists[i]->Draw("same"); // Graphing the energy spectrums for each weight
        WeightedHists[i]->Fit("gaus", "q");
        Double_t mean = WeightedHists[i]->GetFunction("gaus")->GetParameter(1);
        Double_t sigma = WeightedHists[i]->GetFunction("gaus")->GetParameter(2);
        Double_t resolution = 99999.;
        if(mean != 0) resolution = sigma/mean;
        resolutions.push_back(resolution);

        if( resolution < minimum_resolution )
        {
            minimum_resolution = resolution;
            optimal_weight = weights[i];
            optimal_index = i;
        } 
    }

    // Graph of resolution vs weight
    TCanvas* c_Weights = new TCanvas("c_Weights", "", 1200, 800);
    TGraph* gr_Weights = new TGraph(weights.size(), weights.data(), resolutions.data());
    gr_Weights->Draw("ap");
    gr_Weights->SetMarkerStyle(kOpenCircle);
    gr_Weights->SetMarkerColor(kBlue);
    gr_Weights->GetXaxis()->SetTitle("Weight");
    gr_Weights->GetYaxis()->SetTitle("Resolution");
    gr_Weights->SetTitle("");

    std::cout<<"Optimal weight is "<<optimal_weight<<std::endl;
    return WeightedHists[optimal_index];
}
Double_t Resolution(std::string particle = "e-", Double_t energy = 1.0)
{   
    Bool_t ECal_weight = kTRUE; // Weighting procedure used for hadrons
    if(particle == "e-" || particle == "pi0") ECal_weight = kFALSE;
    if(ECal_weight) std::cout<<"Using weighting for ECal."<<std::endl;
    else std::cout<<"Not using weighting for ECal."<<std::endl;

    TString file_name;
    file_name.Form("%s_QGSP/%s_%0.0fGeV.root", particle.c_str(), particle.c_str(), energy); // Output file from Geant4 simulation. Change to whatever name you have.
    std::cout<<"Opening "<<file_name<<std::endl;
    TFile* data_file = new TFile(file_name);

    TH1D* h_TotalEdep = new TH1D("h_TotalEdep", "", 600, 0, Beam_MaxEnergy[energy]);
    TTree* TotalTree = (TTree*) data_file->Get("EdepTotal"); // Tree that holds total energy for HCal and ECal
    TTree* HCalTree = (TTree*) data_file->Get("HCalLayers"); // Tree that holds individual tile information for HCal
    Int_t num_events = (Int_t) TotalTree->GetEntries();
    std::cout<<"Number of events: "<<num_events<<std::endl;

    if(ECal_weight) h_TotalEdep = (TH1D*) ECalWeightingProcess(TotalTree, HCalTree, energy)->Clone();
    else
    {
        Double_t ECalEdep, HCalTileEdep;
        Int_t ECal_EventID, HCal_EventID, HCal_LayerID;

        TotalTree->SetBranchAddress("ECal_Edep_Total", &ECalEdep); // Total energy in ECal per event
        TotalTree->SetBranchAddress("eventID", &ECal_EventID);
        
        HCalTree->SetBranchAddress("HCal_Edep_Tile", &HCalTileEdep); // Energy per tile in each event
        HCalTree->SetBranchAddress("HCal_Layerid", &HCal_LayerID); // Layer number that the tile belongs to
        HCalTree->SetBranchAddress("eventID", &HCal_EventID);

        Double_t ECalEdep_event[num_events];
        Double_t HCalEdep_event[num_events];

        for(Int_t i = 0; i < num_events; i++)
        {
            ECalEdep_event[i] = 0.;
            HCalEdep_event[i] = 0.;
        }

        for(Int_t i = 0; i < num_events; i++)
        {
            TotalTree->GetEntry(i);
            ECalEdep_event[ECal_EventID] += ECalEdep;
            
            for(Int_t itile = 0; itile < num_towers*num_total_layers; itile++)
            {
                HCalTree->GetEntry(itile + i*num_towers*num_total_layers); // 51 layers in each tower, 1 tile per layer in each tower. 36 towers. -> 36*51 tiles to read
                HCalTileEdep *= gRandom->Gaus(1., 0.2); // Smearing 
                if(HCalTileEdep < .5) HCalTileEdep = 0.; // 0.5 MeV cut on tile 
                HCalEdep_event[HCal_EventID] += HCalTileEdep;
            }
        }

        for(Int_t i = 0; i < num_events; i++)
        {
            Double_t total_edep = ECalEdep_event[i] + HCalEdep_event[i]; // Making sure HCal and ECal event ids are the same
            h_TotalEdep->Fill(total_edep);
        }
    }

    h_TotalEdep->SetTitle(""); // Set name and title since they changed if using the ECal weighting process
    h_TotalEdep->SetName("h_TotalEdep");
    
    TCanvas* c_Resolution = new TCanvas("c_Resolution", "", 1000, 1000);
    h_TotalEdep->Draw();
    h_TotalEdep->GetXaxis()->SetTitle("Edep (MeV)");
    h_TotalEdep->GetYaxis()->SetTitle("Number of Events");

    TF1* f_gaus = new TF1("f_gaus", "gaus", 0, Beam_MaxEnergy[energy]);
    h_TotalEdep->Fit(f_gaus, "");

    Double_t mean = f_gaus->GetParameter(1);
    Double_t sigma = f_gaus->GetParameter(2);
    Double_t resolution = 0.;
    if(mean != 0) resolution = sigma/mean;

    std::cout<<"Resolution is "<<resolution<<std::endl; // Use this resolution. Put it in Resolution_Plot.cpp

    TString info_text;
    info_text = Form("%s at %0.0f GeV", particle.c_str(), energy);

    h_TotalEdep->SetTitle(info_text);

    TString res_text;
    res_text = Form("Resolution = %0.5f", resolution);

    TLatex info_caption;
    info_caption.SetTextFont(62);
    info_caption.SetTextSize(.04);
    info_caption.SetNDC(kTRUE);
    info_caption.DrawLatex(.15, .85, info_text);

    TLatex res_caption;
    res_caption.SetTextFont(62);
    res_caption.SetTextSize(.04);
    res_caption.SetNDC(kTRUE);
    res_caption.DrawLatex(.15, .8, res_text);

    return resolution;
}
