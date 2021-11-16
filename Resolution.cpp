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


std::map<Double_t, Double_t> Beam_MaxEnergy = {{1.0, 100}, {2.0, 100}, {3.0, 150}, {5.0, 200}, {10.0, 400}, {20.0, 700}, {30.0, 1000}, {40.0, 1500}, {50.0, 2000}, {60., 2500}, {70., 3000}, {80., 3500}, {90., 4000}, {100., 4500}};
const Int_t num_total_layers = 51;
const Int_t num_tail = 3;


TH1D* ECalWeighting(TTree* tree, Double_t energy, Double_t weight)
{
    Double_t ECalEdep, HCalEdep;
    tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    tree->SetBranchAddress("HCal_Edep_Total", &HCalEdep);
    Int_t num_events = (Int_t) tree->GetEntries();
    char hist_name[100];
    sprintf(hist_name, "h_TotalEdep%f", weight);

    TH1D* h_TotalEdep = new TH1D(hist_name, "", 600, 0, Beam_MaxEnergy[energy]);
    for(Int_t i = 0; i < num_events; i++)
    {
        tree->GetEntry(i);
        if(ECalEdep < 0.183) ECalEdep = 0.;
        h_TotalEdep->Fill(ECalEdep/weight + HCalEdep);        
    }

    return h_TotalEdep;
}
Double_t ECalWeightingProcess(TTree* tree, Double_t energy)
{
    std::vector<Double_t> weights;
    Double_t w0 = 0.6;
    for(Int_t i = 0; i < 20; i++)
    {
        weights.push_back(w0 + i*0.05);
    }

    std::vector<TH1D*> WeightedHists;

    for(auto& weight : weights)
    {
        WeightedHists.push_back(ECalWeighting(tree, energy, weight));
    }

    TCanvas* c_WeightedHists = new TCanvas("c_WeightedHists", "", 1000, 1000);
    c_WeightedHists->Divide(4,5);

    Double_t minimum_resolution = 99999;
    Double_t optimal_weight = -1.;
    std::vector<Double_t> resolutions;

    for(Int_t i = 0; i < WeightedHists.size(); i++)
    {
        c_WeightedHists->cd(i+1);
        WeightedHists[i]->Draw("same");
        WeightedHists[i]->Fit("gaus", "q");
        Double_t mean = WeightedHists[i]->GetFunction("gaus")->GetParameter(1);
        Double_t sigma = WeightedHists[i]->GetFunction("gaus")->GetParameter(2);
        Double_t resolution;
        if(mean != 0) resolution = sigma/mean;
        else resolution = 99999.;
        resolutions.push_back(resolution);

        if( resolution < minimum_resolution)
        {
            minimum_resolution = resolution;
            optimal_weight = weights[i];
        } 
    }

    TCanvas* c_Weights = new TCanvas("c_Weights", "", 1000, 1000);
    TGraph* gr_Weights = new TGraph(weights.size(), weights.data(), resolutions.data());
    gr_Weights->Draw("ap");
    gr_Weights->SetMarkerStyle(kOpenCircle);
    gr_Weights->SetMarkerColor(kBlue);
    gr_Weights->GetXaxis()->SetTitle("Weight");
    gr_Weights->GetYaxis()->SetTitle("Resolution");
    gr_Weights->SetTitle("");


    std::cout<<"Optimal weight is "<<optimal_weight<<std::endl;
    return optimal_weight;
}

TH1D* TailCatcher(TTree* total_tree, TTree* HCal_tree, Double_t energy, Double_t optimal_weight)
{
    Double_t ECalEdep, HCalTileEdep, test;
    Int_t ECal_EventID, HCal_EventID, HCal_LayerID;

    total_tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    total_tree->SetBranchAddress("HCal_Edep_Total", &test);

    total_tree->SetBranchAddress("eventID", &ECal_EventID);
    
    HCal_tree->SetBranchAddress("HCal_Edep_Tile", &HCalTileEdep);
    HCal_tree->SetBranchAddress("HCal_Layerid", &HCal_LayerID);
    HCal_tree->SetBranchAddress("eventID", &HCal_EventID);

    Int_t num_events = (Int_t) total_tree->GetEntries();
    char hist_name[100];
    sprintf(hist_name, "h_TailCatcher");
    Double_t TotalTailCatch = 0.;
    TH1D* h_TotalEdep = new TH1D(hist_name, "", 600, 0, Beam_MaxEnergy[energy]);
    
    Double_t ECalEdep_array[num_events];
    Double_t HCalEdep_array[num_events];
    Double_t TailCatch_array[num_events];

    for(Int_t i = 0; i < num_events; i++)
    {
        ECalEdep_array[i] = 0.;
        HCalEdep_array[i] = 0.;
        TailCatch_array[i] = 0.;
    }

    for(Int_t i = 0; i < num_events; i++)
    {
        total_tree->GetEntry(i);
        if(ECalEdep < 0.126) ECalEdep = 0.;
        ECalEdep_array[ECal_EventID] = ECalEdep;

        for(Int_t itile = 0; itile < 36*num_total_layers; itile++)
        {
            HCal_tree->GetEntry(itile + i*36*num_total_layers);
            HCalEdep_array[HCal_EventID] += HCalTileEdep;
            if( (HCal_LayerID + 1) > num_total_layers - num_tail) 
            {
                TailCatch_array[HCal_EventID] += HCalTileEdep;
            }
        }
    }

    for(Int_t i = 0; i < num_events; i++)
    {
        Double_t fraction = 1.;
        if( (ECalEdep_array[i]/optimal_weight + HCalEdep_array[i]) != 0.) fraction = TailCatch_array[i] / (ECalEdep_array[i]/optimal_weight + HCalEdep_array[i]);
        Double_t reweighted_HCal = HCalEdep_array[i];
        if(fraction < 0.01)
        {
            h_TotalEdep->Fill(ECalEdep_array[i]/optimal_weight + HCalEdep_array[i]);
        }
    }

    return h_TotalEdep;
}
Double_t Resolution(std::string particle = "e-", Double_t energy = 1.0)
{   
    Bool_t EnableTailCatcher = kFALSE; // Tail Catcher used for hadrons
    Bool_t ECal_weight = kTRUE; // Weighting procedure used for hadrons
    if(particle == "e-" || particle == "pi0") ECal_weight = kFALSE;
    if(ECal_weight) std::cout<<"Using weighting for ECal."<<std::endl;
    else std::cout<<"Not using weighting for ECal."<<std::endl;

    TString file_name;
    //file_name.Form("test.root");
    file_name.Form("%s_QGSP/%s_%0.0fGeV.root", particle.c_str(), particle.c_str(), energy);
    // file_name.Form("%s_%0.0fGeV.root", particle.c_str(), energy);
    std::cout<<"Opening "<<file_name<<std::endl;
    TFile* data_file = new TFile(file_name);

    TH1D* h_TotalEdep = new TH1D("h_TotalEdep", "", 600, 0, Beam_MaxEnergy[energy]);
    TTree* Total_tree = (TTree*) data_file->Get("EdepTotal");
    TTree* HCal_tree = (TTree*) data_file->Get("HCalLayers");
    Int_t num_events = (Int_t) Total_tree->GetEntries();
    std::cout<<"Number of events: "<<num_events<<std::endl;

    Double_t optimal_weight = 1.0;
    Double_t HCal_reweight = 0.;
    if(ECal_weight) optimal_weight = ECalWeightingProcess(Total_tree, energy);
    //if(EnableTailCatcher && ECal_weight) HCal_reweight = TailCatcherProcess(Total_tree, HCal_tree, energy, optimal_weight);


    Double_t ECalEdep, HCalEdep;
    Total_tree->SetBranchAddress("ECal_Edep_Total", &ECalEdep);
    Total_tree->SetBranchAddress("HCal_Edep_Total", &HCalEdep);

    if(EnableTailCatcher) h_TotalEdep = (TH1D*) TailCatcher(Total_tree, HCal_tree, energy, optimal_weight);
    else
    {
        for(Int_t i = 0; i < num_events; i++)
        {
            Total_tree->GetEntry(i);
            if(ECalEdep < 0.183) ECalEdep = 0.;
            Double_t ECal_energy = ECalEdep/optimal_weight;
            Double_t total_energy;
            if(1==2) total_energy = ECal_energy + HCal_reweight;
            else total_energy = ECal_energy + HCalEdep;
            h_TotalEdep->Fill(total_energy);
        }
    }
    h_TotalEdep->SetTitle("");
    
    TCanvas* c_Resolution = new TCanvas("c_Resolution", "", 1000, 1000);
    h_TotalEdep->Draw();
    h_TotalEdep->GetXaxis()->SetTitle("Edep (MeV)");
    h_TotalEdep->GetYaxis()->SetTitle("Number of Events");

    TF1* f_gaus = new TF1("f_gaus", "gaus", 0, Beam_MaxEnergy[energy]);
    h_TotalEdep->Fit(f_gaus, "");

    Double_t mean = f_gaus->GetParameter(1);
    Double_t sigma = f_gaus->GetParameter(2);
    Double_t resolution;
    if(mean != 0) resolution = sigma/mean;
    else resolution = 0.;

    std::cout<<"Resolution is "<<resolution<<std::endl;

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
