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
#include "TRandom.h"
#include <map>
#include <vector>
#include <ROOT/TProcessExecutor.hxx>
#include <tuple>

std::map<Double_t, Double_t> Beam_MaxEnergy = {{1.0, 100}, {2.0, 100}, {3.0, 150}, {5.0, 200}, {10.0, 400}, {20.0, 700}, {30.0, 1000}, {40.0, 1500}, {50.0, 2000}, {60., 2500}, {70., 3000}, {80., 3500}, {90., 4000}, {100., 4500}};

double threshold = 6.1;

double reco(double E0)
{
  double E = E0;
  if(E0 < threshold) return 0;
  else{ 
    double a = 0.1;
    double b = 0.0015;
    double sigma =  E* TMath::Sqrt(a*a/E +b*b) ;
    double random = gRandom->Gaus(E*0.03,sigma); 
    E = random;
    return E;
  }
}

void Tuning(Int_t id = 0, Int_t smear = 0)
{
  std::string particle = "e-";
  std::string dd_particle = "electron";
  Double_t energy = 10.;

  TString file_name;
  Double_t max_energy = 0.;
  switch(id)
  {
    case 0:
      file_name.Form("%s_QGSP/%s_%0.0fGeV_20deg.root", particle.c_str(), particle.c_str(), energy);
      max_energy = Beam_MaxEnergy[energy];
      break;
    case 1:
      file_name.Form("../reconstruction_benchmarks/benchmarks/clustering/sim_endcap_%s.root", dd_particle.c_str());
      if(smear) max_energy = Beam_MaxEnergy[energy];
      else max_energy = energy*1e3;
      break;
    case 2:
      file_name.Form("../reconstruction_benchmarks/benchmarks/clustering/rec_endcap_%s.root", dd_particle.c_str());
      max_energy = Beam_MaxEnergy[energy]*3;
      break;
  }
  std::cout<<"Opening "<<file_name<<std::endl;
  TFile* data_file = new TFile(file_name);

  TH1D* h_TotalEdep = new TH1D("h_TotalEdep", "", 600, max_energy/5, max_energy);
  TTree* Total_tree = NULL;
  switch(id)
  {
    case 0:
      Total_tree = (TTree*) data_file->Get("EdepTotal");
      break;
    default:
      Total_tree = (TTree*) data_file->Get("events");
  }
  Int_t num_events = (Int_t) Total_tree->GetEntries();
  std::cout<<"Number of events: "<<num_events<<std::endl;

  Double_t ECalEdepD[1000], HCalEdepD[1000];
  Float_t ECalEdepF[1000], HCalEdepF[1000];
  switch(id)
  {
    case 0:
      Total_tree->SetBranchAddress("ECal_Edep_Total", &ECalEdepD[0]);
      Total_tree->SetBranchAddress("HCal_Edep_Total", &HCalEdepD[0]);
      break;
    case 1:
      Total_tree->SetBranchAddress("EcalEndcapPHits.energyDeposit", ECalEdepD);
      Total_tree->SetBranchAddress("HcalEndcapPHits.energyDeposit", HCalEdepD);
      break;
    case 2:
      Total_tree->SetBranchAddress("EcalEndcapPHitsReco.energy", ECalEdepF);
      Total_tree->SetBranchAddress("HcalEndcapPHitsReco.energy", HCalEdepF);
      break;
  }

  for(Int_t i = 0; i < num_events; i++)
  {
    Total_tree->GetEntry(i);
    Double_t ECal_energy = 0.;
    Double_t HCal_energy = 0.;
    switch(id)
    {
      case 0:
        ECal_energy = ECalEdepD[0];
        HCal_energy = HCalEdepD[0];
        break;
      case 1:
        for(Int_t j=0; j< Total_tree->GetLeaf("EcalEndcapPHits.energyDeposit")->GetLen(); j++)
          ECal_energy += ECalEdepD[j]*1e3;
        for(Int_t j=0; j< Total_tree->GetLeaf("HcalEndcapPHits.energyDeposit")->GetLen(); j++)
          HCal_energy += HCalEdepD[j]*1e3;
        if(smear) ECal_energy = reco(ECal_energy);
        break;
      case 2:
        for(Int_t j=0; j< Total_tree->GetLeaf("EcalEndcapPHitsReco.energy")->GetLen(); j++)
          ECal_energy += ECalEdepF[j]*1e3;
        for(Int_t j=0; j< Total_tree->GetLeaf("HcalEndcapPHitsReco.energy")->GetLen(); j++)
          HCal_energy += HCalEdepF[j]*1e3;
        break;
    }
    if(ECal_energy < 0.183) ECal_energy = 0.;
    Double_t total_energy = ECal_energy;  // + HCal_energy;
    h_TotalEdep->Fill(total_energy);
  }
  h_TotalEdep->SetTitle("");

  TCanvas* c_Resolution = new TCanvas("c_Resolution", "", 1000, 1000);
  h_TotalEdep->Draw();
  h_TotalEdep->GetXaxis()->SetTitle("Edep (MeV)");
  h_TotalEdep->GetYaxis()->SetTitle("Number of Events");

  TF1* f_gaus = new TF1("f_gaus", "gaus", max_energy/5, max_energy);
  h_TotalEdep->Fit(f_gaus, "");

  Double_t mean = f_gaus->GetParameter(1);
  Double_t sigma = f_gaus->GetParameter(2);
  Double_t resolution;
  if(mean != 0) resolution = sigma/mean;
  else resolution = 0.;

  std::cout<<"Resolution is "<<resolution<<std::endl;

  TString info_text;
  info_text = Form("%s at %0.0f GeV", particle.c_str(), energy);
  switch(id)
  {
    case 0:
      info_text += " in Geant4";
      break;
    case 1:
      info_text += " in DD4hep";
      if(smear)
        info_text += " w/ smearing";
      break;
    case 2:
      info_text += " in Juggler";
      break;
  }
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
}
